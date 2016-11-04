#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
#  Copyright Kitware Inc. and Epidemico Inc.
#
#  Licensed under the Apache License, Version 2.0 ( the "License" );
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##############################################################################

import gaia.formats as formats
import gdal
import os
import osr
import ogr
import numpy as np
import sys
from scipy.stats import gaussian_kde
from gaia import types
from gaia.core import GaiaException
from gaia.gaia_process import GaiaProcess
from gaia.geo.geo_inputs import RasterFileIO, VectorFileIO


class SimpleGridDensityProcess(GaiaProcess):
    """
    Process for calculating the number of points within each grid cell
    of a specified size.
    """

    #: Tuple of required inputs; name, type , max # of each; None = no max
    required_inputs = [
        {'description': 'Feature dataset',
         'type': types.VECTOR,
         'max': 1
         }
    ]

    #: Required arguments, data types as dict
    optional_args = [
        {
            'name': 'width',
            'title': 'Output image width',
            'description': 'Output image width in pixels',
            'type': int
        }
    ]

    default_output = formats.RASTER
    width = 1024

    def __init__(self, **kwargs):
        """
        Create an instance of the DensityComputationsProcess class
        :param resolution: grid cell resolution of the output
        :param kwargs: Optional keyword arguments
        """
        super(SimpleGridDensityProcess, self).__init__(**kwargs)

        if not self.output:
            self.output = RasterFileIO(name='result', uri=self.get_outpath())

    def calculate_density(self):
        """
        Calculate a simple estimate of grid density by summing points
        within each grid cell.
        """

        vec_driver = ogr.GetDriverByName('GeoJSON')
        input = self.inputs[0]
        strjson = input.read(format=formats.JSON)
        crs = input.get_epsg()
        dataSource = vec_driver.Open(strjson, 0)

        # Open the source file, and exit if doesn't exist
        if dataSource is None:
            print 'Could not open file ' + self.inputs[0].uri
            sys.exit(1)

        if os.path.exists(self.output.uri):
            vec_driver.DeleteDataSource(self.output.uri)
        else:
            self.output.create_output_dir(self.output.uri)

        # Get the layer
        layer = dataSource.GetLayer()
        gb = GridBounds(layer, width=self.width, crs=crs)

        # Caculate the cell size in x and y direction
        csx = (gb.xmax - gb.xmin) / self.width
        csy = (gb.ymax - gb.ymin) / gb.height

        rows = []
        i = gb.ymax
        while i > gb.ymin:
            j = gb.xmin
            cols = []
            while j < gb.xmax:
                # Set a spatial filter
                layer.SetSpatialFilterRect(j, (i-csy), (j+csx), i)
                # And count the features
                cols.append(layer.GetFeatureCount())
                j += csx
            rows.append(cols)
            i -= csy

        array = np.array(rows)
        ncols = array.shape[1]
        nrows = array.shape[0]

        originX = gb.xmin
        originY = gb.ymin

        # Convert the results to geoTiff raster
        driver = gdal.GetDriverByName('GTiff')
        outRaster = driver.Create(
            self.output.uri, ncols, nrows, 1, gdal.GDT_Byte
        )

        # Set layer geo transform
        outRaster.SetGeoTransform((originX, csx, 0, originY, 0, -csy))
        outband = outRaster.GetRasterBand(1)
        # Add colors to the raster image
        outband.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
        ct = create_color_table(array.max())
        outband.SetColorTable(ct)
        outband.SetNoDataValue(0)
        outband.FlushCache()

        # Write the layer
        outband.WriteArray(array)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromEPSG(4326)
        # Set layer projection
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()
        outband = None

    def compute(self):
        """
        Run the process
        """
        self.calculate_density()


class KernelDensityProcess(GaiaProcess):

    #: Tuple of required inputs; name, type , max # of each; None = no max
    required_inputs = [
        {'description': 'Point dataset',
         'type': types.VECTOR,
         'max': 1
         }
    ]

    optional_args = [
        {
            'description':
                'Weight attribute (attribute should have integer values)',
            'name': 'weight',
            'title': 'Weight',
            'type': str
        },
        {
            'description':
                'Bandwidth estimator (scott, silverman, or a numeric value)',
            'name': 'bandwidth',
            'title': 'Bandwidth Estimator',
            'type': str
        },
        {
            'description':
                'Width of output in pixels',
            'name': 'resolution',
            'title': 'Resolutionr',
            'type': int
        }
    ]

    bandwidth = 'scott'
    width = 1024

    default_output = formats.RASTER

    def __init__(self, **kwargs):
        """
        Create an instance of the KernelDensityProcess class
        """
        super(KernelDensityProcess, self).__init__(**kwargs)

        try:
            self.bandwidth = float(self.bandwidth)
        except ValueError:
            if self.bandwidth not in ('scott', 'silverman'):
                raise GaiaException('Invalid bandwidth')

        if not self.output:
            self.output = RasterFileIO(name='result', uri=self.get_outpath())

    def compute(self):
        vec_driver = ogr.GetDriverByName('GeoJSON')
        input = self.inputs[0]
        strjson = input.read(format=formats.JSON)
        crs = input.get_epsg()
        data_source = vec_driver.Open(strjson, 0)
        layer = data_source.GetLayer()

        # Buffer result to avoid edge effects
        points = []
        for i in range(layer.GetFeatureCount()):
            feature = layer.GetNextFeature()
            geom = feature.GetGeometryRef()
            weight = (1 if not hasattr(self, 'weight') else
                      feature.GetField(self.weight))
            for j in range(weight):
                points.append((geom.GetX(), geom.GetY(), weight))

        xtrain = np.vstack([[c[1] for c in points],
                            [c[0] for c in points]]).T

        # Set up the data grid for the contour plot
        gb = GridBounds(layer, self.width, crs)

        xgrid, ygrid = construct_grids(gb.xres, gb.yres,
                                       self.width, gb.height,
                                       gb.xmin, gb.ymin)
        x, y = np.meshgrid(xgrid, ygrid)
        xy = np.vstack([y.ravel(), x.ravel()]).T

        kde = gaussian_kde(xtrain.T, bw_method=self.bandwidth)
        density = kde(xy.T)
        density = density.reshape(x.shape)
        density = (density/(density.max()/255.0)).astype(np.int32)

        self.output.create_output_dir(self.output.uri)
        driver = gdal.GetDriverByName('GTiff')
        outDs = driver.Create(self.output.uri, density.shape[1],
                              density.shape[0], 1, gdal.GDT_UInt16)
        outband = outDs.GetRasterBand(1)
        # write the data
        outband.WriteArray(density[::-1], 0, 0)
        ct = create_color_table(int(density.max()))
        outband.SetColorTable(ct)
        outband.SetNoDataValue(0)
        # # flush data to disk, set the NoData value and calculate stats
        outband.FlushCache()

        gt = [gb.xmin, gb.xres, 0, gb.ymax, 0, -gb.yres]
        outDs.SetGeoTransform(gt)
        outDs.SetProjection(gb.projection.ExportToWkt())
        outDs.FlushCache()
        del outDs


class GridBounds(object):
    """
    Helper class to determine bounds and pixel height of a raster image,
    based on the extent of a vector layer, specified pixel width, and
    projection CRS.
    """

    xmin = -180.0
    xmax = 180.0
    ymin = -90.0
    ymax = 90.0
    xres = 1024
    yres = 1024
    projection = None

    def __init__(self, layer, width=1024, crs=4326):
        """
        Initialize the object based on input parameters
        :param layer: Vector layer
        :param width: Desired width of output raster in pixels
        :param crs: EPSG code of the layer
        """
        self.projection = osr.SpatialReference()
        self.projection.ImportFromEPSG(int(crs))
        unit = self.projection.GetAttrValue('Unit')
        maxextent = [-180.0, 180.0, -90.0, 90.0]

        if unit != 'degree':
            wg84source = osr.SpatialReference()
            wg84source.ImportFromEPSG(4326)

            epsg = self.projection.GetAttrValue("AUTHORITY", 1)
            if epsg in ('900913', '3857', '102100', '102113', '54004', '41001'):
                maxextent = [-20037508.34, 20037508.34,
                             -20037508.34, 20037508.34]
            else:
                transform = osr.CoordinateTransformation(wg84source,
                                                         self.projection)
                ul = ogr.CreateGeometryFromWkt("POINT (-180.0 90.0)")
                lr = ogr.CreateGeometryFromWkt("POINT (180.0 -90.0)")
                for point in (ul, lr):
                    point.Transform(transform)
                maxextent = [ul.GetX(), lr.GetX(), ul.GetY(), lr.GetY()]

        extent = layer.GetExtent()
        x_buffer = (extent[1]-extent[0])/2
        y_buffer = (extent[3]-extent[2])/2

        # Set the bounding box
        self.xmin = max(extent[0]-x_buffer, maxextent[0])
        self.ymin = max(extent[2]-y_buffer, maxextent[2])
        self.xmax = min(extent[1]+x_buffer, maxextent[1])
        self.ymax = min(extent[3]+y_buffer, maxextent[3])

        self.xres = (self.xmax-self.xmin)/float(width)
        self.height = int(round((self.ymax-self.ymin)/self.xres))
        self.yres = (self.ymax-self.ymin)/float(self.height)


def create_color_table(max, min=0, intervals=None):
    """
    Create and return a color table.
    :param max_value: Maximum value of raster grid
    :param min_value: Optional minimum value (default=0)
    :param intervals: Optional list of ColorRamp values as lists
    :return: gdal.ColorTable
    """
    ct = gdal.ColorTable()
    if not intervals:
        intvl = int(max/4)
        ct.CreateColorRamp(min, (0, 0, 255), intvl-1, (0, 255, 0))
        ct.CreateColorRamp(intvl, (0, 255, 0), intvl*2-1, (127, 127, 0))
        ct.CreateColorRamp(intvl*2, (127, 127, 0), max, (255, 0, 0))
    else:
        for i in intervals:
            ct.CreateColorRamp(i[0], i[1], i[2], i[3])
    return ct


def construct_grids(xres, yres, columns, rows, minx, miny):
    """
    Create a grid of coordinates
    :param xres: Resolution of longitude columns
    :param yres: Resolution of latitude rows
    :param columns: Number of columns
    :param rows: Number of rows
    :param minx: Minimum X coordinate
    :param miny: Minimum Y coordinate
    :return: Grids of X and Y coordinate values
    """
    # x,y coordinates for corner cells
    xmin = minx + xres
    xmax = xmin + (columns * xres)
    ymin = miny + yres
    ymax = ymin + (rows * yres)

    # x coordinates of the grid cells
    xgrid = np.arange(xmin, xmax, xres)
    # y coordinates of the grid cells
    ygrid = np.arange(ymin, ymax, yres)

    return xgrid, ygrid


PLUGIN_CLASS_EXPORTS = [
    SimpleGridDensityProcess,
    KernelDensityProcess
]
