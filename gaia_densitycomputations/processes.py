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
from gaia import types

import gaia.formats as formats
import gdal
import os
import osr
import ogr
import numpy as np
import itertools
import geopandas
import sys

from gaia.inputs import GaiaIO
from gaia.gaia_process import GaiaProcess
from gaia_densitycomputations import config
from gaia.geo.geo_inputs import RasterFileIO
from skimage.graph import route_through_array
import matplotlib.pyplot as plt
from math import sqrt, ceil


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
    required_args = [
        {
            'name': 'resolution',
            'title': 'Output image size',
            'description': 'Output image size (cols,rows), ex: "200,100"',
            'type': str
        }
    ]

    default_output = formats.RASTER

    def __init__(self, resolution, **kwargs):
        """
        Create an instance of the DensityComputationsProcess class
        :param resolution: grid cell resolution of the output
        :param kwargs: Optional keyword arguments
        """
        self.resolution = [int(n) for n in resolution.split(',')]
        super(SimpleGridDensityProcess, self).__init__(**kwargs)

        if not self.output:
            self.output = RasterFileIO(name='result', uri=self.get_outpath())

    def calculate_density(self):
        """
        Calculate a simple estimate of grid density by summing points
        within each grid cell.
        """

        shpDriver = ogr.GetDriverByName('GeoJSON')
        strjson = self.inputs[0].read().to_json()
        dataSource = shpDriver.Open(strjson, 0)

        # Open the source file, and exit if doesn't exist
        if dataSource is None:
            print 'Could not open file ' + self.inputs[0].uri
            sys.exit(1)

        if os.path.exists(self.output.uri):
            shpDriver.DeleteDataSource(self.output.uri)
        else:
            self.output.create_output_dir(self.output.uri)

        # Get the layer
        layer = dataSource.GetLayer()
        # Get the layer extent
        extent = layer.GetExtent()

        # Open the layer
        # Set the bounding box
        xmin = extent[0]
        ymin = extent[2]
        xmax = extent[1]
        ymax = extent[3]

        # Number of columns and rows
        nbrColumns = self.resolution[0]
        nbrRows = self.resolution[1]

        # Caculate the cell size in x and y direction
        csx = (xmax - xmin) / nbrColumns
        csy = (ymax - ymin) / nbrRows

        rows = []
        i = ymax
        while i > ymin:
            j = xmin
            cols = []
            while j < xmax:
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

        originX = extent[0]
        originY = extent[3]

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
        ct = gdal.ColorTable()
        ct.CreateColorRamp(0, (0, 0, 255), 14, (0, 255, 0))
        ct.CreateColorRamp(15, (0, 255, 0), 30, (127, 127, 0))
        ct.CreateColorRamp(30, (127, 127, 0), 50, (255, 0, 0))
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

PLUGIN_CLASS_EXPORTS = [
    SimpleGridDensityProcess,
]
