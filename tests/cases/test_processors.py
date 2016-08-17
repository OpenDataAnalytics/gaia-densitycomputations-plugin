#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
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
###############################################################################
import json
import os
import unittest
import pysal
import gdal
from gaia import formats
from gaia.geo.geo_inputs import RasterFileIO
from gaia_densitycomputations.processes import DensityComputationsProcess

testfile_path = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '../data')


class TestDensityComputationsProcessors(unittest.TestCase):

    def test_process_densitycomputations(self):
        """
        Test DensityComputationsProcess for raster inputs
        """

        uri = os.path.join(testfile_path, 'ports_and_harbours.geojson')
        resolution = {
            'nCol': 200,
            'nRow': 100
        }
        pixelWidth = 10
        pixelHeight = 10
        outputWidth = 1000

        process = DensityComputationsProcess(inputs=[{ "uri": uri, "resolution": resolution, pixelWidth: 10, pixelHeight: 10, outputWidth: 1000 }])
        try:
            process.compute()
            expected_layer = process.output.read()
            # Get layer stats
            expected_results = expected_layer.GetRasterBand(1).GetStatistics(0, 1)

            actual_layer = gdal.Open(os.path.join(
                                    testfile_path,
                                    'densitycomputations_process_results.tif'), gdal.GA_Update)
            actual_results = actual_layer.GetRasterBand(1).GetStatistics(0, 1)

            self.assertEquals(expected_results,
                              actual_results)
        finally:
            if process:
                process.purge()
