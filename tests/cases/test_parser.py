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
import gdal
import unittest
import numpy as np
from gaia import formats
from gaia.parser import deserialize
from gaia.core import config

testfile_path = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '../data')


class TestDensityComputationsViaParser(unittest.TestCase):
    """Tests for the Gaia Least Cost plugin via Parser"""

    def test_process_densitycomputations(self):
        """Test Density Computations Process"""
        with open(os.path.join(testfile_path,
                               'densitycomputations.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            expected_layer = process.output.read()
            # Get layer stats
            expected_results = \
                expected_layer.GetRasterBand(1).GetStatistics(0, 1)

            actual_layer = gdal.Open(os.path.join(
                testfile_path,
                'densitycomputations_process_results.tif'), gdal.GA_Update)
            actual_results = actual_layer.GetRasterBand(1).GetStatistics(0, 1)

            expected_results_rounded = np.around(expected_results, decimals=2)
            actual_results_rounded = np.around(actual_results, decimals=2)
            self.assertEquals(np.all(expected_results_rounded),
                              np.all(actual_results_rounded))
        finally:
            if process:
                process.purge()
