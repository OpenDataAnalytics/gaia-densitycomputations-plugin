## gaia-density-computations-plugin

This is a plugin for Gaia (https://github.com/OpenDataAnalytics/gaia) that provides
a simple estimate of point density by creating a grid of square cells of a specified size
and then summing the points that lie within each grid.  The result is a raster dataset
with cell values equaling the number of points within that cell.

#### Documentation

Documentation for this plugin can be found at http://gaia-densitycomputations-plugin.readthedocs.org.

Documentation for this Gaia can be found at http://gaia.readthedocs.org.

#### Installation

  - pip install -e .
  - pip install -r requirements.txt

#### Inputs
  - A vector dataset of points

#### License

Copyright 2015 Kitware Inc.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0


Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
