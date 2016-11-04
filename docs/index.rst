Gaia Density Computations Plugin
================================

This is a plugin for Gaia (https://github.com/OpenDataAnalytics/gaia) that
calculates the point density of a dataset given a raster grid of a specified size (in columns and rows),
using one of the following processes:

    * `SimpleDensityProcess <gaia_densitycomputations.html#gaia_densitycomputations.processes.SimpleGridDensityProcess>`__

      * Calculate point density by adding the number of points within each grid cell.

      * Required input: a vector dataset containing points

      * Optional arguments:

        * width: Width of output image in pixels (default=1024)

      * `Example <gaia_processes.html#Simple-Grid-Density-Process>`__

    * `KernelDensityProcess <gaia.geo.html#gaia.geo.processes_vector.BufferProcess>`__

      * Calculates point density using a Gaussian kernel density function.

      * Required input: a vector dataset containing points

      * Optional arguments:

        * width: Width of output image in pixels (default=1024)

        * weight: Column of vector dataset, with integer values.  Points will be multiplied by this value.

        * bandwidth: The method used to calculate the kernel density estimator bandwidth.  Valid values are 'scott' (default), 'silverman', or a numerical value.

      * `Example <gaia_processes.html#Kernel-Density-Process>`__


Installation
-----------------

- git clone https://github.com/OpenDataAnalytics/gaia-densitycomputations-plugin.git
- cd gaia-densitycomputations-plugin
- pip install -e .
- pip install -r requirements


Testing
-----------------

- pip install -r requirements-dev.txt
- python -m unittest discover


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   gaia_densitycomputations


.. _Gaia: http://www.github.com/opendataanalytics/gaia
.. _Kitware: http://www.kitware.com
.. _Epidemico: http://epidemico.com
