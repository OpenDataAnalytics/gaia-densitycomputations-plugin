Gaia Least Cost Path Plugin
================================

This is a plugin for Gaia (https://github.com/OpenDataAnalytics/gaia) that
calculates the point density of a dataset given a raster grid of a specified size (in columns and rows).
The current implementation performs a simple calculation of adding the number of points within each grid cell.


An example of how to use this plugin can be found `here <gaia_processes.html>`__.

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
