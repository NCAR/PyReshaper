"""
The PyReshaper

A tool for converting NetCDF time-slice files into time-series format

:AUTHORS: Kevin Paul, John Dennis, Sheri Mickelson, Haiying Xu
:COPYRIGHT: 2015, University Corporation for Atmospheric Research
:LICENSE: See the LICENSE.rst file for details

Send questions and comments to Kevin Paul (kpaul@ucar.edu).
"""

from __future__ import absolute_import

from .version import __version__

from . import reshaper
from . import specification
