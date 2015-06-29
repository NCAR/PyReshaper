#!/usr/bin/env python
"""
PyReshaper -- Setup Script

Copyright 2015, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

from setuptools import setup

exec(open('source/pyreshaper/version.py').read())

setup(name='PyReshaper',
      version=__version__,
      description='Python Time-Slice to Time-Series NetCDF Converter',
      author='Kevin Paul',
      author_email='kpaul@ucar.edu',
      url='https://github.com/NCAR-CISL-ASAP/PyReshaper',
      download_url='https://github.com/NCAR-CISL-ASAP/PyReshaper/tarball/v' + __version__,
      license='https://github.com/NCAR-CISL-ASAP/PyReshaper/blob/master/LICENSE.txt',
      packages=['pyreshaper'],
      package_dir={'pyreshaper': 'source/pyreshaper'},
      package_data={'pyreshaper': ['LICENSE.txt']},
      scripts=['bin/slice2series'],
      install_requires=['Nio', 'mpi4py>=1.3', 'asaptools>=0.4']
      )
