#!/usr/bin/env python
"""
PyReshaper -- Setup Script


Copyright 2019 University Corporation for Atmospheric Research

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from setuptools import setup

__version__ = None
exec(open('source/pyreshaper/version.py').read())

setup(name='PyReshaper',
      version=__version__,
      description='Python Time-Slice to Time-Series NetCDF Converter',
      author='Kevin Paul',
      author_email='kpaul@ucar.edu',
      url='https://github.com/NCAR/PyReshaper',
      download_url='https://github.com/NCAR/PyReshaper/tarball/v' + __version__,
      license='https://github.com/NCAR/PyReshaper/blob/master/LICENSE.txt',
      packages=['pyreshaper'],
      package_dir={'pyreshaper': 'source/pyreshaper'},
      package_data={'pyreshaper': ['LICENSE.txt']},
      scripts=['scripts/s2smake', 'scripts/s2srun'],
      install_requires=['mpi4py', 'asaptools']
      )
