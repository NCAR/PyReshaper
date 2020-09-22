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

__version__ = '1.0.7'

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(name='PyReshaper',
      version=__version__,
      author='Kevin Paul',
      author_email='kpaul@ucar.edu',
      description='Python Time-Slice to Time-Series NetCDF Converter',
      long_description=long_description,
      long_description_content_type="text/x-rst",
      url='https://github.com/NCAR/PyReshaper',
      packages=['pyreshaper'],
      package_dir={'pyreshaper': 'source/pyreshaper'},
      classifiers=[
          "Programming Language :: Python :: 2.7",
          "License :: OSI Approved :: Apache Software License"
          "Operating System :: OS Independent",
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Atmospheric Science",
          "Topic :: Utilities"
      ],
      python_requires='>=2.7,<3.0',
      entry_points="""
          [console_scripts]
          s2smake=pyreshaper.cli.s2smake:main
          s2srun=pyreshaper.cli.s2srun:main
          """,
      install_requires=['mpi4py', 'asaptools']
      )
