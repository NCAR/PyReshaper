#!/usr/bin/env python

from distutils.core import setup

setup(name='PyReshaper',
      version='0.9.1',
      description='Python Time-Slice to Time-Series NetCDF Converter',
      author='Kevin Paul',
      author_email='kpaul@ucar.edu',
      url='https://github.com/NCAR-CISL-ASAP/PyReshaper',
      packages=['pyreshaper'],
      package_data={'pyreshaper': ['LICENSE.txt']},
      scripts=['bin/slice2series'],
      requires=['Nio', 'mpi4py']
     )
