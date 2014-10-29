#!/usr/bin/env python

from distutils.core import setup

setup(name='PyReshaper',
      version='0.9.0',
      description='Python Time-Slice to Time-Series NetCDF Converter',
      author='Kevin Paul',
      author_email='kpaul@ucar.edu',
      url='https://wiki.ucar.edu/display/~kpaul/PyReshaper',
      packages=['pyreshaper'],
      scripts=['bin/slice2series'],
      requires=['Nio', 'mpi4py']
     )
