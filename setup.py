#!/usr/bin/env python

from distutils.core import setup
import shutil

# Copy license file into package directory
shutil.copy('LICENSE.txt', 'source/pyreshaper')

setup(name='PyReshaper',
      version='0.9.2',
      description='Python Time-Slice to Time-Series NetCDF Converter',
      author='Kevin Paul',
      author_email='kpaul@ucar.edu',
      url='https://github.com/NCAR-CISL-ASAP/PyReshaper',
      download_url='https://github.com/NCAR-CISL-ASAP/PyReshaper/tarball/v0.9.2',
      license='https://github.com/NCAR-CISL-ASAP/PyReshaper/blob/master/LICENSE.txt',
      packages=['pyreshaper'],
      package_dir={'pyreshaper': 'source/pyreshaper'},
      package_data={'pyreshaper': ['LICENSE.txt']},
      scripts=['bin/slice2series'],
      requires=['Nio', 'mpi4py']
     )
