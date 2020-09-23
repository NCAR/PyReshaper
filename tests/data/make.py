#!/usr/bin/env python
"""
Copyright 2020, University Corporation for Atmospheric Research
See LICENSE.txt for details
"""

import numpy as np

from pyreshaper import iobackend

from . import config


def generate_data(backend='netCDF4'):
    """
    Generate dataset for testing purposes
    """
    iobackend.set_backend(backend)

    # Test Data Generation
    for i in range(len(config.slices) + 1):

        # Open the file for writing
        fname = config.slices[i] if i < len(config.slices) else 'metafile.nc'
        fobj = iobackend.NCFile(fname, mode='w')

        # Write attributes to file
        for name, value in config.fattrs.items():
            fobj.setncattr(name, value)

        # Create the dimensions in the file
        fobj.create_dimension('lat', config.nlat)
        fobj.create_dimension('lon', config.nlon)
        fobj.create_dimension('time', None)
        fobj.create_dimension('strlen', config.nchar)

        # Create the coordinate variables & add attributes
        lat = fobj.create_variable('lat', 'f', ('lat',))
        lon = fobj.create_variable('lon', 'f', ('lon',))
        time = fobj.create_variable('time', 'f', ('time',))

        # Set the coordinate variable attributes
        lat.setncattr('long_name', 'latitude')
        lat.setncattr('units', 'degrees_north')

        lon.setncattr('long_name', 'longitude')
        lon.setncattr('units', 'degrees_east')

        time.setncattr('long_name', 'time')
        time.setncattr('units', 'days since 01-01-0001')
        time.setncattr('calendar', 'noleap')

        # Set the values of the coordinate variables
        lat[:] = np.linspace(-90, 90, config.nlat, dtype=np.float32)
        lon[:] = np.linspace(-180, 180, config.nlon, endpoint=False, dtype=np.float32)
        time[:] = np.arange(i * config.ntime, (i + 1) * config.ntime, dtype=np.float32)

        # Create the scalar variables
        for n in range(len(config.scalars)):
            vname = config.scalars[n]
            v = fobj.create_variable(vname, 'd', tuple())
            v.setncattr('long_name', 'scalar{0}'.format(n))
            v.setncattr('units', '[{0}]'.format(vname))
            v.assign_value(np.float64(n * 10))

        # Create the time-invariant metadata variables
        all_timvars = config.timvars + ([] if i < len(config.slices) else config.xtimvars)
        for n in range(len(all_timvars)):
            vname = all_timvars[n]
            v = fobj.create_variable(vname, 'd', ('lat', 'lon'))
            v.setncattr('long_name', 'time-invariant metadata {0}'.format(n))
            v.setncattr('units', '[{0}]'.format(vname))
            v[:] = np.ones((config.nlat, config.nlon), dtype=np.float64) * n

        # Create the time-variant character variables
        for n in range(len(config.chvars)):
            vname = config.chvars[n]
            v = fobj.create_variable(vname, 'c', ('time', 'strlen'))
            v.setncattr('long_name', 'character array {0}'.format(n))
            vdata = [str((n + 1) * m) * (m + 1) for m in range(config.ntime)]
            v[:] = (
                np.array(vdata, dtype='S{}'.format(config.nchar))
                .view('S1')
                .reshape(config.ntime, config.nchar)
            )

        # Create the time-variant metadata variables
        for n in range(len(config.tvmvars)):
            vname = config.tvmvars[n]
            v = fobj.create_variable(vname, 'd', ('time', 'lat', 'lon'))
            v.setncattr('long_name', 'time-variant metadata {0}'.format(n))
            v.setncattr('units', '[{0}]'.format(vname))
            v[:] = np.ones((config.ntime, config.nlat, config.nlon), dtype=np.float64) * n

        # Create the time-series variables
        for n in range(len(config.tsvars)):
            vname = config.tsvars[n]
            v = fobj.create_variable(vname, 'd', ('time', 'lat', 'lon'), fill_value=1e36)
            v.setncattr('long_name', 'time-series variable {0}'.format(n))
            v.setncattr('units', '[{0}]'.format(vname))
            v.setncattr('missing_value', 1e36)
            vdata = np.ones((config.ntime, config.nlat, config.nlon), dtype=np.float64) * n
            vmask = np.random.choice(
                [True, False], config.ntime * config.nlat * config.nlon
            ).reshape(config.ntime, config.nlat, config.nlon)
            v[:] = np.ma.MaskedArray(vdata, mask=vmask)


if __name__ == '__main__':
    generate_data(backend='netCDF4')
