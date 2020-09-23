"""
Copyright 2020, University Corporation for Atmospheric Research
See LICENSE.txt for details
"""

import os

import numpy as np

from pyreshaper import iobackend

from .data import config


def check_outfile(infiles, prefix, tsvar, suffix, metadata, once, **kwds):
    """
    Check that a PyReshaper generated output file is correct
    """

    assertions = {}

    def _assert(key, value):
        assertions[key] = value

    outfile = '{0}{1}{2}'.format(prefix, tsvar, suffix)
    _assert('{0!r} exists'.format(outfile), os.path.exists(outfile))
    if not os.path.exists(outfile):
        return assertions
    ncout = iobackend.NCFile(outfile)

    if 'meta1d' in kwds and kwds['meta1d'] is True:
        metadata.append('time')

    if 'metafile' in kwds and kwds['metafile']:
        metafile = iobackend.NCFile('metafile.nc')
        _assert(
            '{0}: Extra time-invariant metadata found'.format(outfile),
            set(config.xtimvars).issubset(set(ncout.variables.keys())),
        )
        for v in config.xtimvars:
            _assert(
                '{0}: Extra time-invariant metadata dimensions'.format(outfile),
                ncout.variables[v].dimensions == ('lat', 'lon'),
            )
    else:
        metafile = None

    series_step = 0
    for infile in infiles:
        _assert('{0!r} exists'.format(infile), os.path.exists(infile))
        if not os.path.exists(infile):
            return assertions

        ncinp = iobackend.NCFile(infile)
        nsteps = ncinp.dimensions['time']
        if infile == infiles[0]:
            scvars = [v for v in ncinp.variables if ncinp.variables[v].dimensions == ()]
            tivars = [v for v in ncinp.variables if 'time' not in ncinp.variables[v].dimensions]
            tsvars = [
                v
                for v in ncinp.variables
                if 'time' in ncinp.variables[v].dimensions and v not in metadata
            ]
            if once:
                tsvars.append('once')

            outdims = {
                'lat': ncinp.dimensions['lat'],
                'lon': ncinp.dimensions['lon'],
                'strlen': ncinp.dimensions['strlen'],
            }

            outmeta = [v for v in ncinp.variables if v not in tsvars]

            _assert('{0}: variable {1!r} found in input'.format(outfile, tsvar), tsvar in tsvars)
            _assert(
                '{0}: global attribute names equal'.format(outfile), ncout.ncattrs == ncinp.ncattrs
            )
            for a in set(ncout.ncattrs).intersection(set(ncinp.ncattrs)):
                _assert(
                    '{0}: global attribute {1} values equal'.format(outfile, a),
                    ncout.getncattr(a) == ncinp.getncattr(a),
                )
            for d, v in outdims.items():
                _assert('{0}: {1!r} in dimensions'.format(outfile, d), d in ncout.dimensions)
                _assert('{0}: dimensions[{1!r}]'.format(outfile, d), ncout.dimensions[d] == v)
            _assert("{0}: 'time' in dimensions".format(outfile), 'time' in ncout.dimensions)
            _assert("{0}: 'time' unlimited".format(outfile), ncout.unlimited('time'))
            if once:
                all_vars = outmeta if tsvar == 'once' else [tsvar]
            else:
                all_vars = [tsvar] + outmeta
            if metafile:
                all_vars += config.xtimvars
            _assert(
                '{0}: variable names same'.format(outfile),
                set(ncout.variables.keys()) == set(all_vars),
            )
            for v in all_vars:
                if v in scvars:
                    expected = ()
                elif v in ncinp.dimensions:
                    expected = (v,)
                elif v in tivars + config.xtimvars:
                    expected = ('lat', 'lon')
                elif v in config.chvars:
                    expected = ('time', 'strlen')
                else:
                    expected = ('time', 'lat', 'lon')
                _assert(
                    '{0}: {1}.dimemsions equal'.format(outfile, v),
                    ncout.variables[v].dimensions == expected,
                )

        for v in all_vars:
            if v in config.xtimvars:
                expected = metafile.variables[v].get_value()
            else:
                expected = ncinp.variables[v].get_value()
            if v == 'time':
                oslice = slice(series_step, series_step + nsteps)
                actual = ncout.variables[v][oslice]
            elif 'time' in ncout.variables[v].dimensions:
                oslice = [slice(None)] * (2 if v in config.chvars else 3)
                oslice[0] = slice(series_step, series_step + nsteps)
                actual = ncout.variables[v][tuple(oslice)]
            else:
                actual = ncout.variables[v].get_value()
            _assert(('{0}: {1!r} values equal').format(outfile, v), np.all(actual == expected))

        series_step += nsteps
        ncinp.close()
    if metafile:
        metafile.close()
    ncout.close()

    return assertions


def check_var_in(var, fname):
    ncf = iobackend.NCFile(fname)
    value = var in ncf.variables
    ncf.close()
    return value
