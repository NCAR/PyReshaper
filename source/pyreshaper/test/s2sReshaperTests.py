"""
Parallel Tests for the Reshaper class

Copyright 2015, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

import unittest

import sys
from glob import glob
from cStringIO import StringIO
from os import linesep as eol
from os import remove
from os.path import exists
from mpi4py import MPI

import Nio
import numpy as np

from pyreshaper.reshaper import Slice2SeriesReshaper, create_reshaper
from pyreshaper.specification import Specifier

import mkTestData

MPI_COMM_WORLD = MPI.COMM_WORLD


class S2SReshaperTests(unittest.TestCase):

    def setUp(self):

        # Parallel Management - Just for Tests
        self.rank = MPI_COMM_WORLD.Get_rank()
        self.size = MPI_COMM_WORLD.Get_size()

        # Test Data Generation
        self._clean_directory()
        if self.rank == 0:
            mkTestData.generate_data()
        MPI_COMM_WORLD.Barrier()

    def tearDown(self):
        self._clean_directory()

    def _clean_directory(self):
        if self.rank == 0:
            for ncfile in glob('*.nc'):
                remove(ncfile)
        MPI_COMM_WORLD.Barrier()

    def _test_header(self, testname):
        if self.rank == 0:
            hline = '-' * 70
            print hline
            print testname
            print hline

    def _convert_header(self, infiles, prefix, suffix, metadata,
                        ncfmt, clevel, serial, verbosity, wmode, once,
                        print_diags=False):
        nfiles = len(infiles)
        ncvers = '3' if ncfmt == 'netcdf' else ('4c' if ncfmt == 'netcdf4c'
                                                else '4')
        self._test_header(("convert() - {} infile(s), NC{}-CL{}, serial={},{}"
                           "            verbosity={}, wmode={!r}, once={}"
                           "").format(nfiles, ncvers, clevel, serial, eol,
                                      verbosity, wmode, once))

    def _assertion(self, name, actual, expected,
                   data=None, show=True, assertion=None):
        rknm = '[{}/{}] {}'.format(self.rank, self.size, name)
        spcr = ' ' * len(rknm)
        msg = eol + rknm
        if data:
            msg += ' - Input:    {}'.format(data) + eol + spcr
        msg += ' - Actual:   {}'.format(actual) + eol + spcr
        msg += ' - Expected: {}'.format(expected)
        if show:
            print msg
        if assertion:
            assertion(actual, expected, msg)
        else:
            self.assertEqual(actual, expected, msg)

    def _check_outfile(self, infiles, prefix, tsvar, suffix,
                       metadata, once, **kwds):
        if self.rank == 0:

            outfile = '{}{}{}'.format(prefix, tsvar, suffix)
            ncout = Nio.open_file(outfile, 'r')

            series_step = 0
            for infile in infiles:
                ncinp = Nio.open_file(infile, 'r')
                nsteps = ncinp.dimensions['time']
                if infile == infiles[0]:
                    scalars = [v for v in ncinp.variables
                               if ncinp.variables[v].dimensions == ()]
                    tivars = [v for v in ncinp.variables
                              if 'time' not in ncinp.variables[v].dimensions]
                    tsvars = [v for v in ncinp.variables
                              if 'time' in ncinp.variables[v].dimensions and
                              v not in metadata]
                    if once:
                        tsvars.append('once')

                    outdims = {'lat': ncinp.dimensions['lat'],
                               'lon': ncinp.dimensions['lon']}

                    outmeta = [v for v in ncinp.variables if v not in tsvars]

                    self._assertion(("{}: variable "
                                     "{!r} exists").format(outfile, tsvar),
                                    tsvar in tsvars, True)
                    self._assertion("exists({!r})".format(outfile),
                                    exists(outfile), True)
                    self._assertion("{}: attributes equal".format(outfile),
                                    ncout.attributes, ncinp.attributes,
                                    assertion=self.assertDictEqual)
                    for d, v in outdims.iteritems():
                        self._assertion("{}: {!r} in dimensions".format(outfile, d),
                                        d in ncout.dimensions, True)
                        self._assertion("{}: dimensions[{!r}]".format(outfile, d),
                                        ncout.dimensions[d], v)
                    self._assertion("{}: 'time' in dimensions".format(outfile),
                                    'time' in ncout.dimensions, True)
                    self._assertion("{}: 'time' unlimited".format(outfile),
                                    ncout.unlimited('time'), True)
                    if once:
                        all_vars = outmeta if tsvar == 'once' else [tsvar]
                    else:
                        all_vars = [tsvar] + outmeta
                    self._assertion("{}: variable set".format(outfile),
                                    set(ncout.variables.keys()), set(all_vars),
                                    assertion=self.assertSetEqual)
                    for v in all_vars:
                        if v in scalars:
                            expected = ()
                        elif v in ncinp.dimensions:
                            expected = (v,)
                        elif v in tivars:
                            expected = ('lat', 'lon')
                        else:
                            expected = ('time', 'lat', 'lon')
                        self._assertion("{}: {}.dimemsions equal".format(outfile, v),
                                        ncout.variables[v].dimensions, expected)

                for v in all_vars:
                    expected = ncinp.variables[v].get_value()
                    if v == 'time':
                        oslice = slice(series_step, series_step + nsteps)
                        actual = ncout.variables[v][oslice]
                    elif 'time' in ncout.variables[v].dimensions:
                        oslice = [slice(None)] * 3
                        oslice[0] = slice(series_step, series_step + nsteps)
                        actual = ncout.variables[v][tuple(oslice)]
                    else:
                        actual = ncout.variables[v].get_value()
                    self._assertion(("{}: {!r} values "
                                     "equal").format(outfile, v),
                                    np.all(actual == expected), True)

                series_step += nsteps
                ncinp.close()
            ncout.close()
        MPI_COMM_WORLD.Barrier()

    def _run_convert(self, infiles, prefix, suffix, metadata,
                     ncfmt, clevel, serial, verbosity, wmode, once,
                     print_diags=False):
        if not (serial and self.rank > 0):
            spec = Specifier(infiles=infiles, ncfmt=ncfmt, compression=clevel,
                             prefix=prefix, suffix=suffix, metadata=metadata)
            rshpr = create_reshaper(spec, serial=serial,
                                    verbosity=verbosity,
                                    wmode=wmode, once=once)
            rshpr.convert()
            if print_diags:
                rshpr.print_diagnostics()
        MPI_COMM_WORLD.Barrier()

    def _run_convert_assert_no_output(self, infiles, prefix, suffix, metadata,
                                      ncfmt, clevel, serial, verbosity, wmode,
                                      once, print_diags=False):
        oldout = sys.stdout
        newout = StringIO()
        sys.stdout = newout
        self._run_convert(infiles, prefix, suffix, metadata, ncfmt, clevel,
                          serial, 0, wmode, once, print_diags=False)
        actual = newout.getvalue()
        self._assertion("stdout empty", actual, '')
        sys.stdout = oldout

    def _test_create_reshaper(self, serial, verbosity, wmode):
        self._test_header(("create_reshaper(serial={}, verbosity={}, "
                           "wmode={!r})").format(serial, verbosity, wmode))
        if not (serial and self.rank > 0):
            spec = Specifier(infiles=mkTestData.slices, ncfmt='netcdf',
                             compression=0, prefix='output.', suffix='.nc',
                             metadata=[])
            rshpr = create_reshaper(spec, serial=serial, verbosity=verbosity,
                                    wmode=wmode)
            self._assertion("type(reshaper)", type(rshpr),
                            Slice2SeriesReshaper)

    def test_create_reshaper_serial_V0_W(self):
        self._test_create_reshaper(serial=True, verbosity=0, wmode='w')

    def test_create_reshaper_serial_V1_W(self):
        self._test_create_reshaper(serial=True, verbosity=1, wmode='w')

    def test_create_reshaper_serial_V2_W(self):
        self._test_create_reshaper(serial=True, verbosity=2, wmode='w')

    def test_create_reshaper_serial_V1_O(self):
        self._test_create_reshaper(serial=True, verbosity=1, wmode='o')

    def test_create_reshaper_serial_V1_S(self):
        self._test_create_reshaper(serial=True, verbosity=1, wmode='s')

    def test_create_reshaper_serial_V1_A(self):
        self._test_create_reshaper(serial=True, verbosity=1, wmode='a')

    def test_create_reshaper_parallel_V1_W(self):
        self._test_create_reshaper(serial=False, verbosity=1, wmode='w')

    def test_convert_All_NC3_CL0_SER_V0_W(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': True, 'verbosity': 0, 'wmode': 'w', 'once': False,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert_assert_no_output(**args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_1_NC3_CL0_SER_V0_W(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices[0:1], 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': True, 'verbosity': 0, 'wmode': 'w', 'once': False,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert_assert_no_output(**args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC4_CL1_SER_V0_W(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf4', 'clevel': 1,
                'serial': True, 'verbosity': 0, 'wmode': 'w', 'once': False,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert_assert_no_output(**args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V1_W(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'wmode': 'w', 'once': False,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert(**args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V1_W_ONCE(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'wmode': 'w', 'once': True,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert(**args)
        self._check_outfile(tsvar='once', **args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V1_O(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'wmode': 'o', 'once': False,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert(**args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V1_O_ONCE(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'wmode': 'o', 'once': True,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert(**args)
        self._check_outfile(tsvar='once', **args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V1_S(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'wmode': 's', 'once': False,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert(**args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V1_S_ONCE(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'infiles': mkTestData.slices, 'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'wmode': 's', 'once': True,
                'print_diags': False}
        self._convert_header(**args)
        self._run_convert(**args)
        self._check_outfile(tsvar='once', **args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V3_A(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'once': False,
                'print_diags': False}
        self._convert_header(infiles=mkTestData.slices, wmode='a', **args)
        self._run_convert(infiles=[mkTestData.slices[0]], wmode='w', **args)
        for infile in mkTestData.slices[1:]:
            self._run_convert(infiles=[infile], wmode='a', **args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(infiles=mkTestData.slices, tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V3_A_ONCE(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'once': True,
                'print_diags': False}
        self._convert_header(infiles=mkTestData.slices, wmode='a', **args)
        self._run_convert(infiles=[mkTestData.slices[0]], wmode='w', **args)
        for infile in mkTestData.slices[1:]:
            self._run_convert(infiles=[infile], wmode='a', **args)
        self._check_outfile(infiles=mkTestData.slices, tsvar='once', **args)
        for tsvar in mkTestData.tsvars:
            self._check_outfile(infiles=mkTestData.slices, tsvar=tsvar, **args)

    def test_convert_All_NC3_CL0_PAR_V3_A_MISSING(self):
        mdata = [v for v in mkTestData.tvmvars]
        mdata.append('time')
        args = {'prefix': 'out.', 'suffix': '.nc',
                'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                'serial': False, 'verbosity': 1, 'once': False,
                'print_diags': False}
        missingvar = mkTestData.tsvars[2]
        self._convert_header(infiles=mkTestData.slices, wmode='a', **args)
        self._run_convert(infiles=[mkTestData.slices[0]], wmode='w', **args)
        self._run_convert(infiles=[mkTestData.slices[1]], wmode='a', **args)
        remove(args['prefix'] + missingvar + args['suffix'])
        for infile in mkTestData.slices[2:]:
            self._run_convert(infiles=[infile], wmode='a', **args)
        for tsvar in mkTestData.tsvars:
            if tsvar == missingvar:
                self._check_outfile(infiles=mkTestData.slices[2:], tsvar=tsvar, **args)
            else:
                self._check_outfile(infiles=mkTestData.slices, tsvar=tsvar, **args)


if __name__ == "__main__":
    hline = '=' * 70
    if MPI_COMM_WORLD.Get_rank() == 0:
        print hline
        print 'STANDARD OUTPUT FROM ALL TESTS:'
        print hline
    MPI_COMM_WORLD.Barrier()

    mystream = StringIO()
    tests = unittest.TestLoader().loadTestsFromTestCase(S2SReshaperTests)
    unittest.TextTestRunner(stream=mystream).run(tests)
    MPI_COMM_WORLD.Barrier()

    results = MPI_COMM_WORLD.gather(mystream.getvalue())
    if MPI_COMM_WORLD.Get_rank() == 0:
        for rank, result in enumerate(results):
            print hline
            print 'TESTS RESULTS FOR RANK ' + str(rank) + ':'
            print hline
            print str(result)
