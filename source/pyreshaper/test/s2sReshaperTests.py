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

MPI_COMM_WORLD = MPI.COMM_WORLD


class S2SReshaperTests(unittest.TestCase):

    def setUp(self):

        # Parallel Management - Just for Tests
        self.rank = MPI_COMM_WORLD.Get_rank()
        self.size = MPI_COMM_WORLD.Get_size()

        # Test Data Generation
        self._clean_directory()
        self.nlat = 19
        self.nlon = 36
        self.ntime = 10
        self.slices = ['input{}.nc'.format(i) for i in xrange(5)]
        self.scalars = ['scalar{}'.format(i) for i in xrange(2)]
        self.timvars = ['tim{}'.format(i) for i in xrange(2)]
        self.tvmvars = ['tvm{}'.format(i) for i in xrange(2)]
        self.tsvars = ['tsvar{}'.format(i) for i in xrange(4)]
        self.fattrs = {'attr1': 'attribute one',
                       'attr2': 'attribute two'}
        if self.rank == 0:
            for i in xrange(len(self.slices)):

                # Open the file for writing
                fname = self.slices[i]
                fobj = Nio.open_file(fname, 'w')

                # Write attributes to file
                for name, value in self.fattrs.iteritems():
                    setattr(fobj, name, value)

                # Create the dimensions in the file
                fobj.create_dimension('lat', self.nlat)
                fobj.create_dimension('lon', self.nlon)
                fobj.create_dimension('time', None)

                # Create the coordinate variables & add attributes
                lat = fobj.create_variable('lat', 'f', ('lat',))
                lon = fobj.create_variable('lon', 'f', ('lon',))
                time = fobj.create_variable('time', 'f', ('time',))

                # Set the coordinate variable attributes
                setattr(lat, 'long_name', 'latitude')
                setattr(lon, 'long_name', 'longitude')
                setattr(time, 'long_name', 'time')
                setattr(lat, 'units', 'degrees north')
                setattr(lon, 'units', 'degrees east')
                setattr(time, 'units', 'days from 01-01-0001')

                # Set the values of the coordinate variables
                lat[:] = np.linspace(-90, 90, self.nlat, dtype=np.float32)
                lon[:] = np.linspace(-180, 180, self.nlon, endpoint=False, dtype=np.float32)
                time[:] = np.arange(i * self.ntime, (i + 1) * self.ntime, dtype=np.float32)

                # Create the scalar variables
                for n in xrange(len(self.scalars)):
                    vname = self.scalars[n]
                    v = fobj.create_variable(vname, 'd', ())
                    setattr(v, 'long_name', 'scalar{}'.format(n))
                    setattr(v, 'units', '[{}]'.format(vname))
                    v.assign_value(np.float64(n * 10))

                # Create the time-invariant metadata variables
                for n in xrange(len(self.timvars)):
                    vname = self.timvars[n]
                    v = fobj.create_variable(vname, 'd', ('lat', 'lon'))
                    setattr(v, 'long_name', 'time-invariant metadata {}'.format(n))
                    setattr(v, 'units', '[{}]'.format(vname))
                    v[:] = np.ones((self.nlat, self.nlon), dtype=np.float64)

                # Create the time-variant metadata variables
                for n in xrange(len(self.tvmvars)):
                    vname = self.tvmvars[n]
                    v = fobj.create_variable(vname, 'd', ('time', 'lat', 'lon'))
                    setattr(v, 'long_name', 'time-variant metadata {}'.format(n))
                    setattr(v, 'units', '[{}]'.format(vname))
                    v[:] = np.ones((self.ntime, self.nlat, self.nlon), dtype=np.float64)

                # Create the time-series variables
                for n in xrange(len(self.tsvars)):
                    vname = self.tsvars[n]
                    v = fobj.create_variable(vname, 'd', ('time', 'lat', 'lon'))
                    setattr(v, 'long_name', 'time-series variable {}'.format(n))
                    setattr(v, 'units', '[{}]'.format(vname))
                    v[:] = np.ones((self.ntime, self.nlat, self.nlon), dtype=np.float64)

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

    def _test_create_reshaper(self, serial, verbosity, wmode):
        self._test_header(("create_reshaper(serial={}, verbosity={}, "
                           "wmode={!r})").format(serial, verbosity, wmode))
        if not (serial and self.rank > 0):
            spec = Specifier(infiles=self.slices, ncfmt='netcdf',
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

    def _check_outfiles(self, infiles, prefix, suffix, metadata, once, numfact=1):
        if self.rank == 0:

            nsteps = 0
            for infile in infiles:
                ncinp = Nio.open_file(infile, 'r')
                nsteps += ncinp.dimensions['time']
                ncinp.close()
            nsteps *= numfact

            ncinp = Nio.open_file(infiles[0], 'r')

            scalars = [v for v in ncinp.variables
                       if ncinp.variables[v].dimensions == ()]
            tivars = [v for v in ncinp.variables
                      if 'time' not in ncinp.variables[v].dimensions]
            tsvars = [v for v in ncinp.variables
                      if 'time' in ncinp.variables[v].dimensions and v not in metadata]

            if once:
                tsvars.append('once')

            outfiles = ['{}{}{}'.format(prefix, v, suffix) for v in tsvars]

            outdims = {'time': nsteps,
                       'lat': ncinp.dimensions['lat'],
                       'lon': ncinp.dimensions['lon']}

            outmeta = [v for v in ncinp.variables if v not in tsvars]

            for tsvar, outfile in zip(tsvars, outfiles):

                self._assertion("exists({})".format(outfile),
                                exists(outfile), True)

                ncout = Nio.open_file(outfile, 'r')

                self._assertion("{}: attributes equal".format(outfile),
                                ncout.attributes, ncinp.attributes,
                                assertion=self.assertDictEqual)

                for d, v in outdims.iteritems():
                    self._assertion("{}: {} in dimensions".format(outfile, d),
                                    d in ncout.dimensions, True)

                    self._assertion("{}: dimensions[{}]".format(outfile, d),
                                    ncout.dimensions[d], v)

                self._assertion("{}: time unlimited".format(outfile),
                                ncout.unlimited('time'), True)

                if once:
                    all_vars = outmeta if tsvar == 'once' else [tsvar]
                else:
                    all_vars = [tsvar] + outmeta

                self._assertion("{}: variable list".format(outfile),
                                set(ncout.variables.keys()), set(all_vars),
                                assertion=self.assertSetEqual)

                for v in all_vars:
                    self._assertion("{}: {} in variables".format(outfile, v),
                                    v in ncout.variables, True)

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

                ncout.close()

            ncinp.close()

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

    def _test_convert(self, infiles, prefix, suffix, metadata,
                      ncfmt, clevel, serial, verbosity, wmode, once,
                      print_diags=False, numfact=1):
        nfiles = len(infiles)
        ncvers = '3' if ncfmt == 'netcdf' else ('4c' if ncfmt == 'netcdf4c' else '4')
        self._test_header(("convert() - {} infile(s), NC{}-CL{}, serial={},{}"
                           "            verbosity={}, wmode={!r}, once={}"
                           "").format(nfiles, ncvers, clevel, serial, eol,
                                      verbosity, wmode, once))
        self._run_convert(infiles, prefix, suffix, metadata, ncfmt, clevel,
                          serial, verbosity, wmode, once, print_diags)
        self._check_outfiles(infiles, prefix, suffix, metadata, once, numfact)

    def _test_convert_no_output(self, infiles, prefix, suffix, metadata,
                                ncfmt, clevel, serial, wmode, once, numfact=1):
        oldout = sys.stdout
        newout = StringIO()
        sys.stdout = newout
        self._run_convert(infiles, prefix, suffix, metadata, ncfmt, clevel,
                          serial, 0, wmode, once, print_diags=False)
        actual = newout.getvalue()
        self._assertion("stdout empty", actual, '')
        sys.stdout = oldout
        self._check_outfiles(infiles, prefix, suffix, metadata, once, numfact)

    def testReshaperConvert_All_NC3_CL0_SER_V0_W(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        self._test_convert_no_output(infiles=infiles, prefix='out.', suffix='.nc',
                                     metadata=mdata, ncfmt='netcdf', clevel=0,
                                     serial=True, wmode='w', once=False)

    def testReshaperConvert_1_NC3_CL0_SER_V0_W(self):
        infiles = self.slices[0:1]
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        self._test_convert_no_output(infiles=infiles, prefix='out.', suffix='.nc',
                                     metadata=mdata, ncfmt='netcdf', clevel=0,
                                     serial=True, wmode='w', once=False)

    def testReshaperConvert_All_NC4_CL1_SER_V0_W(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        self._test_convert_no_output(infiles=infiles, prefix='out.', suffix='.nc',
                                     metadata=mdata, ncfmt='netcdf4', clevel=1,
                                     serial=True, wmode='w', once=False)

    def testReshaperConvert_All_NC3_CL0_PAR_V1_W(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        self._test_convert(infiles=infiles, prefix='out.', suffix='.nc',
                           metadata=mdata, ncfmt='netcdf', clevel=0,
                           serial=False, verbosity=1, wmode='w', once=False,
                           print_diags=True)

    def testReshaperConvert_All_NC3_CL0_PAR_V1_W_ONCE(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        self._test_convert(infiles=infiles, prefix='out.', suffix='.nc',
                           metadata=mdata, ncfmt='netcdf', clevel=0,
                           serial=False, verbosity=1, wmode='w', once=True,
                           print_diags=True)

    def testReshaperConvert_All_NC3_CL0_PAR_V1_O(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        convert_args = {'infiles': infiles, 'prefix': 'out.', 'suffix': '.nc',
                        'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                        'serial': False, 'verbosity': 1, 'once': False}
        self._run_convert(wmode='w', print_diags=False, **convert_args)
        self._test_convert(wmode='o', print_diags=True, **convert_args)

    def testReshaperConvert_All_NC3_CL0_PAR_V1_O_ONCE(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        convert_args = {'infiles': infiles, 'prefix': 'out.', 'suffix': '.nc',
                        'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                        'serial': False, 'verbosity': 1, 'once': True}
        self._run_convert(wmode='w', print_diags=False, **convert_args)
        self._test_convert(wmode='o', print_diags=True, **convert_args)

    def testReshaperConvert_All_NC3_CL0_PAR_V1_S(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        convert_args = {'infiles': infiles, 'prefix': 'out.', 'suffix': '.nc',
                        'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                        'serial': False, 'verbosity': 1, 'once': False}
        self._run_convert(wmode='w', print_diags=False, **convert_args)
        self._test_convert(wmode='s', print_diags=True, **convert_args)

    def testReshaperConvert_All_NC3_CL0_PAR_V1_S_ONCE(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        convert_args = {'infiles': infiles, 'prefix': 'out.', 'suffix': '.nc',
                        'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                        'serial': False, 'verbosity': 1, 'once': True}
        self._run_convert(wmode='w', print_diags=False, **convert_args)
        self._test_convert(wmode='s', print_diags=True, **convert_args)

    def testReshaperConvert_All_NC3_CL0_PAR_V3_A(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        convert_args = {'infiles': infiles, 'prefix': 'out.', 'suffix': '.nc',
                        'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                        'serial': False, 'verbosity': 3, 'once': False}
        self._run_convert(wmode='w', print_diags=False, **convert_args)
        self._test_convert(wmode='a', print_diags=True, numfact=2, **convert_args)

    def testReshaperConvert_All_NC3_CL0_PAR_V3_A_ONCE(self):
        infiles = self.slices
        mdata = [v for v in self.tvmvars]
        mdata.append('time')
        convert_args = {'infiles': infiles, 'prefix': 'out.', 'suffix': '.nc',
                        'metadata': mdata, 'ncfmt': 'netcdf', 'clevel': 0,
                        'serial': False, 'verbosity': 3, 'once': True}
        self._run_convert(wmode='w', print_diags=False, **convert_args)
        self._test_convert(wmode='a', print_diags=True, numfact=2, **convert_args)

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
