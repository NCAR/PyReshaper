"""
Parallel Tests for the Reshaper class

Copyright 2015, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

import unittest

from os import linesep as eol
from os import remove
from os.path import exists
from mpi4py import MPI

import Nio

from pyreshaper.reshaper import Slice2SeriesReshaper, create_reshaper
from pyreshaper.specification import Specifier
import makeTestData

MPI_COMM_WORLD = MPI.COMM_WORLD


class S2SReshaperParTests(unittest.TestCase):

    def setUp(self):
        self.rank = MPI_COMM_WORLD.Get_rank()
        self.size = MPI_COMM_WORLD.Get_size()
        self.nlat = 19
        self.nlon = 36
        self.ntime = 10
        self.infiles = ['input{}.nc'.format(i) for i in xrange(5)]
        self.scalars = ['scalar{}'.format(i) for i in xrange(2)]
        self.timvars = ['tim{}'.format(i) for i in xrange(2)]
        self.tvmvars = ['tvm{}'.format(i) for i in xrange(2)]
        self.tsvars = ['tsvar{}'.format(i) for i in xrange(4)]
        self.fattrs = {'attr1': 'attribute one', 'attr2': 'attribute two'}
        if self.rank == 0:
            makeTestData.make_data(nlat=self.nlat,
                                   nlon=self.nlon,
                                   ntime=self.ntime,
                                   slices=self.infiles,
                                   scalars=self.scalars,
                                   timvars=self.timvars,
                                   tvmvars=self.tvmvars,
                                   tsvars=self.tsvars,
                                   fattrs=self.fattrs)
        MPI_COMM_WORLD.Barrier()
        self.ncfmt = 'netcdf'
        self.compression = 0
        self.prefix = 'output.'
        self.suffix = '.nc'
        self.metadata = [v for v in self.tvmvars]
        self.metadata.append('time')
        self.spec = Specifier(
            infiles=self.infiles, ncfmt=self.ncfmt, compression=self.compression,
            prefix=self.prefix, suffix=self.suffix, metadata=self.metadata)
        self.rshpr = create_reshaper(self.spec, serial=False, verbosity=3,
                                     wmode='w')
        self.outfiles = ['{}{}{}'.format(self.prefix, v, self.suffix)
                         for v in self.tsvars]

    def tearDown(self):
        if self.rank == 0:
            for infile in self.infiles:
                if exists(infile):
                    remove(infile)
            for outfile in self.outfiles:
                if exists(outfile):
                    remove(outfile)
        MPI_COMM_WORLD.Barrier()

    def _info_msg(self, name, data, actual, expected, show=True):
        rknm = '[{}/{}] {}'.format(self.rank, self.size, name)
        spcr = ' ' * len(rknm)
        msg = ''.join([eol,
                       rknm, ' - Input: ', str(data), eol,
                       spcr, ' - Actual:   ', str(actual), eol,
                       spcr, ' - Expected: ', str(expected)])
        if show:
            print msg
        return msg

    def testCreateReshaperType(self):
        actual = type(self.rshpr)
        expected = Slice2SeriesReshaper
        msg = self._info_msg("type(reshaper)",
                             None, actual, expected)
        self.assertEqual(actual, expected, msg)

    def testReshaperConvert(self):
        self.rshpr.convert()
        self.rshpr.print_diagnostics()
        if self.rank == 0:
            for outfile, tsvar in zip(self.outfiles, self.tsvars):

                actual = exists(outfile)
                expected = True
                msg = self._info_msg("exists({})".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                ncout = Nio.open_file(outfile, 'r')

                actual = 'lat' in ncout.dimensions
                expected = True
                msg = self._info_msg("{}: lat in dimensions".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                actual = 'lon' in ncout.dimensions
                expected = True
                msg = self._info_msg("{}: lon in dimensions".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                actual = 'time' in ncout.dimensions
                expected = True
                msg = self._info_msg("{}: time in dimensions".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                actual = ncout.dimensions['lat']
                expected = self.nlat
                msg = self._info_msg("{}: dimensions[lat]".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                actual = ncout.dimensions['lon']
                expected = self.nlon
                msg = self._info_msg("{}: dimensions[lon]".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                nsteps = len(self.infiles) * self.ntime

                actual = ncout.dimensions['time']
                expected = nsteps
                msg = self._info_msg("{}: dimensions[time]".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                actual = ncout.unlimited('time')
                expected = True
                msg = self._info_msg("{}: time unlimited".format(outfile),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                for v in self.scalars:
                    actual = v in ncout.variables
                    expected = True
                    msg = self._info_msg("{}: var: {}".format(outfile, v),
                                         None, actual, expected)
                    self.assertEqual(actual, expected, msg)

                    actual = ncout.variables[v].dimensions
                    expected = ()
                    msg = self._info_msg("{}: dim: {}".format(outfile, v),
                                         None, actual, expected)
                    self.assertTupleEqual(actual, expected, msg)

                for v in self.timvars:
                    actual = v in ncout.variables
                    expected = True
                    msg = self._info_msg("{}: var: {}".format(outfile, v),
                                         None, actual, expected)
                    self.assertEqual(actual, expected, msg)

                    actual = ncout.variables[v].dimensions
                    expected = ('lat', 'lon')
                    msg = self._info_msg("{}: dim: {}".format(outfile, v),
                                         None, actual, expected)
                    self.assertTupleEqual(actual, expected, msg)

                for v in self.tvmvars:
                    actual = v in ncout.variables
                    expected = True
                    msg = self._info_msg("{}: var: {}".format(outfile, v),
                                         None, actual, expected)
                    self.assertEqual(actual, expected, msg)

                    actual = ncout.variables[v].dimensions
                    expected = ('time', 'lat', 'lon')
                    msg = self._info_msg("{}: dim: {}".format(outfile, v),
                                         None, actual, expected)
                    self.assertTupleEqual(actual, expected, msg)

                actual = tsvar in ncout.variables
                expected = True
                msg = self._info_msg("{}: var: {}".format(outfile, tsvar),
                                     None, actual, expected)
                self.assertEqual(actual, expected, msg)

                actual = ncout.variables[tsvar].dimensions
                expected = ('time', 'lat', 'lon')
                msg = self._info_msg("{}: dim: {}".format(outfile, tsvar),
                                     None, actual, expected)
                self.assertTupleEqual(actual, expected, msg)

                actual = set(ncout.variables.keys())
                expected = set(['lat', 'lon', 'time', tsvar])
                expected.update(self.scalars)
                expected.update(self.timvars)
                expected.update(self.tvmvars)
                msg = self._info_msg("{}: variable list".format(outfile),
                                     None, actual, expected)
                self.assertSetEqual(actual, expected, msg)

                ncout.close()

        MPI_COMM_WORLD.Barrier()


if __name__ == "__main__":
    hline = '=' * 70
    if MPI_COMM_WORLD.Get_rank() == 0:
        print hline
        print 'STANDARD OUTPUT FROM ALL TESTS:'
        print hline
    MPI_COMM_WORLD.Barrier()

    from cStringIO import StringIO
    mystream = StringIO()
    tests = unittest.TestLoader().loadTestsFromTestCase(S2SReshaperParTests)
    unittest.TextTestRunner(stream=mystream).run(tests)
    MPI_COMM_WORLD.Barrier()

    results = MPI_COMM_WORLD.gather(mystream.getvalue())
    if MPI_COMM_WORLD.Get_rank() == 0:
        for rank, result in enumerate(results):
            print hline
            print 'TESTS RESULTS FOR RANK ' + str(rank) + ':'
            print hline
            print str(result)
