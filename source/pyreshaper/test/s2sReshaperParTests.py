"""
Parallel Tests for the Reshaper class
"""

import unittest

from os import linesep as eol
from pyreshaper.reshaper import Slice2SeriesReshaper, create_reshaper
from pyreshaper.specification import Slice2SeriesSpecifier

from mpi4py import MPI
MPI_COMM_WORLD = MPI.COMM_WORLD


class S2SReshaperParTests(unittest.TestCase):

    def setUp(self):
        self.infiles = ['input1.nc', 'input2.nc', 'input3.nc']
        self.ncfmt = 'netcdf'
        self.prefix = 'output.'
        self.suffix = '.nc'
        self.metadata = ['meta1', 'meta2']

        self.spec = Slice2SeriesSpecifier(
            infiles=self.infiles, ncfmt=self.ncfmt, prefix=self.prefix,
            suffix=self.suffix, metadata=self.metadata)

        self.size = MPI_COMM_WORLD.Get_size()
        self.rank = MPI_COMM_WORLD.Get_rank()

    def tearDown(self):
        pass

    def _info_msg(self, name, data, actual, expected):
        rknm = ''.join(
            ['[', str(self.rank), '/', str(self.size), '] ', str(name)])
        spcr = ' ' * len(rknm)
        msg = ''.join([eol,
                       rknm, ' - Input: ', str(data), eol,
                       spcr, ' - Actual:   ', str(actual), eol,
                       spcr, ' - Expected: ', str(expected)])
        return msg

    def testCreateReshaperType(self):
        actual = create_reshaper(
            self.spec, serial=False, verbosity=2, once=False)
        expected = Slice2SeriesReshaper
        msg = self._info_msg("type(create_reshaper)",
                             None, type(actual), expected)
        self.assertEqual(type(actual), expected, msg)


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
