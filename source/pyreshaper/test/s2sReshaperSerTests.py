"""
Parallel Tests for the Reshaper class
"""

import unittest

from os import linesep as eol
from pyreshaper.reshaper import Slice2SeriesReshaper, create_reshaper
from pyreshaper.specification import Specifier


class S2SReshaperParTests(unittest.TestCase):

    def setUp(self):
        self.infiles = ['input1.nc', 'input2.nc', 'input3.nc']
        self.ncfmt = 'netcdf'
        self.prefix = 'output.'
        self.suffix = '.nc'
        self.metadata = ['meta1', 'meta2']

        self.spec = Specifier(
            infiles=self.infiles, ncfmt=self.ncfmt, prefix=self.prefix,
            suffix=self.suffix, metadata=self.metadata)

    def tearDown(self):
        pass

    def _info_msg(self, name, data, actual, expected):
        spcr = ' ' * len(name)
        msg = ''.join([eol,
                       name, ' - Input: ', str(data), eol,
                       spcr, ' - Actual:   ', str(actual), eol,
                       spcr, ' - Expected: ', str(expected)])
        return msg

    def testCreateReshaperType(self):
        actual = create_reshaper(
            self.spec, serial=True, verbosity=2, once=False)
        expected = Slice2SeriesReshaper
        msg = self._info_msg("type(create_reshaper)",
                             None, type(actual), expected)
        self.assertEqual(type(actual), expected, msg)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
