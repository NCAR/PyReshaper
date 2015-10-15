"""
Copyright 2015, University Corporation for Atmospheric Research
See LICENSE.txt for details
"""

import imp
import sys
import unittest

from StringIO import StringIO

s2smake = imp.load_source('s2smake', '../../../scripts/s2smake')


class s2smakeTest(unittest.TestCase):

    def setUp(self):
        self.saved_stdout = sys.stdout
        self.saved_stderr = sys.stderr
        self.stdout = StringIO()
        self.stderr = StringIO()
        sys.stdout = self.stdout
        sys.stderr = self.stderr

    def tearDown(self):
        sys.stdout = self.saved_stdout
        sys.stderr = self.saved_stderr

    def test_CLI_empty(self):
        argv = []
        self.assertRaises(SystemExit, s2smake.cli, argv)
        outstr = self.stdout.getvalue().strip()
        print outstr
        print '*' * 40
        errstr = self.stderr.getvalue().strip()
        print errstr

    def test_CLI_help(self):
        argv = ['-h']
        self.assertRaises(SystemExit, s2smake.cli, argv)


if __name__ == "__main__":
    unittest.main()
