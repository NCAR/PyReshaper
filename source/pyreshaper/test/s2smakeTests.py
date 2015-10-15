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
        print '_' * 7
        print 'STDOUT:'
        print '-' * 7
        print self.outstr
        print '_' * 7
        print 'STDERR:'
        print '-' * 7
        print self.errstr

    def test_CLI_empty(self):
        argv = []
        self.assertRaises(SystemExit, s2smake.cli, argv)
        self.outstr = self.stdout.getvalue().strip()
        self.assertEqual(self.outstr, '', 'STDOUT Unexpected')
        self.errstr = self.stderr.getvalue().strip()

    def test_CLI_help(self):
        argv = ['-h']
        self.assertRaises(SystemExit, s2smake.cli, argv)
        self.outstr = self.stdout.getvalue().strip()
        self.errstr = self.stderr.getvalue().strip()
        self.assertEqual(self.errstr, '', 'STDERR Unexpected')

    def test_CLI_defaults(self):
        argv = ['s2smakeTests.py']
        args = s2smake.cli(argv)
        self.assertEqual(args.compression_level, 1,
                         'Default compression level is not 1')
        self.assertListEqual(args.infiles, argv,
                             'Default infiles list is not {}'.format(argv))
        self.assertListEqual(args.metadata, [],
                             'Default metadata list is not []')
        self.assertEqual(args.netcdf_format, 'netcdf4',
                         'Default NetCDF format is not "netcdf4"')
        self.assertEqual(args.output_prefix, 'tseries.',
                         'Default output prefix is not "tseries."')
        self.assertEqual(args.output_suffix, '.nc',
                         'Default output suffix is not ".nc"')
        self.assertEqual(args.specfile, 'input.s2s',
                         'Default output suffix is not ".nc"')
        self.outstr = self.stdout.getvalue().strip()
        self.assertEqual(self.outstr, '', 'STDOUT Unexpected')
        self.errstr = self.stderr.getvalue().strip()
        self.assertEqual(self.errstr, '', 'STDERR Unexpected')

    def test_CLI_set_all_short(self):
        clevel = 3
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'mysepc.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        infiles = ['s2smakeTests.py', 'specificationTests.py']

        argv = ['-c', str(clevel), '-f', ncfmt]
        for md in metadata:
            argv.extend(['-m', md])
        argv.extend(['-o', specfile, '-p', prefix, '-s', suffix])
        argv.extend(infiles)
        args = s2smake.cli(argv)

        self.assertEqual(args.compression_level, clevel,
                         'Default compression level is not {!r}'.format(clevel))
        self.assertListEqual(args.infiles, infiles,
                             'Default infiles list is not {}'.format(infiles))
        self.assertListEqual(args.metadata, metadata,
                             'Default metadata list is not {}'.format(metadata))
        self.assertEqual(args.netcdf_format, ncfmt,
                         'Default NetCDF format is not {!r}'.format(ncfmt))
        self.assertEqual(args.output_prefix, prefix,
                         'Default output prefix is not {!r}'.format(prefix))
        self.assertEqual(args.output_suffix, suffix,
                         'Default output suffix is not {!r}'.format(suffix))
        self.assertEqual(args.specfile, specfile,
                         'Default output suffix is not {!r}'.format(specfile))
        self.outstr = self.stdout.getvalue().strip()
        self.assertEqual(self.outstr, '', 'STDOUT Unexpected')
        self.errstr = self.stderr.getvalue().strip()
        self.assertEqual(self.errstr, '', 'STDERR Unexpected')

    def test_CLI_set_all_long(self):
        clevel = 3
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'mysepc.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        infiles = ['s2smakeTests.py', 'specificationTests.py']

        argv = ['--compression_level', str(clevel), '--netcdf_format', ncfmt]
        for md in metadata:
            argv.extend(['--metadata', md])
        argv.extend(['--specfile', specfile, '--output_prefix', prefix,
                     '--output_suffix', suffix])
        argv.extend(infiles)
        args = s2smake.cli(argv)

        self.assertEqual(args.compression_level, clevel,
                         'Default compression level is not {!r}'.format(clevel))
        self.assertListEqual(args.infiles, infiles,
                             'Default infiles list is not {}'.format(infiles))
        self.assertListEqual(args.metadata, metadata,
                             'Default metadata list is not {}'.format(metadata))
        self.assertEqual(args.netcdf_format, ncfmt,
                         'Default NetCDF format is not {!r}'.format(ncfmt))
        self.assertEqual(args.output_prefix, prefix,
                         'Default output prefix is not {!r}'.format(prefix))
        self.assertEqual(args.output_suffix, suffix,
                         'Default output suffix is not {!r}'.format(suffix))
        self.assertEqual(args.specfile, specfile,
                         'Default output suffix is not {!r}'.format(specfile))
        self.outstr = self.stdout.getvalue().strip()
        self.assertEqual(self.outstr, '', 'STDOUT Unexpected')
        self.errstr = self.stderr.getvalue().strip()
        self.assertEqual(self.errstr, '', 'STDERR Unexpected')

if __name__ == "__main__":
    unittest.main()
