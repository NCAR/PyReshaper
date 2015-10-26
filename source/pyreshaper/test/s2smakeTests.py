"""
Copyright 2015, University Corporation for Atmospheric Research
See LICENSE.txt for details
"""

import os
import imp
import sys
import unittest
import cPickle as pickle

from pyreshaper.specification import Specifier

s2smake = imp.load_source('s2smake', '../../../scripts/s2smake')


class s2smakeTest(unittest.TestCase):

    def test_CLI_empty(self):
        argv = []
        self.assertRaises(SystemExit, s2smake.cli, argv)

    def test_CLI_help(self):
        argv = ['-h']
        self.assertRaises(SystemExit, s2smake.cli, argv)

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

    def test_CLI_set_all_short(self):
        clevel = 3
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'myspec.s2s'
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

    def test_CLI_set_all_long(self):
        clevel = 3
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'myspec.s2s'
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

    def test_main_defaults(self):
        argv = ['s2smakeTests.py']
        specfile = 'input.s2s'
        if os.path.exists(specfile):
            os.remove(specfile)
        s2smake.main(argv)
        self.assertTrue(os.path.exists(specfile), 'Specfile not found')
        spec = pickle.load(open(specfile, 'rb'))
        os.remove(specfile)
        self.assertTrue(isinstance(spec, Specifier),
                        'Specfile does not contain a Specifier')
        self.assertEqual(spec.compression_level, 1,
                         'Default compression level is not 1')
        self.assertListEqual(spec.input_file_list, argv,
                             'Default infiles list is not {}'.format(argv))
        self.assertListEqual(spec.time_variant_metadata, [],
                             'Default metadata list is not []')
        self.assertEqual(spec.netcdf_format, 'netcdf4',
                         'Default NetCDF format is not "netcdf4"')
        self.assertEqual(spec.output_file_prefix, os.path.abspath('tseries.'),
                         'Default output prefix is not "tseries."')
        self.assertEqual(spec.output_file_suffix, '.nc',
                         'Default output suffix is not ".nc"')

    def test_main_set_all_short(self):
        clevel = 3
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'myspec.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        infiles = ['s2smakeTests.py', 'specificationTests.py']

        argv = ['-c', str(clevel), '-f', ncfmt]
        for md in metadata:
            argv.extend(['-m', md])
        argv.extend(['-o', specfile, '-p', prefix, '-s', suffix])
        argv.extend(infiles)

        if os.path.exists(specfile):
            os.remove(specfile)
        s2smake.main(argv)

        self.assertTrue(os.path.exists(specfile), 'Specfile not found')
        spec = pickle.load(open(specfile, 'rb'))
        os.remove(specfile)

        self.assertTrue(isinstance(spec, Specifier),
                        'Specfile does not contain a Specifier')
        self.assertEqual(spec.compression_level, clevel,
                         'Default compression level is not {!r}'.format(clevel))
        self.assertListEqual(spec.input_file_list, infiles,
                             'Default infiles list is not {}'.format(infiles))
        self.assertListEqual(spec.time_variant_metadata, metadata,
                             'Default metadata list is not {}'.format(metadata))
        self.assertEqual(spec.netcdf_format, ncfmt,
                         'Default NetCDF format is not {!r}'.format(ncfmt))
        self.assertEqual(spec.output_file_prefix, os.path.abspath(prefix),
                         'Default output prefix is not {!r}'.format(prefix))
        self.assertEqual(spec.output_file_suffix, suffix + '.nc',
                         'Default output suffix is not {!r}'.format(suffix))

    def test_main_set_all_long(self):
        clevel = 3
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'myspec.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        infiles = ['s2smakeTests.py', 'specificationTests.py']

        argv = ['--compression_level', str(clevel), '--netcdf_format', ncfmt]
        for md in metadata:
            argv.extend(['--metadata', md])
        argv.extend(['--specfile', specfile, '--output_prefix', prefix,
                     '--output_suffix', suffix])
        argv.extend(infiles)

        if os.path.exists(specfile):
            os.remove(specfile)
        s2smake.main(argv)

        self.assertTrue(os.path.exists(specfile), 'Specfile not found')
        spec = pickle.load(open(specfile, 'rb'))
        os.remove(specfile)

        self.assertTrue(isinstance(spec, Specifier),
                        'Specfile does not contain a Specifier')
        self.assertEqual(spec.compression_level, clevel,
                         'Default compression level is not {!r}'.format(clevel))
        self.assertListEqual(spec.input_file_list, infiles,
                             'Default infiles list is not {}'.format(infiles))
        self.assertListEqual(spec.time_variant_metadata, metadata,
                             'Default metadata list is not {}'.format(metadata))
        self.assertEqual(spec.netcdf_format, ncfmt,
                         'Default NetCDF format is not {!r}'.format(ncfmt))
        self.assertEqual(spec.output_file_prefix, os.path.abspath(prefix),
                         'Default output prefix is not {!r}'.format(prefix))
        self.assertEqual(spec.output_file_suffix, suffix + '.nc',
                         'Default output suffix is not {!r}'.format(suffix))

if __name__ == "__main__":
    unittest.main()
