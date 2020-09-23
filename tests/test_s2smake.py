"""
Copyright 2020, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

import os
import pickle

import pytest

from pyreshaper.cli import s2smake
from pyreshaper.specification import Specifier

cwd = os.path.dirname(os.path.realpath(__file__))


class CLITests:
    def test_empty(self):
        argv = []
        with pytest.raises(ValueError):
            s2smake.cli(argv)

    def test_help(self):
        argv = ['-h']
        with pytest.raises(SystemExit):
            s2smake.cli(argv)

    def test_defaults(self):
        argv = ['test_s2smake.py']
        opts, args = s2smake.cli(argv)
        assert opts.compression_level == 1
        assert opts.least_significant_digit is None
        for i1, i2 in zip(args, argv):
            assert i1 == i2
        assert len(opts.metadata) == 0
        assert opts.netcdf_format == 'netcdf4'
        assert opts.output_prefix == 'tseries.'
        assert opts.output_suffix == '.nc'
        assert opts.time_series is None
        assert opts.specfile == 'input.s2s'
        assert not opts.meta1d
        assert opts.metafile is None
        assert opts.exclude == []

    def test_set_all_short(self):
        clevel = 3
        lsigfig = 2
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'myspec.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        xlist = ['x', 'y', 'z']
        infiles = ['test_s2smake.py', 'test_specification.py']

        argv = ['-1']
        argv.extend(['-c', str(clevel), '-d', str(lsigfig), '-f', ncfmt])
        for md in metadata:
            argv.extend(['-m', md])
        for x in xlist:
            argv.extend(['-x', x])
        argv.extend(['-o', specfile, '-p', prefix, '-s', suffix])
        argv.extend(infiles)
        opts, args = s2smake.cli(argv)

        assert opts.compression_level == clevel
        assert opts.least_significant_digit == lsigfig
        assert len(args) == len(infiles)
        for i1, i2 in zip(args, infiles):
            assert i1 == i2
        assert len(opts.metadata) == len(metadata)
        for i1, i2 in zip(opts.metadata, metadata):
            assert i1 == i2
        assert opts.netcdf_format == ncfmt
        assert opts.output_prefix == prefix
        assert opts.output_suffix == suffix
        assert opts.time_series is None
        assert opts.specfile == specfile
        assert opts.meta1d
        assert opts.exclude == xlist

    def test_set_all_long(self):
        clevel = 3
        lsigfig = 2
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        xlist = ['x', 'y', 'z']
        specfile = 'myspec.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        tseries = ['tsvar1', 'tsvar2']
        infiles = ['test_s2smake.py', 'test_specification.py']

        argv = ['--meta1d']
        argv.extend(
            [
                '--compression_level',
                str(clevel),
                '--least_significant_digit',
                str(lsigfig),
                '--netcdf_format',
                ncfmt,
            ]
        )
        for md in metadata:
            argv.extend(['--metadata', md])
        for x in xlist:
            argv.extend(['--exclude', x])
        argv.extend(['--specfile', specfile, '--output_prefix', prefix, '--output_suffix', suffix])
        for ts in tseries:
            argv.extend(['--time_series', ts])
        argv.extend(infiles)
        opts, args = s2smake.cli(argv)

        assert opts.compression_level == clevel
        assert opts.least_significant_digit == lsigfig
        assert len(args) == len(infiles)
        for i1, i2 in zip(args, infiles):
            assert i1 == i2
        for i1, i2 in zip(opts.metadata, metadata):
            assert i1 == i2
        assert opts.netcdf_format == ncfmt
        assert opts.output_prefix == prefix
        assert opts.output_suffix == suffix
        for i1, i2 in zip(opts.time_series, tseries):
            assert i1 == i2
        assert opts.specfile == specfile
        assert opts.meta1d
        assert opts.exclude == xlist


class MainTests(object):
    def test_defaults(self):
        argv = [cwd + '/test_s2smake.py']
        specfile = 'input.s2s'
        if os.path.exists(specfile):
            os.remove(specfile)
        s2smake.main(argv)
        assert os.path.exists(specfile), 'Specfile not found'
        spec = pickle.load(open(specfile, 'rb'))
        os.remove(specfile)
        assert isinstance(spec, Specifier)
        assert spec.compression_level == 1
        assert spec.least_significant_digit is None
        assert len(spec.input_file_list) == len(argv)
        for i1, i2 in zip(spec.input_file_list, argv):
            assert i1 == i2
        assert len(spec.time_variant_metadata) == 0
        assert spec.netcdf_format == 'netcdf4'
        assert spec.output_file_prefix == os.path.abspath('tseries.')
        assert spec.output_file_suffix == '.nc'
        assert spec.time_series is None
        assert spec.assume_1d_time_variant_metadata is False

    def test_set_all_short(self):
        clevel = 3
        lsigfig = 2
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        xlist = ['x', 'y', 'z']
        specfile = 'myspec.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        infiles = [cwd + f for f in ['/test_s2smake.py', '/test_specification.py']]

        argv = ['-1', '-c', str(clevel), '-d', str(lsigfig), '-f', ncfmt]
        for md in metadata:
            argv.extend(['-m', md])
        for x in xlist:
            argv.extend(['-x', x])
        argv.extend(['-o', specfile, '-p', prefix, '-s', suffix])
        argv.extend(infiles)

        if os.path.exists(specfile):
            os.remove(specfile)
        s2smake.main(argv)

        assert os.path.exists(specfile), 'Specfile not found'
        spec = pickle.load(open(specfile, 'rb'))
        os.remove(specfile)

        assert isinstance(spec, Specifier)
        assert spec.compression_level == clevel
        assert spec.least_significant_digit == lsigfig
        assert len(spec.input_file_list) == len(infiles)
        for i1, i2 in zip(spec.input_file_list, infiles):
            assert i1 == i2
        assert len(spec.time_variant_metadata) == len(metadata)
        for i1, i2 in zip(spec.time_variant_metadata, metadata):
            assert i1 == i2
        assert spec.netcdf_format == ncfmt
        assert spec.output_file_prefix == os.path.abspath(prefix)
        assert spec.output_file_suffix == suffix + '.nc'
        assert spec.time_series is None
        assert spec.assume_1d_time_variant_metadata is True
        assert spec.exclude_list == xlist

    def test_set_all_long(self):
        clevel = 3
        lsigfig = 2
        ncfmt = 'netcdf'
        metadata = ['meta1', 'meta2']
        specfile = 'myspec.s2s'
        prefix = 'prefix.'
        suffix = '.suffix'
        tseries = ['tsvar1', 'tsvar2']
        infiles = [cwd + f for f in ['/test_s2smake.py', '/test_specification.py']]

        argv = [
            '--meta1d',
            '--compression_level',
            str(clevel),
            '--least_significant_digit',
            str(lsigfig),
            '--netcdf_format',
            ncfmt,
        ]
        for md in metadata:
            argv.extend(['--metadata', md])
        argv.extend(['--specfile', specfile, '--output_prefix', prefix, '--output_suffix', suffix])
        for ts in tseries:
            argv.extend(['--time_series', ts])
        argv.extend(infiles)

        if os.path.exists(specfile):
            os.remove(specfile)
        s2smake.main(argv)

        assert os.path.exists(specfile), 'Specfile not found'
        spec = pickle.load(open(specfile, 'rb'))
        os.remove(specfile)

        assert isinstance(spec, Specifier)
        assert spec.compression_level == clevel
        assert spec.least_significant_digit == lsigfig
        assert len(spec.input_file_list) == len(infiles)
        for i1, i2 in zip(spec.input_file_list, infiles):
            assert i1 == i2
        assert len(spec.time_variant_metadata) == len(metadata)
        for i1, i2 in zip(spec.time_variant_metadata, metadata):
            assert i1 == i2
        assert spec.netcdf_format == ncfmt
        assert spec.output_file_prefix == os.path.abspath(prefix)
        assert spec.output_file_suffix == suffix + '.nc'
        for i1, i2 in zip(spec.time_series, tseries):
            assert i1 == i2
        assert spec.assume_1d_time_variant_metadata is True
