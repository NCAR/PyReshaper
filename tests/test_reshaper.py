"""
Copyright 2020, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

import inspect
import sys
import unittest
from glob import glob
from io import StringIO
from os import linesep as eol, remove
from os.path import exists

from pyreshaper.reshaper import Reshaper, create_reshaper
from pyreshaper.specification import Specifier

from .checks import check_outfile, check_var_in
from .data import config, make


class CommonTestsBase(object):
    def setUp(self):

        # Default arguments for testing
        self.spec_args = {
            'infiles': config.slices,
            'ncfmt': 'netcdf4',
            'compression': 0,
            'prefix': 'out.',
            'suffix': '.nc',
            'timeseries': None,
            'metadata': [v for v in config.tvmvars] + ['time'] + [v for v in config.chvars],
            'meta1d': False,
            'metafile': None,
        }
        self.create_args = {
            'serial': False,
            'verbosity': 1,
            'wmode': 'w',
            'once': False,
            'simplecomm': None,
        }
        self.convert_args = {'output_limit': 0, 'rchunks': None, 'wchunks': None}

        # Test Data Generation
        self.clean()
        make.generate_data()

    def clean(self):
        for ncfile in glob('*.nc'):
            remove(ncfile)

    def header(self):
        mf = len(config.slices)
        mt = len(config.tsvars)
        nf = len(self.spec_args['infiles'])
        nt = mt if self.spec_args['timeseries'] is None else len(self.spec_args['timeseries'])

        hline = '-' * 100
        hdrstr = [
            '',
            hline,
            '{}.{}:'.format(self.__class__.__name__, inspect.stack()[1][3]),
            '',
            '   specifier({}/{} infile(s), {}/{} TSV(s), ncfmt={ncfmt}, compression={compression}, meta1d={meta1d}, backend={backend}, metafile={metafile})'.format(
                nf, mf, nt, mt, **self.spec_args
            ),
            '   create(serial={serial}, verbosity={verbosity}, wmode={wmode}, once={once}, simplecomm={simplecomm})'.format(
                **self.create_args
            ),
            '   convert(output_limit={output_limit}, rchunks={rchunks}, wchunks={wchunks})'.format(
                **self.convert_args
            ),
            hline,
            '',
        ]
        print(eol.join(hdrstr))

    def check(self, tsvar):
        args = {}
        args.update(self.spec_args)
        args.update(self.create_args)
        args.update(self.convert_args)
        assertions_dict = check_outfile(tsvar=tsvar, **args)
        failed_assertions = [key for key, value in assertions_dict.items() if value is False]
        assert_msgs = ['Output file check for variable {0!r}:'.format(tsvar)]
        assert_msgs.extend(['   {0}'.format(assrt) for assrt in failed_assertions])
        assert len(failed_assertions) == 0, eol.join(assert_msgs)

    def convert(self, print_diags=False):
        if self.create_args.get('verbosity', 1) == 0:
            oldout = sys.stdout
            newout = StringIO()
            sys.stdout = newout

        spec = Specifier(**self.spec_args)
        rshpr = create_reshaper(spec, **self.create_args)
        assert type(rshpr) == Reshaper, 'type(reshaper) failure'
        rshpr.convert(**self.convert_args)

        if self.create_args.get('verbosity', 1) == 0:
            actual = newout.getvalue()
            assert actual == '', 'stdout should be empty'
            sys.stdout = oldout

        if print_diags:
            rshpr.print_diagnostics()

    def test_defaults(self):
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_I1(self):
        self.spec_args['infiles'] = config.slices[1:2]
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_V0(self):
        self.create_args['verbosity'] = 0
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_V3(self):
        self.create_args['verbosity'] = 3
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_TSV2(self):
        self.spec_args['timeseries'] = config.tsvars[1:3] + ['tsvarX']
        self.header()
        self.convert()
        for tsvar in self.spec_args['timeseries']:
            if tsvar in config.tsvars:
                self.check(tsvar)
            else:
                fname = self.spec_args['prefix'] + tsvar + self.spec_args['suffix']
                assert not exists(fname), 'File {0!r} should not exist'.format(fname)

    def test_exclude(self):
        self.spec_args['exclude_list'] = config.timvars[0:1]
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            fname = self.spec_args['prefix'] + tsvar + self.spec_args['suffix']
            for timvar in config.timvars:
                if timvar in self.spec_args['exclude_list']:
                    xassert = self.assertFalse
                else:
                    xassert = self.assertTrue
                xassert(check_var_in(timvar, fname))

    def test_NC3(self):
        self.spec_args['ncfmt'] = 'netcdf'
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_ser(self):
        self.create_args['serial'] = True
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_CL1(self):
        self.spec_args['compression'] = 1
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_read_chunking(self):
        self.convert_args['rchunks'] = {'lat': 1, 'time': config.ntime}
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_write_chunking(self):
        self.convert_args['wchunks'] = {'lat': 1, 'time': config.ntime}
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_meta1d(self):
        self.spec_args['meta1d'] = True
        self.spec_args['metadata'] = [v for v in config.tvmvars]
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_metafile(self):
        self.spec_args['metafile'] = 'metafile.nc'
        self.spec_args['metadata'] = [v for v in config.tvmvars]
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_once(self):
        self.create_args['once'] = True
        self.header()
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_overwrite(self):
        self.create_args['wmode'] = 'o'
        self.header()
        self.create_args['verbosity'] = 0
        self.convert()
        self.create_args['verbosity'] = 1
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_skip(self):
        self.create_args['wmode'] = 's'
        self.header()
        self.create_args['verbosity'] = 0
        self.convert()
        self.create_args['verbosity'] = 1
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_append(self):
        self.create_args['wmode'] = 'a'
        self.header()
        self.create_args['wmode'] = 'w'
        self.spec_args['infiles'] = config.slices[0:2]
        self.convert()
        self.create_args['wmode'] = 'a'
        self.spec_args['infiles'] = config.slices[2:]
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_append_missing(self):
        missing = config.tsvars[2]
        self.create_args['wmode'] = 'a'
        self.header()
        self.create_args['wmode'] = 'w'
        self.spec_args['infiles'] = config.slices[0:2]
        self.convert()
        remove(self.spec_args['prefix'] + missing + self.spec_args['suffix'])

        self.create_args['wmode'] = 'a'
        self.spec_args['infiles'] = config.slices[2:]
        self.convert()
        for tsvar in config.tsvars:
            if tsvar == missing:
                self.spec_args['infiles'] = config.slices[2:]
            else:
                self.spec_args['infiles'] = config.slices
            self.check(tsvar)


class NetCDF4Tests(unittest.TestCase, CommonTestsBase):
    """NetCDF4 Python Tests"""

    def setUp(self):
        CommonTestsBase.setUp(self)
        self.spec_args['backend'] = 'netCDF4'

    def tearDown(self):
        self.clean()


class NioTests(unittest.TestCase, CommonTestsBase):
    """PyNIO Tests"""

    def setUp(self):
        CommonTestsBase.setUp(self)
        self.spec_args['backend'] = 'Nio'

    def tearDown(self):
        self.clean()
