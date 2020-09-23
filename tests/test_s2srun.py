"""
Copyright 2020, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

import unittest
from glob import glob
from os import linesep as eol, remove
from os.path import exists

import pytest

from pyreshaper.cli import s2srun
from pyreshaper.specification import Specifier

from .checks import check_outfile
from .data import config, make


class CLITests:
    @pytest.fixture(autouse=True)
    def init(self):
        self.run_args = {
            'serial': False,
            'rchunks': {'y': 5},
            'wchunks': {'x': 4, 'y': 10},
            'limit': 0,
            'verbosity': 1,
            'write_mode': 'w',
            'once': False,
            'specfile': 'input.s2s',
        }

    def longargs(self):
        argv = [
            '--verbosity',
            str(self.run_args['verbosity']),
            '--write_mode',
            self.run_args['write_mode'],
        ]
        if self.run_args['limit'] > 0:
            argv.extend(['--limit', str(self.run_args['limit'])])
        if self.run_args['once']:
            argv.append('--once')
        if self.run_args['serial']:
            argv.append('--serial')
        if self.run_args['rchunks']:
            rchunks = []
            for d in self.run_args['rchunks']:
                c = '{},{}'.format(d, self.run_args['rchunks'][d])
                rchunks.extend(['--read_chunk', c])
            argv.extend(rchunks)
        if self.run_args['wchunks']:
            wchunks = []
            for d in self.run_args['wchunks']:
                c = '{},{}'.format(d, self.run_args['wchunks'][d])
                wchunks.extend(['--write_chunk', c])
            argv.extend(wchunks)
        argv.append('input.s2s')
        return argv

    def shortargs(self):
        long_to_short = {
            '--verbosity': '-v',
            '--write_mode': '-m',
            '--limit': '-l',
            '--once': '-l',
            '--serial': '-s',
            '--read_chunk': '-c',
            '--write_chunk': '-w',
        }
        return [long_to_short[a] if a in long_to_short else a for a in self.longargs()]

    def cliassert(self, args):
        opts, specfile = s2srun.cli(args)
        assert opts.once == self.run_args['once'], 'Once-file incorrect'
        assert opts.rchunks == self.run_args['rchunks'], 'Read chunks incorrect'
        assert opts.wchunks == self.run_args['wchunks'], 'Write chunks incorrect'
        assert opts.limit == self.run_args['limit'], 'Output limit incorrect'
        assert opts.write_mode == self.run_args['write_mode'], 'Write mode incorrect'
        assert opts.serial == self.run_args['serial'], 'Serial mode incorrect'
        assert opts.verbosity == self.run_args['verbosity'], 'Verbosity incorrect'
        assert specfile == 'input.s2s', 'Specfile name incorrect'

    def test_empty(self):
        with pytest.raises(ValueError):
            s2srun.cli([])

    def test_help(self):
        with pytest.raises(SystemExit):
            s2srun.cli(['-h'])

    def test_short(self):
        self.cliassert(self.shortargs())

    def test_long(self):
        self.cliassert(self.longargs())


class NetCDF4Tests(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def init(self):
        self.spec_args = {
            'infiles': config.slices,
            'ncfmt': 'netcdf4',
            'compression': 0,
            'prefix': 'out.',
            'suffix': '.nc',
            'timeseries': None,
            'metadata': [v for v in config.tvmvars] + ['time'],
            'meta1d': False,
            'backend': 'netCDF4',
        }

        self.run_args = {
            'serial': False,
            'chunks': [],
            'limit': 0,
            'verbosity': 1,
            'write_mode': 'w',
            'once': False,
            'specfile': 'input.s2s',
        }

        self.clean()
        make.generate_data()

    def tearDown(self):
        self.clean()

    def clean(self):
        for ncfile in glob('*.nc'):
            remove(ncfile)

    def check(self, tsvar):
        args = {}
        args.update(self.spec_args)
        args.update(self.run_args)
        assertions_dict = check_outfile(tsvar=tsvar, **args)
        failed_assertions = [key for key, value in assertions_dict.items() if value is False]
        assert_msgs = ['Output file check for variable {0!r}:'.format(tsvar)]
        assert_msgs.extend(['   {0}'.format(assrt) for assrt in failed_assertions])
        assert len(failed_assertions) == 0, eol.join(assert_msgs)

    def runargs(self):
        argv = [
            '-v',
            str(self.run_args['verbosity']),
            '-m',
            self.run_args['write_mode'],
            '-l',
            str(self.run_args['limit']),
        ]
        if self.run_args['once']:
            argv.append('-1')
        if self.run_args['serial']:
            argv.append('-s')
        if len(self.run_args['chunks']) > 0:
            chunks = []
            for c in self.run_args['chunks']:
                chunks.extend(['-c', c])
            argv.extend(chunks)
        argv.append('input.s2s')
        return argv

    def convert(self):
        Specifier(**self.spec_args).write('input.s2s')
        s2srun.main(self.runargs())
        remove('input.s2s')

    def test_defaults(self):
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_I1(self):
        self.spec_args['infiles'] = config.slices[1:2]
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_V0(self):
        self.run_args['verbosity'] = 0
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_V3(self):
        self.run_args['verbosity'] = 3
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_TSV2(self):
        self.spec_args['timeseries'] = config.tsvars[1:3] + ['tsvarX']
        self.convert()
        for tsvar in self.spec_args['timeseries']:
            if tsvar in config.tsvars:
                self.check(tsvar)
            else:
                fname = self.spec_args['prefix'] + tsvar + self.spec_args['suffix']
                assert not exists(fname), 'File {0!r} should not exist'.format(fname)

    def test_NC3(self):
        self.spec_args['ncfmt'] = 'netcdf'
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_ser(self):
        self.run_args['serial'] = True
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_CL1(self):
        self.spec_args['compression'] = 1
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_CL1_LSF3(self):
        self.spec_args['least_significant_digit'] = 3
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_meta1d(self):
        self.spec_args['meta1d'] = True
        self.spec_args['metadata'] = [v for v in config.tvmvars]
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_once(self):
        self.run_args['once'] = True
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_overwrite(self):
        self.run_args['write_mode'] = 'o'
        self.run_args['verbosity'] = 0
        self.convert()
        self.run_args['verbosity'] = 1
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_skip(self):
        self.run_args['write_mode'] = 's'
        self.run_args['verbosity'] = 0
        self.convert()
        self.run_args['verbosity'] = 1
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_append(self):
        self.run_args['write_mode'] = 'a'
        self.run_args['write_mode'] = 'w'
        self.spec_args['infiles'] = config.slices[0:2]
        self.convert()
        self.run_args['write_mode'] = 'a'
        self.spec_args['infiles'] = config.slices[2:]
        self.convert()
        for tsvar in config.tsvars:
            self.check(tsvar)

    def test_append_missing(self):
        missing = config.tsvars[2]
        self.run_args['write_mode'] = 'a'
        self.run_args['write_mode'] = 'w'
        self.spec_args['infiles'] = config.slices[0:2]
        self.convert()
        remove(self.spec_args['prefix'] + missing + self.spec_args['suffix'])
        self.run_args['write_mode'] = 'a'
        self.spec_args['infiles'] = config.slices[2:]
        self.convert()
        for tsvar in config.tsvars:
            if tsvar == missing:
                self.spec_args['infiles'] = config.slices[2:]
            else:
                self.spec_args['infiles'] = config.slices
            self.check(tsvar)


class NioTests(NetCDF4Tests):
    def setUp(self):
        NetCDF4Tests.setUp(self)
        self.spec_args['backend'] = 'Nio'

    def tearDown(self):
        self.clean()
