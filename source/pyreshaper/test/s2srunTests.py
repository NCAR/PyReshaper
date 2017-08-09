"""
Parallel Tests for the Reshaper class

Copyright 2017, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

import unittest

import imp
import sys
import inspect
from glob import glob
from cStringIO import StringIO
from os import linesep as eol
from os import remove
from os.path import exists
from mpi4py import MPI
import pickle

from pyreshaper.specification import Specifier

import mkTestData

s2srun = imp.load_source('s2srun', '../../../scripts/s2srun')

MPI_COMM_WORLD = MPI.COMM_WORLD  # @UndefinedVariable


#=======================================================================================================================
# CLITests
#=======================================================================================================================
class CLITests(unittest.TestCase):

    def test_empty(self):
        argv = []
        self.assertRaises(ValueError, s2srun.cli, argv)

    def test_help(self):
        argv = ['-h']
        self.assertRaises(SystemExit, s2srun.cli, argv)

    def test_defaults(self):
        argv = ['s2srunTests.py']
        opts, opts_specfile = s2srun.cli(argv)
        self.assertFalse(opts.once, 'Once-file is set')
        self.assertEqual(opts.limit, 0, 'Output limit is set')
        self.assertEqual(opts.write_mode, 'w', 'Write mode is not "w"')
        self.assertFalse(opts.serial, 'Serial mode is set')
        self.assertEqual(opts.verbosity, 1, 'Verbosity is not 1')
        self.assertEqual(opts_specfile, argv[0], 'Specfile name is not set')

    def test_short(self):
        once = True
        limit = 3
        write_mode = 'x'
        serial = True
        verbosity = 5
        specfile = 'myspec.s2s'

        argv = []
        if once:
            argv.append('-1')
        argv.extend(['-l', str(limit)])
        argv.extend(['-m', str(write_mode)])
        if serial:
            argv.append('-s')
        argv.extend(['-v', str(verbosity)])
        argv.append(specfile)
        opts, opts_specfile = s2srun.cli(argv)

        self.assertEqual(opts.once, once, 'Once-file is not {0!r}'.format(once))
        self.assertEqual(opts.limit, limit, 'Output limit is not {0!r}'.format(limit))
        self.assertEqual(opts.write_mode, write_mode, 'Write mode is not {0!r}'.format(write_mode))
        self.assertEqual(opts.serial, serial, 'Serial mode is not {0!r}'.format(serial))
        self.assertEqual(opts.verbosity, verbosity, 'Verbosity is not {0!r}'.format(verbosity))
        self.assertEqual(opts_specfile, specfile, 'Specfile name is not {0!r}'.format(specfile))

    def test_long(self):
        once = True
        limit = 3
        write_mode = 'x'
        serial = True
        verbosity = 5
        specfile = 'myspec.s2s'

        argv = []
        if once:
            argv.append('--once')
        argv.extend(['--limit', str(limit)])
        argv.extend(['--write_mode', str(write_mode)])
        if serial:
            argv.append('--serial')
        argv.extend(['--verbosity', str(verbosity)])
        argv.append(specfile)
        opts, opts_specfile = s2srun.cli(argv)

        self.assertEqual(opts.once, once, 'Once-file is not {0!r}'.format(once))
        self.assertEqual(opts.limit, limit, 'Output limit is not {0!r}'.format(limit))
        self.assertEqual(opts.write_mode, write_mode, 'Write mode is not {0!r}'.format(write_mode))
        self.assertEqual(opts.serial, serial, 'Serial mode is not {0!r}'.format(serial))
        self.assertEqual(opts.verbosity, verbosity, 'Verbosity is not {0!r}'.format(verbosity))
        self.assertEqual(opts_specfile, specfile, 'Specfile name is not {0!r}'.format(specfile))


#=======================================================================================================================
# MainTests
#=======================================================================================================================
class MainTests(unittest.TestCase):

    def setUp(self):

        # Parallel Management - Just for Tests
        self.rank = MPI_COMM_WORLD.Get_rank()
        self.size = MPI_COMM_WORLD.Get_size()
        
        # Default arguments for testing
        self.spec_args = {'infiles': mkTestData.slices,
                          'ncfmt': 'netcdf4',
                          'compression': 0,
                          'prefix': 'out.',
                          'suffix': '.nc',
                          'timeseries': None,
                          'metadata': [v for v in mkTestData.tvmvars] + ['time'],
                          'meta1d': False,
                          'backend': 'netCDF4'}
        self.create_args = {'serial': False,
                            'verbosity': 1,
                            'wmode': 'w',
                            'once': False,
                            'simplecomm': None}
        self.convert_args = {'output_limit': 0,
                             'chunks': None}
        
        # Test Data Generation
        self.clean()
        if self.rank == 0:
            mkTestData.generate_data()
        MPI_COMM_WORLD.Barrier()

    def tearDown(self):
        self.clean()

    def clean(self):
        if self.rank == 0:
            for ncfile in glob('*.nc'):
                remove(ncfile)
        MPI_COMM_WORLD.Barrier()

    def header(self, testname):
        if self.rank == 0:
            nf = len(self.spec_args['infiles'])
            nt = len(mkTestData.tsvars) if self.spec_args['timeseries'] is None else len(self.spec_args['timeseries'])

            hline = '-' * 70
            hdrstr = [hline, '{}:'.format(testname), '',
                      '   specifier({} infile(s), {} TSV(s), ncfmt={ncfmt},'.format(nf, nt, **self.spec_args),
                      '             compression={compression}, meta1d={meta1d}, backend={backend})'.format(**self.spec_args),
                      '   create(serial={serial}, verbosity={verbosity}, wmode={wmode},'.format(**self.create_args),
                      '          once={once}, simplecomm={simplecomm})'.format(**self.create_args),
                      '   convert(output_limit={output_limit}, chunks={chunks})'.format(**self.convert_args), hline]
            print eol.join(hdrstr)

    def check(self, tsvar):
        args = {}
        args.update(self.spec_args)
        args.update(self.create_args)
        args.update(self.convert_args)        
        assertions_dict = mkTestData.check_outfile(tsvar=tsvar, **args)
        failed_assertions = [key for key, value in assertions_dict.iteritems() if value is False]
        assert_msgs = ['Output file check for variable {0!r}:'.format(tsvar)]
        assert_msgs.extend(['   {0}'.format(assrt) for assrt in failed_assertions])
        self.assertEqual(len(failed_assertions), 0, eol.join(assert_msgs))

    def convert(self):
        if not (self.create_args.get('serial', False) and self.rank > 0):
            spec = Specifier(**self.spec_args)
            specfile = 'input.s2s'
            pickle.dump(spec, open(specfile, 'wb'))

            argv = ['-v', str(self.create_args.get('verbosity', 1)), '-m', self.create_args.get('wmode', 'w')]
            if self.create_args.get('once', False):
                argv.append('-1')
            argv.append(specfile)
            s2srun.main(argv)
            remove(specfile)
        MPI_COMM_WORLD.Barrier()

    def test_defaults(self):
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_I1(self):
        self.spec_args['infiles'] = mkTestData.slices[1:2]
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_V0(self):
        self.create_args['verbosity'] = 0
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_V3(self):
        self.create_args['verbosity'] = 3
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_TSV2(self):
        self.spec_args['timeseries'] = mkTestData.tsvars[1:3] + ['tsvarX']
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in self.spec_args['timeseries']:
                if tsvar in mkTestData.tsvars:
                    self.check(tsvar)
                else:
                    fname = self.spec_args['prefix'] + tsvar + self.spec_args['suffix']
                    self.assertFalse(exists(fname), 'File {0!r} should not exist'.format(fname))
        MPI_COMM_WORLD.Barrier()

    def test_NC3(self):
        self.spec_args['ncfmt'] = 'netcdf'
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_Nio(self):
        self.spec_args['backend'] = 'Nio'
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()
        
    def test_netCDF4(self):
        self.spec_args['backend'] = 'netCDF4'
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_ser(self):
        self.create_args['serial'] = True
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_CL1(self):
        self.spec_args['compression'] = 1
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_meta1d(self):
        self.spec_args['meta1d'] = True
        self.spec_args['metadata'] = [v for v in mkTestData.tvmvars]
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_once(self):
        self.create_args['once'] = True
        self.header(inspect.currentframe().f_code.co_name)
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_overwrite(self):
        self.create_args['wmode'] = 'o'
        self.header(inspect.currentframe().f_code.co_name)
        self.create_args['verbosity'] = 0
        self.convert()
        self.create_args['verbosity'] = 1
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_skip(self):
        self.create_args['wmode'] = 's'
        self.header(inspect.currentframe().f_code.co_name)
        self.create_args['verbosity'] = 0
        self.convert()
        self.create_args['verbosity'] = 1
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_append(self):
        self.create_args['wmode'] = 'a'
        self.header(inspect.currentframe().f_code.co_name)
        self.create_args['wmode'] = 'w'
        self.spec_args['infiles'] = mkTestData.slices[0:2]
        self.convert()
        self.create_args['wmode'] = 'a'
        self.spec_args['infiles'] = mkTestData.slices[2:]
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()

    def test_append_missing(self):
        missing = mkTestData.tsvars[2]
        self.create_args['wmode'] = 'a'
        self.header(inspect.currentframe().f_code.co_name)
        self.create_args['wmode'] = 'w'
        self.spec_args['infiles'] = mkTestData.slices[0:2]
        self.convert()
        if self.rank == 0:
            remove(self.spec_args['prefix'] + missing + self.spec_args['suffix'])
        MPI_COMM_WORLD.Barrier()
        self.create_args['wmode'] = 'a'
        self.spec_args['infiles'] = mkTestData.slices[2:]
        self.convert()
        if self.rank == 0:
            for tsvar in mkTestData.tsvars:
                if tsvar == missing:
                    self.spec_args['infiles'] = mkTestData.slices[2:]
                else:
                    self.spec_args['infiles'] = mkTestData.slices
                self.check(tsvar)
        MPI_COMM_WORLD.Barrier()


#=======================================================================================================================
# CLI
#=======================================================================================================================
if __name__ == "__main__":
    hline = '=' * 70
    if MPI_COMM_WORLD.Get_rank() == 0:
        print hline
        print 'STANDARD OUTPUT FROM ALL TESTS:'
        print hline
    MPI_COMM_WORLD.Barrier()

    if MPI_COMM_WORLD.Get_rank() == 0:
        mystream = StringIO()
        clitests = unittest.TestLoader().loadTestsFromTestCase(CLITests)
        unittest.TextTestRunner(stream=mystream).run(clitests)
        print hline
        print 'CLI TESTS RESULTS:'
        print hline
        print mystream.getvalue()
    MPI_COMM_WORLD.Barrier()

    mystream = StringIO()
    maintests = unittest.TestLoader().loadTestsFromTestCase(MainTests)
    unittest.TextTestRunner(stream=mystream).run(maintests)
    MPI_COMM_WORLD.Barrier()

    results = MPI_COMM_WORLD.gather(mystream.getvalue())
    if MPI_COMM_WORLD.Get_rank() == 0:
        for rank, result in enumerate(results):
            print hline
            print 'MAIN TESTS RESULTS FOR RANK ' + str(rank) + ':'
            print hline
            print str(result)
