#!/usr/bin/env python
#==============================================================================
#
#  runtests.py
#
#  Command-line script to run the Yellowstone test suite (i.e., the "bakeoff").
#  This script can run any number of ways.
#
#==============================================================================

# Builtin Modules
import os
import sys
import shutil
import argparse
import cPickle as pickle

# Package Modules
from utilities import testtools as tt
from utilities import runtools as rt

#==============================================================================
# Command-Line Interface Definition
#==============================================================================
_DESC_ = """This program is designed to run yellowstone-specific
            tests of the PyReshaper.  Each named test (or all tests if
            the -a or --all option is used) will be given a run
            directory in the "rundirs" directory with the same
            name as the test itself.  The run script will be placed
            in this run directory, as will be placed the run output
            error file.  All output data files will be placed in the
            output subdirectory."""

_PARSER_ = argparse.ArgumentParser(description=_DESC_)
_PARSER_.add_argument('-a', '--all', default=False,
                      action='store_true', dest='all_tests',
                      help='True or False, indicating whether to run all '
                           'tests [Default: False]')
_PARSER_.add_argument('-c', '--code', default='STDD0002', type=str,
                      help='The name of the project code for charging in '
                           'parallel runs (ignored if running in serial) '
                           '[Default: STDD0002]')
_PARSER_.add_argument('-i', '--infofile', default=None, type=str,
                      help='Location of the testinfo.json database file '
                           '[Default: None]')
_PARSER_.add_argument('-f', '--format', default='netcdf4c',
                      type=str, dest='ncformat',
                      help='The NetCDF file format to use for the output '
                           'data produced by the test.  [Default: netcdf4c]')
_PARSER_.add_argument('-l', '--list', default=False,
                      action='store_true', dest='list_tests',
                      help='True or False, indicating whether to list all '
                           'tests, instead of running tests. [Default: False]')
_PARSER_.add_argument('-m', '--multiple', default=False,
                      action='store_true', dest='multispec',
                      help='True or False, indications whether the tests '
                           'should be run from a single Reshaper submission '
                           '(i.e., multiple Specifiers in one run) '
                           '[Default: False]')
_PARSER_.add_argument('-n', '--nodes', default=0, type=int,
                      help='The integer number of nodes to request in parallel'
                           ' runs (0 means run in serial) [Default: 0]')
_PARSER_.add_argument('-o', '--overwrite', default=False,
                      action='store_true', dest='overwrite',
                      help='True or False, indicating whether to force '
                           'deleting any existing test or run directories, '
                           'if found [Default: False]')
_PARSER_.add_argument('-q', '--queue', default='economy', type=str,
                      help='The name of the queue to request in parallel runs '
                           '(ignored if running in serial) '
                           '[Default: economy]')
_PARSER_.add_argument('-s', '--skip_existing', default=False,
                      action='store_true',
                      help='Whether to skip time-series generation for '
                           'variables with existing output files. '
                           '[Default: False]')
_PARSER_.add_argument('--statsfile', default=None, type=str,
                      help='Location of the teststats.json database file '
                           '[Default: None]')
_PARSER_.add_argument('-t', '--tiling', default=16, type=int,
                      help='The integer number of processes per node to '
                           'request in parallel runs (ignored if running '
                           'in serial) [Default: 16]')
_PARSER_.add_argument('-w', '--wtime', default=240, type=int,
                      help='The number of minutes to request for the wall '
                           'clock in parallel runs (ignored if running in '
                           'serial) [Default: 240]')
_PARSER_.add_argument('-z', '--analyze', default=False,
                      action='store_true',
                      help='Whether to analyze (generate statistics for) the '
                           'test input, rather than run the tests. '
                           '[Default: False]')
_PARSER_.add_argument('test', type=str, nargs='*',
                      help='Name of test to run or analyze.')


#==============================================================================
# Main Function
#==============================================================================
def runtests(args):
    """
    Run a set of tests

    Parameters:
        args (argparse.Namespace): A namespace of command-line parsed arguments
            describing the tests to run and how to run them
    """

    # Check for tests to run
    if len(args.test) == 0 and not args.all_tests and not args.list_tests:
        _PARSER_.print_help()
        sys.exit(1)

    # Create/read the testing info and stats files
    testdb = tt.TestDB(dbname=args.infofile, stname=args.statsfile)

    # List tests if only listing
    if args.list_tests:
        testdb.print_tests()
        sys.exit(1)

    # Generate the list of tests to run/analyze
    if args.all_tests:
        test_list = testdb.get_database().keys()
    else:
        test_list = [t for t in args.test if t in testdb.get_database()]

    # Analyze test input, if requested (overwrite forces re-analysis)
    if args.analyze:
        testdb.analyze(tests=test_list, force=args.overwrite)
        testdb.save_statistics(stname=args.statsfile)
        sys.exit(0)

    # Run the requested tests
    else:
        for test_name in test_list:

            print 'Running test: {0!s}'.format(test_name)

            # Set the test directory
            if args.nodes > 0:
                runtype = 'par{0!s}x{1!s}'.format(args.nodes, args.tiling)
            else:
                runtype = 'ser'
            testdir = os.path.abspath(os.path.join('rundirs', str(test_name), runtype))

            # If the test directory doesn't exist, make it and move into it
            if os.path.exists(testdir):
                if args.overwrite:
                    shutil.rmtree(testdir)
                else:
                    print "   Skipping existing"
                    continue
            if not os.path.exists(testdir):
                os.makedirs(testdir)
            os.chdir(testdir)

            # Set the output directory
            outputdir = os.path.join(testdir, 'output')

            # If the output directory doesn't exists, create it
            if not os.path.exists(outputdir):
                os.mkdir(outputdir)

            # Create the specifier and write to file (specfile)
            testspec = testdb.create_specifier(test_name=str(test_name),
                                               ncfmt=args.ncformat,
                                               outdir=outputdir)
            testspecfile = str(test_name) + '.spec'
            pickle.dump(testspec, open(testspecfile, 'wb'))

            # Generate the command and arguments
            if args.nodes > 0:
                runcmdargs = ['poe', 'slice2series', '-v3']
            else:
                runcmdargs = ['slice2series', '--serial', '-v3']
            runcmdargs.extend(['--specfile', testspecfile])
            runcmd = ' '.join(runcmdargs)

            # Create and start the job
            job = rt.Job(runcmds=[runcmd], nodes=args.nodes, name=str(test_name),
                         tiling=args.tiling, minutes=args.wtime,
                         queue=args.queue, project=args.code)
            job.start()
            if args.nodes <= 0:
                job.wait()


#==============================================================================
# Main Command-line Operation
#==============================================================================
if __name__ == '__main__':
    args = _PARSER_.parse_args()
    runtests(args)
