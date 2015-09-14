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
import argparse

# Package Modules
from utilities import testtools as tt
from utilities import runtools as rt


#==============================================================================
# Command-Line Interface Definition
#==============================================================================
def parse_cli():
    desc = """This program is designed to run yellowstone-specific
              tests of the PyReshaper.  Each named test (or all tests if
              the -a or --all option is used) will be given a run
              directory in the "rundirs" directory with the same
              name as the test itself.  The run script will be placed
              in this run directory, as will be placed the run output
              error file.  All output data files will be placed in the
              output subdirectory."""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-a', '--all', default=False,
                        action='store_true', dest='all',
                        help='True or False, indicating whether to run all tests '
                             '[Default: False]')
    parser.add_argument('-c', '--code', default='STDD0002',
                        help='The name of the project code for charging in '
                             'parallel runs (ignored if running in serial) '
                             '[Default: STDD0002]')
    parser.add_argument('-d', '--database', default=None,
                        help='Location of the testinfo.json file '
                             '[Default: None]')
    parser.add_argument('-f', '--format', default='netcdf4c', dest='ncformat',
                        help='The NetCDF file format to use for the output data '
                             'produced by the test.  [Default: netcdf4c]')
    parser.add_argument('-l', '--list', default=False,
                        action='store_true', dest='list_tests',
                        help='True or False, indicating whether to list all tests, '
                             'instead of running tests. [Default: False]')
    parser.add_argument('-m', '--multiple', default=False,
                        action='store_true', dest='multiple',
                        help='True or False, indicating whether to run the tests '
                             'in a single MultiSpecReshaper instance, which runs '
                             'each parallel test in sequence in the same job, '
                             'instead of running each parallel test individually '
                             'with its own job. [Default: False]')
    parser.add_argument('-o', '--overwrite', default=False,
                        action='store_true', dest='overwrite',
                        help='True or False, indicating whether to force deleting '
                             'any existing test or run directories, if found '
                             '[Default: False]')
    parser.add_argument('-q', '--queue', default='economy',
                        help='The name of the queue to request in parallel runs '
                             '(ignored if running in serial) '
                             '[Default: economy]')
    parser.add_argument('-n', '--nodes', default=0, type='int',
                        help='The integer number of nodes to request in parallel'
                             ' runs (0 means run in serial) [Default: 0]')
    parser.add_argument('-S', '--skip_existing', default=False,
                        action='store_true', dest='skip_existing',
                        help='Whether to skip time-series generation for '
                             'variables with existing output files. '
                             '[Default: False]')
    parser.add_argument('-t', '--tiling', default=16, type='int',
                        help='The integer number of processes per node to request '
                             'in parallel runs (ignored if running in serial) '
                             '[Default: 16]')
    parser.add_argument('-w', '--wtime', default=240, type='int',
                        help='The number of minutes to request for the wall clock '
                             'in parallel runs (ignored if running in serial) '
                             '[Default: 240]')

    # Return the parsed CLI options
    return parser.parse_args()


#==============================================================================
# Main Function
#==============================================================================
def runtests(options={}, tests=[]):
    print options
    print tests


#==============================================================================
# Main Command-line Operation
#==============================================================================
if __name__ == '__main__':
    options, arguments = parse_cli()
    print options
    print arguments
    runtests(dict(options), tests=list(arguments))
