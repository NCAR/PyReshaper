#!/usr/bin/env python
#======================================================================
#
#  This script is designed to run cprnc netCDF file comparisons for
#  a given set of PyReshaper tests.
#
#======================================================================

import os
import re
import sys
import glob
import string
import argparse
import json
from subprocess import Popen, PIPE, STDOUT

# Package Modules
from utilities import testtools as tt
from utilities import runtools as rt

#==============================================================================
# Command-Line Interface Definition
#==============================================================================
_DESC_ = 'Check the results of tests found in the rundirs directory.'
_PARSER_ = argparse.ArgumentParser(description=_DESC_)
_PARSER_.add_argument('-c', '--code', default='STDD0002', type=str,
                      help='The name of the project code for charging in '
                           'parallel runs (ignored if running in serial) '
                           '[Default: STDD0002]')
_PARSER_.add_argument('-i', '--infofile', default=None,
                      help='Location of the testinfo.json file '
                           '[Default: None]')
_PARSER_.add_argument('-l', '--list', default=False,
                      action='store_true', dest='list_tests',
                      help='True or False, indicating whether to list all '
                           'tests that have been run with resulting output, '
                           'instead of actually comparing any tests. '
                           '[Default: False]')
_PARSER_.add_argument('-n', '--nodes', default=0, type=int,
                      help='The integer number of nodes to request in parallel'
                           ' runs (0 means run in serial) [Default: 0]')
_PARSER_.add_argument('-q', '--queue', default='economy', type=str,
                      help='The name of the queue to request in parallel runs '
                           '(ignored if running in serial) '
                           '[Default: economy]')
_PARSER_.add_argument('-w', '--wtime', default=240, type=int,
                      help='The number of minutes to request for the wall '
                           'clock in parallel runs (ignored if running in '
                           'serial) [Default: 240]')
_PARSER_.add_argument('-x', '--executable', type=str,
                      default='/glade/p/work/kpaul/installs/intel/12.1.5/cprnc/bin/cprnc',
                      help='The path to the CPRNC executable.')
_PARSER_.add_argument('testdir', type=str, nargs='*',
                      help='Name of a test directory to check')


#==============================================================================
# grep - Search through a file for a given pattern
#==============================================================================
def grep(pattern, filename):
    if not os.path.exists(filename):
        return None
    fobj = open(filename, 'r')
    results = [line.rstrip() for line in fobj if re.search(pattern, line)]
    print results
    fobj.close()
    if len(results) == 0:
        return None
    else:
        return results


#==============================================================================
# Local MPI options and handling
#==============================================================================
class BasicComm(object):

    def __init__(self, serial=False):
        self.serial = serial
        self.rank = 0
        self.size = 1
        if not self.serial:
            from mpi4py import MPI
            self.MPI = MPI
            self.rank = MPI.COMM_WORLD.Get_rank()
            self.size = MPI.COMM_WORLD.Get_size()

    def sync(self):
        if self.serial:
            return
        else:
            self.MPI.COMM_WORLD.Barrier()

    def gather(self, data):
        if self.serial:
            return data
        else:
            return self.MPI.COMM_WORLD.gather(data)


#==============================================================================
# Parse input and set up comparison information
#==============================================================================
def get_comparison_info(options, arguments, comm, testing_database):

    # Current working directory
    cwd = os.getcwd()

    # Comparison info - Organized by full_test_name
    comparison_info = {}

    # Parse each test argument for test directories
    possible_test_dirs = []
    if (options.all_tests):
        possible_test_dirs = glob.glob(
            os.path.join('results', '*', 'ser', '*'))
        possible_test_dirs.extend(
            glob.glob(os.path.join('results', '*', 'par*', '*')))
    else:
        possible_test_dirs = arguments

    # Validate each possible test
    for possible_test_dir in possible_test_dirs:
        if (comm.rank == 0):
            print 'Validating possible test dir:', possible_test_dir

        # Split out the NetCDF format, run type (serial, parallel), and
        # test name from the directory name
        root, ncformat = os.path.split(possible_test_dir)
        root, run_type = os.path.split(root)
        root, test_name = os.path.split(root)
        if (comm.rank == 0):
            print '  Test Name:', test_name
            print '  Run Type:', run_type
            print '  NetCDF Format:', ncformat

        # Check that the test name is in the database
        # Define this as a "Good" test
        good_test = (test_name in testing_database)

        # If still "good", check for results/output directory
        if (good_test):
            if (comm.rank == 0):
                print '  Test found in test info'
            new_results_dir = os.path.join(
                cwd, 'results', test_name, run_type, ncformat, 'output')
            good_test = good_test and os.path.isdir(new_results_dir)
            if (comm.rank == 0):
                print '  New results dir found:', new_results_dir

        # If still "good", check for >0 output/result file
        if (good_test):
            os.chdir(new_results_dir)
            new_results_ls = glob.glob('*.nc')
            os.chdir(cwd)
            good_test = good_test and (len(new_results_ls) > 0)
            if (comm.rank == 0):
                print '  ' + str(len(new_results_ls)) + ' files found in directory'

        # If still "good", check for results dir for comparison
        if (good_test):
            old_results_dir = testing_database[test_name]['results_dir']
            good_test = good_test and os.path.isdir(old_results_dir)
            if (comm.rank == 0):
                print '  Old results dir found:', old_results_dir

        # If still "good", look for missing files and make sure there are
        # some new files to compare against the old
        if (good_test):
            os.chdir(old_results_dir)
            old_results_ls = glob.glob('*.nc')
            os.chdir(cwd)
            missing_tests = set(new_results_ls) - set(old_results_ls)
            if (len(missing_tests) > 0):
                if (comm.rank == 0):
                    print '  Did not find', len(missing_tests),
                    print 'new test files in old results dir:'
                for missing_test in missing_tests:
                    new_results_ls.remove(missing_test)
                    if (comm.rank == 0):
                        print '    ', missing_test
                if (comm.rank == 0):
                    print '  Missing tests have been removed from comparison list'
                good_test = good_test and (len(new_results_ls) > 0)
            else:
                if (comm.rank == 0):
                    print '  All new test files found in old results dir'

        # If still "good", generate the name of the log file and the directory
        # in which comparison (CPRNC) results will be placed
        if (good_test):
            full_test_name = os.path.join(test_name, run_type, ncformat)
            if (comm.rank == 0):
                print '  Full test name:', full_test_name

            log_dir = os.path.join(cwd, 'results', full_test_name)
            log_filename = os.path.join(log_dir, 'check-' + test_name + '.log')
            if (comm.rank == 0):
                print '  Log file will be:', log_filename

            cprnc_out_dir = os.path.join(log_dir, 'compare')
            if (comm.rank == 0):
                print '  CPRNC output directory will be:', cprnc_out_dir

            comparison_info[full_test_name] = {}
            comparison_info[full_test_name][
                'results_filenames'] = new_results_ls
            comparison_info[full_test_name][
                'new_results_dir'] = new_results_dir
            comparison_info[full_test_name][
                'old_results_dir'] = old_results_dir
            comparison_info[full_test_name]['cprnc_out_dir'] = cprnc_out_dir
            comparison_info[full_test_name]['log_filename'] = log_filename
            comparison_info[full_test_name]['log_output'] = []

        if (comm.rank == 0):
            print

    # Print out tests to be checked
    if (comm.rank == 0):
        if len(comparison_info.keys()) == 0:
            print 'No valid tests found.'
            sys.exit(0)
        else:
            print 'Valid tests for checking are:'
            total_num_files = 0
            for full_test_name in comparison_info.keys():
                num_files = len(
                    comparison_info[full_test_name]['results_filenames'])
                total_num_files += num_files
                print '  ', full_test_name, ' (' + str(num_files) + ' files)'
            print 'Total number of file comparisons:', total_num_files
        print

    # Exit now, if only listing tests
    if (options.list_tests):
        sys.exit(0)

    # Add slots for counting failures and total results
    for full_test_name in comparison_info.keys():
        comparison_info[full_test_name]['num_failures'] = 0
        comparison_info[full_test_name]['num_checks'] = len(
            comparison_info[full_test_name]['results_filenames'])

    return comparison_info


#==============================================================================
# Run comparison tests
#==============================================================================
def compare_results(comparison_info, comm, cprnc_exec):

    # Base arguments to running cprnc
    cprnc_args = [cprnc_exec, "-m", "", ""]

    # String indicating IDENTITY (valid comparison)
    ident_str = "files seem to be IDENTICAL"

    # Before checking files, create CPRNC output directories for each test
    for full_test_name in comparison_info.keys():
        cprnc_out_dir = comparison_info[full_test_name]['cprnc_out_dir']
        if (comm.rank == 0):
            if (not os.path.isdir(cprnc_out_dir)):
                os.makedirs(cprnc_out_dir)

    # Wait for for master rank to create the output directory
    comm.sync()

    # Create a flat list of files to check
    files_to_check = []
    for full_test_name in comparison_info.keys():
        for file_name in comparison_info[full_test_name]['results_filenames']:
            check_dict = {'file_name': file_name,
                          'full_test_name': full_test_name}
            files_to_check.append(check_dict)

    # Start comparing output files
    if (comm.rank == 0):
        print 'Checking files:'
        print

    # Loop over all files (designated for this processor)
    for file_to_check in files_to_check[comm.rank::comm.size]:

        full_test_name = file_to_check['full_test_name']
        file_name = file_to_check['file_name']

        new_results_dir = comparison_info[
            full_test_name]['new_results_dir']
        old_results_dir = comparison_info[
            full_test_name]['old_results_dir']
        cprnc_out_dir = comparison_info[
            full_test_name]['cprnc_out_dir']

        # Create the old and new file paths
        new_results_path = os.path.join(new_results_dir, file_name)
        old_results_path = os.path.join(old_results_dir, file_name)

        # CPRNC output file
        cprnc_out_filename = file_name + '.cprnc'
        cprnc_out_path = os.path.join(cprnc_out_dir, cprnc_out_filename)

        # If the output file already exists, read it, otherwise write it
        cprnc_out = ''
        if (os.path.exists(cprnc_out_path)):
            cprnc_out_file = open(cprnc_out_path, 'r')
            cprnc_out = cprnc_out_file.read()
            cprnc_out_file.close()
        else:
            # Set the arguments to CPRNC, run, and catch output
            cprnc_args[2:] = [old_results_path, new_results_path]
            cprnc_proc = Popen(cprnc_args, stdout=PIPE, stderr=STDOUT)
            cprnc_out = cprnc_proc.communicate()[0]

            # Write the output file
            cprnc_out_file = open(cprnc_out_path, 'w')
            cprnc_out_file.write(cprnc_out)
            cprnc_out_file.close()

        # Add to the list of log file output
        log_result = ' GOOD: '
        if (string.rfind(cprnc_out, ident_str) < 0):
            comparison_info[full_test_name]['num_failures'] += 1
            log_result = '* BAD: '
        log_str = log_result + file_name + os.linesep
        comparison_info[full_test_name]['log_output'].append(log_str)

        print '   ', log_result + file_name

    # Now wait for all files to be processed
    comm.sync()

    # Assemble rank info to gather
    rank_log_output = {}
    rank_num_failures = {}
    for full_test_name in comparison_info.keys():
        rank_log_output[full_test_name] = comparison_info[
            full_test_name]['log_output']
        rank_num_failures[full_test_name] = comparison_info[
            full_test_name]['num_failures']

    # Gather all of the statitics and log files from all processors
    all_log_output = comm.gather(rank_log_output)
    all_num_failures = comm.gather(rank_num_failures)

    if comm.rank == 0:

        # Open and write all log files
        print
        print 'Writing log files.'
        print
        for full_test_name in comparison_info.keys():
            log_file = open(
                comparison_info[full_test_name]['log_filename'], 'w')
            for r_log_output in all_log_output:
                log_strings = r_log_output[full_test_name]
                for log_str in log_strings:
                    log_file.write(log_str)
            log_file.close()

        # Print diagnostic info to stdout
        print
        print 'Number of File Comparison Failures:'
        print
        for full_test_name in comparison_info.keys():
            num_test_failures = sum([nf[full_test_name]
                                     for nf in all_num_failures])
            print '  ', full_test_name + ':', num_test_failures,
            print 'failures (of',
            print comparison_info[full_test_name]['num_checks'],
            print 'total comparisons)'

#==============================================================================
# Command-Line Operation
#==============================================================================
if __name__ == '__main__':
    args = _PARSER_.parse_args()

    # Create/read the testing info and stats files
    testdb = tt.TestDB(dbname=args.infofile)
    testdb_dict = testdb.get_database()

    # Get the list of the available test names for comparison
    valid_names = [str(t) for t in testdb_dict.keys()]

    # Find valid tests
    found_tests = {}
    for rdir in glob.glob(os.path.join('results.d', '*', '[ser,par]*', '*')):
        tempdir, ncfmt = os.path.split(rdir)
        tempdir, runtype = os.path.split(tempdir)
        tempdir, test_name = os.path.split(tempdir)
        if test_name in valid_names:

            print 'Found test results for name: {}'.format(test_name)
            print 'Results directory: {}'.format(rdir)

            # Look for the corresponding log file and check for successful run
            logfile = os.path.join(rdir, '{0!s}.log'.format(test_name))
            if not grep(r'Successfully completed.', logfile):
                continue
            print 'Test completed successfully.'

            # Look for the specfile
            specfile = os.path.join(rdir, '{0!s}.spec'.format(test_name))
            if not os.path.exists(specfile):
                continue
            print 'Found specfile: {}'.format(specfile)

            # Look for the output directory
            outdir = os.path.join(rdir, 'output')
            if not os.path.exists(outdir):
                continue
            print 'Found output directory: {}'.format(outdir)
