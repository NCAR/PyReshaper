#!/usr/bin/env python
#======================================================================
#
#  This script is designed to run cprnc netCDF file comparisons for
#  a given set of PyReshaper tests.
#
#======================================================================

import os
import sys
import glob
import string
import optparse
import json

from subprocess import Popen, PIPE, STDOUT


#==============================================================================
# Command-Line Interface Definition
#==============================================================================
def parse_cli():
    usage = 'usage: %prog [options] test_name1 test_name2 ...'
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-a', '--all', default=False,
                      action='store_true', dest='all',
                      help='True or False, indicating whether to run all tests '
                           '[Default: False]')
    parser.add_option('-i', '--testing_database', default=None,
                      help='Location of the testinfo.json file '
                           '[Default: None]')
    parser.add_option('-l', '--list', default=False,
                      action='store_true', dest='list_tests',
                      help='True or False, indicating whether to list all tests '
                           'that have been run with resulting output, instead of '
                           'actually comparing any tests. '
                           '[Default: False]')
    parser.add_option('-s', '--serial', default=False,
                      action='store_true', dest='serial',
                      help='True or False, indicating whether to check tests '
                           'in serial (True), rather than parallel (False). '
                           '[Default: False]')

    # Parse the CLI options and arguments
    (options, arguments) = parser.parse_args()

    # Set serial if listing only
    if options.list_tests:
        options.serial = True

    #  and return
    return (options, arguments)


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
# Get the testing database info
#==============================================================================
def get_testing_database(options):

    # Get the testinfo.json data
    testing_database_filename = ''
    if (options.testing_database == None):
        runtest_dir = os.path.dirname(__file__)
        testing_database_filename = os.path.join(runtest_dir, 'testinfo.json')
    else:
        testing_database_filename = os.path.abspath(options.testing_database)

    # Try opening and reading the testinfo file
    testing_database = {}
    try:
        testing_database_file = open(testing_database_filename, 'r')
        testing_database = dict(json.load(testing_database_file))
        testing_database_file.close()
    except:
        err_msg = 'Problem reading and parsing test info file: ' \
            + str(testing_database_filename)
        raise ValueError(err_msg)

    return testing_database


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
    if (options.all):
        possible_test_dirs = glob.glob(os.path.join('*', 'ser', '*'))
        possible_test_dirs.extend(glob.glob(os.path.join('*', 'par*', '*')))
    else:
        possible_test_dirs = arguments

    # Validate each possible test
    for possible_test_dir in possible_test_dirs:
        if (comm.rank == 0):
            print 'Validating possible test dir:', possible_test_dir
        temp, ncformat = os.path.split(possible_test_dir)
        temp, run_type = os.path.split(temp)
        temp, test_name = os.path.split(temp)
        if (comm.rank == 0):
            print '  Test Name:', test_name
            print '  Run Type:', run_type
            print '  NetCDF Format:', ncformat

        good_test = (test_name in testing_database)

        if (good_test):
            if (comm.rank == 0):
                print '  Test found in test info'
            new_results_dir = os.path.join(
                cwd, test_name, run_type, ncformat, 'output')
            good_test = good_test and os.path.isdir(new_results_dir)
            if (comm.rank == 0):
                print '  New results dir found:', new_results_dir

        if (good_test):
            os.chdir(new_results_dir)
            new_results_ls = glob.glob('*.nc')
            os.chdir(cwd)
            good_test = good_test and (len(new_results_ls) > 0)
            if (comm.rank == 0):
                print '  ' + str(len(new_results_ls)) + ' files found in directory'

        if (good_test):
            old_results_dir = testing_database[test_name]['results_dir']
            good_test = good_test and os.path.isdir(old_results_dir)
            if (comm.rank == 0):
                print '  Old results dir found:', old_results_dir

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

        if (good_test):
            cprnc_out_dir = os.path.join(
                cwd, test_name, run_type, ncformat, 'compare')
            if (comm.rank == 0):
                print '  CPRNC output directory will be:', cprnc_out_dir

            full_test_name = os.path.join(test_name, run_type, ncformat)
            if (comm.rank == 0):
                print '  Full test name:', full_test_name

            log_dir = os.path.join(cwd, full_test_name)
            log_filename = os.path.join(log_dir, 'check-' + test_name + '.log')
            if (comm.rank == 0):
                print '  Log file will be:', log_filename

            comparison_info[full_test_name] = {}
            comparison_info[full_test_name][
                'results_filenames'] = new_results_ls
            comparison_info[full_test_name][
                'new_results_dirs'] = new_results_dir
            comparison_info[full_test_name][
                'old_results_dirs'] = old_results_dir
            comparison_info[full_test_name]['cprnc_out_dirs'] = cprnc_out_dir
            comparison_info[full_test_name]['log_filenames'] = log_filename
            comparison_info[full_test_name]['log_output'] = []

        if (comm.rank == 0):
            print

    # Print out tests to be checked
    if (comm.rank == 0):
        print 'Valid tests for checking are:'
        for full_test_name in comparison_info.keys():
            print '  ', full_test_name
        print

    # Exit now, if only listing tests
    if (options.list_tests):
        sys.exit(0)

    # Reorganize the filenames to allow top-level looping over individual files
    # and initialize comparison statistics
    for full_test_name in comparison_info.keys():
        comparison_info[full_test_name]['num_failures'] = 0
        comparison_info[full_test_name]['num_checks'] = len(
            comparison_info[full_test_name]['results_filenames'])

    return comparison_info


#==============================================================================
# Run comparison tests
#==============================================================================
def compare_results(comparison_info, comm):

    # Base arguments to running cprnc
    cprnc_args = ["/glade/p/work/kpaul/installs/intel/12.1.5/cprnc/bin/cprnc",
                  "-m", "", ""]

    # String indicating IDENTITY (valid comparison)
    ident_str = "files seem to be IDENTICAL"

    # Before checking files, create CPRNC output directories for each test
    for full_test_name in comparison_info.keys():
        cprnc_out_dir = comparison_info[full_test_name]['cprnc_out_dirs']
        if (comm.rank == 0):
            if (not os.path.isdir(cprnc_out_dir)):
                os.makedirs(cprnc_out_dir)

    comm.sync()

    # Create a flat list of files to check
    files_to_check = []
    for full_test_name in comparison_info.keys():
        for file_name in comparison_info[full_test_name]['results_filenames']:
            check_dict = {'file_name': file_name,
                          'full_test_name': full_test_name}
            files_to_check.append(check_dict)

    if (comm.rank == 0):
        print 'Checking files:'
        print

    # Loop over all files (designated for this processor)
    for file_to_check in files_to_check[comm.rank::comm.size]:

        full_test_name = file_to_check['full_test_name']
        file_name = file_to_check['file_name']

        new_results_dir = comparison_info[
            full_test_name]['new_results_dirs']
        old_results_dir = comparison_info[
            full_test_name]['old_results_dirs']
        cprnc_out_dir = comparison_info[
            full_test_name]['cprnc_out_dirs']

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

    # Gather all of the statitics and log files from all processors
    all_log_output = comm.gather(
        comparison_info[full_test_name]['log_output'])
    all_num_failures = comm.gather(
        comparison_info[full_test_name]['num_failures'])

    # Open and write all log files
    if (comm.rank == 0):
        for full_test_name in comparison_info.keys():
            log_file = open(
                comparison_info[full_test_name]['log_filenames'], 'w')
            for rank_log_output in all_log_output:
                log_strings = rank_log_output[full_test_name]
                for log_str in log_strings:
                    log_file.write(log_str)
            log_file.close()

    # Print diagnostic info to stdout
    if (comm.rank == 0):
        print
        print 'Number of File Comparison Failures:'
        print
        for full_test_name in comparison_info.keys():
            num_test_failures = 0
            for rank_num_failures in all_num_failures:
                num_test_failures += rank_num_failures[full_test_name]
            print '  ', full_test_name + ':', num_test_failures,
            print 'failures (of',
            print comparison_info[full_test_name]['num_checks'],
            print 'total comparisons)'


#==============================================================================
# Main Program
#==============================================================================
def main(options, arguments):
    comm = BasicComm(serial=options.serial)
    testing_database = get_testing_database(options)
    comparison_info = get_comparison_info(
        options, arguments, comm, testing_database)
    compare_results(comparison_info, comm)


if __name__ == '__main__':
    options, arguments = parse_cli()
    main(options, arguments)
