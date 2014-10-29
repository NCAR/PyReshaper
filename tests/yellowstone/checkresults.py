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
usage = 'usage: %prog [options] test_name1 test_name2 ...'
parser = optparse.OptionParser(usage=usage)
parser.add_option('-a', '--all', default=False,
                  action='store_true', dest='all',
                  help='True or False, indicating whether to run all tests '
                       '[Default: False]')
parser.add_option('-i', '--test_info', default=None,
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

# Parse the CLI options and assemble the Reshaper inputs
(options, arguments) = parser.parse_args()

#==============================================================================
# Determine MPI options, so that STDOUT can be prevented from worker nodes
#==============================================================================

# Initialize MPI variables
mpi_rank = 0
mpi_size = 1

if (not options.serial):
    from mpi4py import MPI
    mpi_rank = MPI.COMM_WORLD.Get_rank()
    mpi_size = MPI.COMM_WORLD.Get_size()

def sync():
    if (options.serial):
        return
    else:
        MPI.COMM_WORLD.Barrier()

#==============================================================================
# Parse input and set up comparison information
#==============================================================================

# Get the testinfo.json data
test_info_filename = ''
if (options.test_info == None):
    runtest_dir = os.path.dirname(__file__)
    test_info_filename = os.path.join(runtest_dir, 'testinfo.json')
else:
    test_info_filename = os.path.abspath(options.test_info)

# Try opening and reading the testinfo file
test_info = {}
try:
    test_info_file = open(test_info_filename, 'r')
    test_info = dict(json.load(test_info_file))
    test_info_file.close()
except:
    err_msg = 'Problem reading and parsing test info file: ' \
            + str(test_info_filename)
    raise ValueError(err_msg)

# Current working directory
cwd = os.getcwd()

# Comparison info - Organized by full_test_name
full_test_names = []
new_results_dirs = {}
old_results_dirs = {}
results_filenames = {}
cprnc_out_dirs = {}
log_filenames = {}
log_output = {}

# Parse each test argument for test directories
possible_test_dirs = []
if (options.all):
    possible_test_dirs = glob.glob(os.path.join('*', 'ser'))
    possible_test_dirs.extend(glob.glob(os.path.join('*', 'par*')))
else:
    possible_test_dirs = arguments

# Validate each possible test
for possible_test_dir in possible_test_dirs:
    if (mpi_rank == 0):
        print 'Validating possible test dir:', possible_test_dir
    temp, run_type = os.path.split(possible_test_dir)
    temp, test_name = os.path.split(temp)
    if (mpi_rank == 0):
        print '  Test Name:', test_name
        print '  Run Type:', run_type

    good_test = (test_name in test_info)

    if (good_test):
        if (mpi_rank == 0):
            print '  Test found in test info'
        new_results_dir = os.path.join(cwd, test_name, run_type, 'output')
        good_test = good_test and os.path.isdir(new_results_dir)
        if (mpi_rank == 0):
            print '  New results dir found:', new_results_dir

    if (good_test):
        os.chdir(new_results_dir)
        new_results_ls = glob.glob('*.nc')
        os.chdir(cwd)
        good_test = good_test and (len(new_results_ls) > 0)
        if (mpi_rank == 0):
            print '  ' + str(len(new_results_ls)) + ' files found in directory'

    if (good_test):
        old_results_dir = test_info[test_name]['results_dir']
        good_test = good_test and os.path.isdir(old_results_dir)
        if (mpi_rank == 0):
            print '  Old results dir found:', old_results_dir

    if (good_test):
        os.chdir(old_results_dir)
        old_results_ls = glob.glob('*.nc')
        os.chdir(cwd)
        missing_tests = set(new_results_ls) - set(old_results_ls)
        if (len(missing_tests) > 0):
            if (mpi_rank == 0):
                print '  Did not find', len(missing_tests),
                print 'new test files in old results dir:'
            for missing_test in missing_tests:
                new_results_ls.remove(missing_test)
                if (mpi_rank == 0):
                    print '    ', missing_test
            if (mpi_rank == 0):
                print '  Missing tests have been removed from comparison list'
            good_test = good_test and (len(new_results_ls) > 0)
        else:
            if (mpi_rank == 0):
                print '  All new test files found in old results dir'

    if (good_test):
        cprnc_out_dir = os.path.join(cwd, test_name, run_type, 'compare')
        if (mpi_rank == 0):
            print '  CPRNC output directory will be:', cprnc_out_dir

        full_test_name = os.path.join(test_name, run_type)
        if (mpi_rank == 0):
            print '  Full test name:', full_test_name

        log_dir = os.path.join(cwd, full_test_name)
        log_filename = os.path.join(log_dir, 'check-' + test_name + '.log')
        if (mpi_rank == 0):
            print '  Log file will be:', log_filename

        full_test_names.append(full_test_name)
        results_filenames[full_test_name] = new_results_ls
        new_results_dirs[full_test_name] = new_results_dir
        old_results_dirs[full_test_name] = old_results_dir
        cprnc_out_dirs[full_test_name] = cprnc_out_dir
        log_filenames[full_test_name] = log_filename
        log_output[full_test_name] = []

    if (mpi_rank == 0):
        print

# Print out tests to be checked
if (mpi_rank == 0):
    print 'Valid tests for checking are:'
    for full_test_name in full_test_names:
        print '  ', full_test_name
    print

# Exit now, if only listing tests
if (options.list_tests):
    sys.exit(0)

# Reorganize the filenames to allow top-level looping over individual files
# and initialize comparison statistics
check_files = []
num_failures = {}
num_checks = {}
for full_test_name in full_test_names:
    num_failures[full_test_name] = 0
    num_checks[full_test_name] = len(results_filenames[full_test_name])
    for results_file in results_filenames[full_test_name]:
        check_file = {}
        check_file['full_test_name'] = full_test_name
        check_file['file_name'] = results_file
        check_files.append(check_file)

#==============================================================================
# Run comparison tests
#==============================================================================

# Base arguments to running cprnc
cprnc_args = ["/glade/p/work/kpaul/installs/cprnc/bin/cprnc", "-m", "", ""]

# String indicating IDENTITY (valid comparison)
ident_str = "files seem to be IDENTICAL"

# Before checking files, create CPRNC output directories for each test
for full_test_name in full_test_names:
    cprnc_out_dir = cprnc_out_dirs[full_test_name]
    if (mpi_rank == 0):
        if (not os.path.isdir(cprnc_out_dir)):
            os.makedirs(cprnc_out_dir)

sync()
if (mpi_rank == 0):
    print 'Checking files:'
    print

# Loop over all files (designated for this processor)
for check_file in check_files[mpi_rank::mpi_size]:
    file_name = check_file['file_name']
    full_test_name = check_file['full_test_name']

    new_results_dir = new_results_dirs[full_test_name]
    old_results_dir = old_results_dirs[full_test_name]
    cprnc_out_dir = cprnc_out_dirs[full_test_name]

    # Create the old and new file paths
    new_results_path = os.path.join(new_results_dir, file_name)
    old_results_path = os.path.join(old_results_dir, file_name)

    # CPRNC output file
    cprnc_out_filename = file_name + '.cprnc'
    cprnc_out_path = os.path.join(cprnc_out_dir, cprnc_out_filename)

    # If the output file already exists, then read it, otherwise write it
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
        num_failures[full_test_name] += 1
        log_result = '* BAD: '
    log_str = log_result + file_name + os.linesep
    log_output[full_test_name].append(log_str)

    print '   ', log_result + file_name

# Now wait for all files to be processed
if (not options.serial):
    MPI.COMM_WORLD.Barrier()

# Gather all of the statitics and log files from all processors
all_log_output = log_output
all_num_failures = num_failures
if (not options.serial):
    all_log_output = MPI.COMM_WORLD.gather(log_output)
    all_num_failures = MPI.COMM_WORLD.gather(num_failures)

# Open and write all log files
if (mpi_rank == 0):
    for full_test_name in full_test_names:
        log_file = open(log_filenames[full_test_name], 'w')
        for rank_log_output in all_log_output:
            log_strings = rank_log_output[full_test_name]
            for log_str in log_strings:
                log_file.write(log_str)
        log_file.close()

# Print diagnostic info to stdout
if (mpi_rank == 0):
    print
    print 'Number of File Comparison Failures:'
    print
    for full_test_name in full_test_names:
        num_test_failures = 0
        for rank_num_failures in all_num_failures:
            num_test_failures += rank_num_failures[full_test_name]
        print '  ', full_test_name + ':', num_test_failures,
        print 'failures (of', num_checks[full_test_name], 'total comparisons)'
