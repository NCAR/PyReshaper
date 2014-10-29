#!/usr/bin/env python
#==============================================================================
#
#  scripttest.py
#
#  This script is written to make running the "Bake-Off" time-slice to
#  time-series tests on Yellowstone easy and efficient.  This allows code
#  re-use for the various tests.
#
#  This script requires the testinfo.json file be present in the same directory
#  as the runtest script itself.
#
#  This script runs a single test at a time.  It does not use the multi-spec
#  reshaper.  It runs the single-stream script: slice2series
#
#==============================================================================
import optparse
import os
import sys
import stat
import shutil
import glob

from subprocess import Popen, PIPE, STDOUT
import json

#==============================================================================
# Command-Line Interface Definition
#==============================================================================
usage = 'usage: %prog [options] test_name1 test_name2 ...\n\n' \
      + '       This program is designed to run yellowstone-specific\n' \
      + '       tests of the PyReshaper.  Each named test (or all tests if\n' \
      + '       the -a or --all option is used) will be given a run\n' \
      + '       directory in the current working directory with the same\n' \
      + '       name as the test itself.  The run script will be placed\n' \
      + '       in this run directory, as will be placed the run output/\n' \
      + '       error file.  All output data files will be placed in the\n' \
      + '       output subdirectory.'
parser = optparse.OptionParser(usage=usage)
parser.add_option('-a', '--all', default=False,
                  action='store_true', dest='all',
                  help='True or False, indicating whether to run all tests '
                       '[Default: False]')
parser.add_option('-c', '--code', default='STDD0002',
                  help='The name of the project code for charging in '
                       'parallel runs (ignored if running in serial) '
                       '[Default: STDD0002]')
parser.add_option('-f', '--force', default=False,
                  action='store_true', dest='force',
                  help='True or False, indicating whether to force deleting '
                       'any existing test or run directories, if found '
                       '[Default: False]')
parser.add_option('-i', '--test_info', default=None,
                  help='Location of the testinfo.json file '
                       '[Default: None]')
parser.add_option('-l', '--list', default=False,
                  action='store_true', dest='list_tests',
                  help='True or False, indicating whether to list all tests, '
                       'instead of running tests. [Default: False]')
parser.add_option('-o', '--only', default=0, type='int',
                  help='Runs each test limited to the indicated number of '
                       'output files per processor.  A value of 0 means '
                       'no limit (i.e., all output files) [Default: 0]')
parser.add_option('--once', default=False,
                  action='store_true', dest='once_file',
                  help='True or False, indicating whether to write metadata '
                       'to a once file. [Default: False]')
parser.add_option('-p', '--processes', default=0, type='int',
                  help='The integer number of processes to request in parallel'
                       ' runs (0 means run in serial) [Default: 0]')
parser.add_option('-q', '--queue', default='economy',
                  help='The name of the queue to request in parallel runs '
                       '(ignored if running in serial) '
                       '[Default: economy]')
parser.add_option('-t', '--tiling', default=16, type='int',
                  help='The integer number of processes per node to request '
                       'in parallel runs (ignored if running in serial) '
                       '[Default: 16]')
parser.add_option('-w', '--wtime', default=240, type='int',
                  help='The number of minutes to request for the wall clock '
                       'in parallel runs (ignored if running in serial) '
                       '[Default: 240]')

# Parse the CLI options and assemble the Reshaper inputs
(options, arguments) = parser.parse_args()

#==============================================================================
# Initial pre-processing before running test
#==============================================================================

# For convenience, pull out the useful option arguments
num_procs = options.processes
num_procs_per_node = options.tiling
queue_name = options.queue
project_code = options.code
output_limit = options.only

wtime_hours = options.wtime / 60
if (wtime_hours > 99):
    wtime_hours = 99
    print 'Requested number of hours too large.  Limiting to', wtime_hours, '.'
wtime_minutes = options.wtime % 60
wtime_str = '%02d:%02d' % (wtime_hours, wtime_minutes)

# See if there is a user-defined testinfo file, otherwise look for default
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

# List tests only, if option is set
if (options.list_tests):
    print 'Tests found in the Test Info file are:'
    print
    for test_name in test_info:
        print '   ' + str(test_name)
    sys.exit(0)

# Determine which tests to run
test_names = arguments
if (options.all):
    test_names = test_info.keys()
else:
    for test_name in test_names:
        if (test_name not in test_info):
            print str(test_name) + ' not in testinfo file.  Ignoring.'
print 'Tests to be run:',
for test_name in test_names:
    print str(test_name),
print ' '
print


#==============================================================================
# gen_lsf_preamble - write LSF preamble as a string
#==============================================================================
def gen_lsf_preamble(test_name):
    run_script_str = [
        '#BSUB -n ' + str(num_procs),
        '#BSUB -R "span[ptile=' + str(num_procs_per_node) + ']"',
        '#BSUB -q ' + queue_name,
        '#BSUB -a poe',
        '#BSUB -x',
        '#BSUB -o reshaper-' + test_name + '.%J.log',
        '#BSUB -J reshaper-' + test_name,
        '#BSUB -P ' + project_code,
        '#BSUB -W ' + wtime_str,
        '',
        'export MP_TIMEOUT=14400',
        'export MP_PULSE=1800',
        'export MP_DEBUG_NOTIMEOUT=yes']
    return run_script_str


#==============================================================================
# gen_run_script - write bash run script as string
#==============================================================================
def gen_run_script(test_name, output_dir):

    # Start creating the run scripts for each test
    run_script_str = ['#!/bin/bash']

    # If necessary, add the parallel preamble
    if (num_procs > 0):
        run_script_str.extend(gen_lsf_preamble(test_name))

    # Now create the rest of the run script
    run_script_str.extend(['',
        '# NOTE: Your PATH and PYTHONPATH must be properly set',
        '#       before this script will run without error',
        '',
        '# Necessary modules to load',
        'module load python',
        'module load all-python-libs',
        ''])

    # Generate the command line arguments
    run_cmd_str = ''
    if (num_procs > 0):
        run_cmd_str += 'mpirun.lsf '
    run_cmd_str += 'slice2series '
    if (num_procs == 0):
        run_cmd_str += '--serial '
    if (options.once_file):
        run_cmd_str += '--once '
    run_cmd_str += '-v3 '
    if (output_limit > 0):
        run_cmd_str += '-l' + str(output_limit) + ' '

    format_str = str(test_info[test_name]['netcdf_format'])
    run_cmd_str += '-f ' + format_str + ' '

    prefix_str = os.path.join(output_dir,
                              str(test_info[test_name]['output_prefix']))
    run_cmd_str += '-p ' + prefix_str + ' '

    suffix_str = str(test_info[test_name]['output_suffix'])
    run_cmd_str += '-s ' + suffix_str + ' '

    for var_name in test_info[test_name]['metadata']:
        run_cmd_str += '-m ' + str(var_name) + ' '

    input_files_list = []
    for input_glob in test_info[test_name]['input_globs']:
        full_input_glob = os.path.join(str(test_info[test_name]['input_dir']),
                                       input_glob)
        input_files_list.extend(glob.glob(full_input_glob))
    run_cmd_str += ' '.join(input_files_list)

    run_script_str.append(run_cmd_str)

    return run_script_str


#==============================================================================
# For each test, create the run directory, script, output dir, and launch
#==============================================================================
cwd = os.getcwd()
for test_name in test_names:
    print 'Currently preparing test:', test_name

    # Create the test directory
    test_dir = os.path.abspath(test_name)
    if (num_procs > 0):
        test_dir = os.path.join(test_dir, 'par' + str(num_procs))
    else:
        test_dir = os.path.join(test_dir, 'ser')

    # Delete old directory, if forced
    if (os.path.isdir(test_dir)):
        if (options.force):
            shutil.rmtree(test_dir, ignore_errors=True)
            print '  Directory (' + test_dir + ') removed.'
        else:
            err_msg = '  Test directory (' + test_dir + ') already exists.'
            raise RuntimeError(err_msg)

    # Make new test directory
    os.makedirs(test_dir)
    print '  Test directory (' + test_dir + ') created.'

    # Create the output subdirectory
    output_dir = os.path.join(test_dir, 'output')
    os.makedirs(output_dir)
    print '  Output data directory (' + output_dir + ') created.'

    # Generate the run script
    run_script_str = ''
    run_script_filename = 'run-' + test_name + '.sh'
    run_script_abspath = os.path.join(test_dir, run_script_filename)
    run_script_str = gen_run_script(test_name, output_dir)
    run_script_file = open(run_script_abspath, 'w')
    run_script_file.write(os.linesep.join(run_script_str))
    run_script_file.close()
    print '  Run script written to', run_script_abspath

    # Now launch the test
    os.chdir(test_dir)
    print '  Launching test job:',
    if (num_procs == 0):
        # Make the script executable
        os.chmod(run_script_abspath,
                 stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)

        # Launch the serial job as a subprocess
        job = Popen([run_script_abspath], stdout=PIPE, stderr=STDOUT,
                    env=os.environ.copy())
        print 'PID:', str(job.pid)
        sys.stdout.flush()

        # Wait for job to finish and grab job output
        job_output = job.communicate()[0]

        # Write output to log file
        log_file = open('reshaper-' + test_name + '.log', 'w')
        log_file.write(job_output)
        log_file.close()

    else:
        # Open up the run script for input to LSF's bsub
        run_script_file = open(run_script_abspath, 'r')

        # Launch the parallel job with LSF bsub
        job = Popen(['bsub'], stdout=PIPE, stderr=STDOUT,
                    stdin=run_script_file, env=os.environ.copy())

        # Grab the bsub output
        job_output = job.communicate()[0]

        # Close the script file and print submission info
        run_script_file.close()
        print job_output
        sys.stdout.flush()

    os.chdir(cwd)
