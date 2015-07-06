#!/usr/bin/env python
#==============================================================================
#
#  runtestsm.py
#
#  This script is written to make running the "Bake-Off" time-slice to
#  time-series tests on Yellowstone easy and efficient.  This allows code
#  re-use for the various tests.
#
#  This script requires the testinfo.json file be present in the same directory
#  as the runtest script itself.
#
#  This script runs MULTIPLE tests at one time through a single instance of
#  the multi-specifier reshaper class.
#
#==============================================================================
import optparse
import os
import sys
import shutil
import glob

from subprocess import Popen, PIPE, STDOUT
import json

from runtestsi import *


#==============================================================================
# Generate the test run command
#==============================================================================
def create_run_command(options, python_script_name):

    # Generate the command line run_args
    run_cmd = []
    if (options.nodes > 0):
        run_cmd.append('mpirun.lsf')
    run_cmd.extend(['python', python_script_name])

    return ' '.join(run_cmd)


#==============================================================================
# Generate the python script
#==============================================================================
def create_python_script(options, test_name, output_dir, testing_database):

    # Start the python script strings
    python_script_str = [
        '#!/usr/bin/env python',
        'import glob',
        '#',
        '# PYTHONPATH must be set for this to work',
        '#',
        'from pyreshaper import specification',
        'from pyreshaper import reshaper',
        '',
        'spec_list = {}',
        '']

    # Define the Specifiers for each test
    input_dir = str(testing_database[test_name]['input_dir'])
    input_file_list = []
    for input_glob in testing_database[test_name]['input_globs']:
        input_glob_str = str(input_glob)
        full_input_glob = str(os.path.join(input_dir, input_glob))
        input_file_list.extend(glob.glob(full_input_glob))
    netcdf_format = str(testing_database[test_name]['netcdf_format'])
    output_prefix = str(
        os.path.join(output_dir,
                     str(testing_database[test_name]['output_prefix'])))
    output_suffix = str(testing_database[test_name]['output_suffix'])
    time_variant_metadata = map(
        str, testing_database[test_name]['metadata'])

    python_script_str.extend([
        'spec = specification.create_specifier()',
        'spec.input_file_list = ' + repr(input_file_list),
        'spec.netcdf_format = ' + repr(netcdf_format),
        'spec.output_file_prefix = ' + repr(output_prefix),
        'spec.output_file_suffix = ' + repr(output_suffix),
        'spec.time_variant_metadata = ' + repr(time_variant_metadata),
        'spec_list[' + repr(str(test_name)) + '] = spec',
        ''])

    # Now create the multi-spec reshaper and run
    ser = str(num_procs == 0)
    reshaper_params = 'spec_list, serial=' + ser + ', verbosity=2'
    python_script_str.extend([
        'rshpr = reshaper.create_reshaper(' + reshaper_params + ')',
        'rshpr.convert(output_limit=' + str(output_limit) + ')',
        'rshpr.print_diagnostics()'
        ''])

    # Open the python script file and write script
    python_script_name = 'run-' + test_name + '.py'
    python_script_file = open(python_script_name, 'w')
    python_script_file.write(os.linesep.join(python_script_str))
    python_script_file.close()

    return python_script_name


#==============================================================================
# Command-Line Operation
#==============================================================================
if __name__ == '__main__':
    options, arguments = parse_cli()
    testing_database = get_testing_database(options)
    tests_to_run = get_tests_to_run(options, arguments, testing_database)
    for test_name in tests_to_run:
        print 'Currently preparing test:', test_name

        test_dir = get_test_dirname(options, test_name)
        create_dir('test', test_dir)

        output_dir = os.path.join(test_dir, 'output')
        create_dir('output', output_dir)

    python_script = create_python_script(
        options, test_name, output_dir, testing_database)

    run_command = create_run_command(options, python_script)

    run_script = create_run_script(
        options, test_name, test_dir, run_command)

    #run_test(test_name, test_dir, run_script)
