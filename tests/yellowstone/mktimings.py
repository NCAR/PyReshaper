#!/usr/bin/env python
#======================================================================
#
#  This script is designed to extract timing data from the scripttest
#  output.  It writes it to a JSON file, adding to the JSON file if the
#  particular log has not been recorded.  Overwrites if the log is already
#  present
#
#======================================================================

import os
import glob
import datetime
import optparse
import json

#==============================================================================
# Command-Line Interface Definition
#==============================================================================
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage=usage)
parser.add_option('-i', '--test_info', default=None,
                  help='Location of the testinfo.json file '
                       '[Default: None]')

# Parse the CLI options and assemble the Reshaper inputs
(options, arguments) = parser.parse_args()

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


#==============================================================================
# find_shortest_str - Helper Function
#==============================================================================
def find_shortest_str(strng, left, right=os.linesep, loc=0):
    sloc = strng.find(left, loc)
    if (sloc < 0):
        return '-1', 0
    sloc += len(left)
    eloc = strng.find(right, sloc)
    return strng[sloc:eloc].strip(), eloc

#==============================================================================
# Extract timing info and put in JSON file
#==============================================================================

# Open JSON file
json_name = 'timings.json'
json_data = {}
if (os.path.isfile(json_name)):
    json_file = open(json_name)
    json_data = dict(json.load(json_file))
    json_file.close()

# Current working directory
cwd = os.getcwd()

# Parse each test argument for test directories
test_names = glob.glob(os.path.join('*', 'ser'))
test_names.extend(glob.glob(os.path.join('*', 'par*')))

# Extract each possible test
for full_test_name in test_names:
    print
    print 'Extracting times from test dir:', full_test_name
    temp, run_type = os.path.split(full_test_name)
    temp, test_name = os.path.split(temp)
    print '  Test Name:', test_name
    print '  Run Type:', run_type

    common_name = test_info[test_name]['common_name']
    print '  Common Name:', common_name
    ncfmt = test_info[test_name]['netcdf_format']
    method_name = 'pyreshaper4c'
    if (ncfmt == 'netcdf'):
        method_name = 'pyreshaper'
    elif (ncfmt == 'netcdf4'):
        method_name = 'pyreshaper4'

    # Check for this test in the JSON data
    if (common_name not in json_data):
        json_data[common_name] = {'results': {}}

    # Check for method in results
    if (method_name not in json_data[common_name]['results']):
        json_data[common_name]['results'][method_name] = {}

    # Get the number of cores from the run type
    num_cores = 1
    if (run_type[0:3] == 'par'):
        num_cores = int(run_type[3:])

    # Look for log files
    test_dir = os.path.join(cwd, test_name, run_type)
    glob_name = ''.join(['reshaper-', test_name, '.*.log'])
    glob_path = os.path.join(test_dir, glob_name)
    glob_names = glob.glob(glob_path)

    # Continue if nothing to do here
    if (len(glob_names) == 0):
        print '  No log files found. Skipping.'
        continue

    # For each log file, extract the necessary info
    for log_name in glob_names:

        # Get the JOBID from the log filename
        job_id = log_name.rsplit('.', 2)[1]
        print '  Processing Job ID:', job_id

        # Open the log file and read the contents
        log_file = open(log_name)
        log_str = log_file.read()
        log_file.close()

        # Start the search through the file at character 0, and do it in
        # order so that you only need to search through it once...
        loc = 0

        # Look for the use of a metadata "once" file...
        once_file_str, once_loc = find_shortest_str(log_str,
           'Closed "once" ', loc=0)
        used_once_file = False
        if (once_loc > 0):
            used_once_file = True

        # Find the internal timing data from the run output:
        openi_str, loc = find_shortest_str(log_str,
           'Open Input Files: ', loc=loc)
        tot_timing_str, loc = find_shortest_str(log_str,
            'Complete Conversion Process: ', loc=loc)
        openo_str, loc = find_shortest_str(log_str,
            'Open Output Files: ', loc=loc)
        metaTI_str, loc = find_shortest_str(log_str,
            'Write Time-Invariant Metadata: ', loc=loc)
        metaTV_str, loc = find_shortest_str(log_str,
            'Write Time-Variant Metadata: ', loc=loc)
        TS_str, loc = find_shortest_str(log_str,
            'Write Time-Series Variables: ', loc=loc)

        # Get read/write size amounts
        actual_str, nloc = find_shortest_str(log_str,
            'Actual Data: ', loc=loc)
        requested_str, nloc = find_shortest_str(log_str,
            'Requested Data: ', loc=loc)

        # Find the job hosts string
        host_str, loc = find_shortest_str(log_str,
            'Job was executed on host(s) ', loc=loc)
        if (loc == 0):
            print '    Could not find host string.  Skipping.'
            continue
        queue_str, sloc = find_shortest_str(host_str,
            'in queue <', right='>', loc=0)
        if (sloc == 0):
            print '    Could not find queue string.  Skipping.'
            continue
        print '    Queue:', queue_str
        user_str, sloc = find_shortest_str(host_str,
            'as user <', right='>', loc=sloc)
        if (sloc == 0):
            print '    Could not find user string.  Skipping.'
            continue
        print '    User:', user_str
        sys_str, sloc = find_shortest_str(host_str,
            'in cluster <', right='>', loc=sloc)
        if (sloc == 0):
            print '    Could not find system string.  Skipping.'
            continue
        print '    System:', sys_str

        # Get the starting timestamp
        start_str, loc = find_shortest_str(log_str,
            'Started at ', loc=loc)
        if (loc == 0):
            print '    Could not find start timestamp.  Skipping.'
            continue
        start = datetime.datetime.strptime(start_str, '%c')
        print '    Starting Timestamp:', start

        # Get the ending timestamp
        end_str, loc = find_shortest_str(log_str,
            'Results reported at ', loc=loc)
        if (loc == 0):
            print '    Could not find end timestamp.  Skipping.'
            continue
        end = datetime.datetime.strptime(end_str, '%c')
        print '    Ending Timestamp:', end

        # Compute the job number from the end timestamp
        job_num = end.strftime('%y%m%d-%H%M%S')

        # Compute the elapsed time
        elapsed = (end - start).total_seconds()
        print '    Elapsed time (REAL):    ', elapsed, 'sec'

        # Display the internal timing info for comparison
        print '    Elapsed time (INTERNAL):', float(tot_timing_str), 'sec'

        # Find the BSUB submission properties
        nprocs_str, loc = find_shortest_str(log_str,
            '#BSUB -n ', loc=loc)
        if (loc == 0):
            print '    Could not find num procs string.  Skipping.'
            continue
        print '    BSUB Cores:', nprocs_str
        tiling_str, loc = find_shortest_str(log_str,
            '#BSUB -R "span[ptile=', right=']', loc=loc)
        if (loc == 0):
            print '    Could not find num procs string.  Skipping.'
            continue
        print '    BSUB Tiling:', tiling_str
        num_nodes = int(nprocs_str) / int(tiling_str)

        # Write the JSON data
        if (job_num not in json_data[common_name]['results'][method_name]):
            json_data[common_name]['results'][method_name][job_num] = {}

        json_data[common_name]['results'][method_name][job_num]['sys'] = \
            sys_str
        json_data[common_name]['results'][method_name][job_num]['cores'] = \
            num_cores
        json_data[common_name]['results'][method_name][job_num]['nodes'] = \
            num_nodes
        json_data[common_name]['results'][method_name][job_num]['real'] = \
            elapsed
        json_data[common_name]['results'][method_name][job_num]['metadata'] = \
            True
        json_data[common_name]['results'][method_name][job_num]['once'] = \
            used_once_file
        json_data[common_name]['results'][method_name][job_num]['actual'] = \
            float(actual_str)
        json_data[common_name]['results'][method_name][job_num]['request'] = \
            float(requested_str)
        json_data[common_name]['results'][method_name][job_num]['openi'] = \
            float(openi_str)
        json_data[common_name]['results'][method_name][job_num]['openo'] = \
            float(openo_str)
        json_data[common_name]['results'][method_name][job_num]['metaTI'] = \
            float(metaTI_str)
        json_data[common_name]['results'][method_name][job_num]['metaTV'] = \
            float(metaTV_str)
        json_data[common_name]['results'][method_name][job_num]['TS'] = \
            float(TS_str)

# Write the JSON data file
json_file = open(json_name, 'w')
json_file.write(json.dumps(json_data, sort_keys=True,
    indent=4, separators=(',', ': ')))
json_file.close()
