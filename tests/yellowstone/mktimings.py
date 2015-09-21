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
import argparse
import json

# Package Modules
from utilities import testtools as tt

#==============================================================================
# Command-Line Interface Definition
#==============================================================================
_DESC_ = """This program is designed to gather statistics for tests and
            test input defined in the testing database file."""

_PARSER_ = argparse.ArgumentParser(description=_DESC_)
_PARSER_.add_argument('-i', '--infofile', default='testinfo.json', type=str,
                      help='Location of the testinfo.json database file '
                           '[Default: testinfo.json]')
_PARSER_.add_argument('-s', '--statsfile', default='teststats.json', type=str,
                      help='Location of the teststats.json database file '
                           '[Default: teststats.json]')
_PARSER_.add_argument('-t', '--timefile', default='timings.json', type=str,
                      help='Location of the timings.json database file '
                           '[Default: timings.json]')


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
# Command-line Operation
#==============================================================================
if __name__ == '__main__':
    args = _PARSER_.parse_args()

    # Create/read the testing info and stats files
    testdb = tt.TestDB(dbname=args.infofile, stname=args.statsfile).get_database()

    # Get the timings.json data file
    timefn = os.path.abspath(args.timefile)

    # Open JSON file
    timedb = {}
    if (os.path.isfile(timefn)):
        json_file = open(timefn)
        timedb = dict(json.load(json_file))
        json_file.close()

    # Current working directory
    cwd = os.getcwd()

    # Extract each possible test
    for rundir in glob.iglob(os.path.join('results.d', '*', '[ser,par]*', '*')):
        print
        print 'Extracting times from test dir:', rundir
        root, ncfmt = os.path.split(rundir)
        root, run_type = os.path.split(root)
        root, test_name = os.path.split(root)
        print '  Test Name:', test_name
        print '  Run Type:', run_type
        print '  NetCDF Format:', ncfmt

        # Skip if not an individual test
        if test_name not in testdb:
            print '  Test name not found in database. Skipping.'
            continue

        # Prepare the timing database information
        common_name = testdb[test_name]['common_name']
        print '  Common Name:', common_name
        method_name = 'pyreshaper4c'
        if ncfmt == 'netcdf':
            method_name = 'pyreshaper'
        elif ncfmt == 'netcdf4':
            method_name = 'pyreshaper4'

        # Check for this test in the timing database
        if common_name not in timedb:
            timedb[common_name] = {'results': {}}

        # Check for method in results
        if method_name not in timedb[common_name]['results']:
            timedb[common_name]['results'][method_name] = {}

        # Get the number of cores and nodes from the run type
        num_cores = 1
        num_nodes = 0
        if run_type[0:3] == 'par':
            num_nodes, nppn = map(int, run_type[3:].split('x'))
            num_cores = num_nodes * nppn

        # Look for log files
        glob_path = os.path.join(rundir, '{0!s}*.log'.format(test_name))
        glob_names = glob.glob(glob_path)

        # Continue if nothing to do here
        if (len(glob_names) == 0):
            print '  No log files found. Skipping.'
            continue

        # For each log file, extract the necessary info
        for log_name in glob_names:

            # Get the JOBID from the log filename
            lognm = os.path.split(log_name)[1]
            print '  Processing log:', lognm

            # Open the log file and read the contents
            log_file = open(log_name)
            log_str = log_file.read()
            log_file.close()

            # Start the search through the file at character 0, and do it in
            # order so that you only need to search through it once...
            loc = 0

            # Look for the use of a metadata "once" file...
            once_file_str, once_loc = find_shortest_str(
                log_str, 'Closed "once" ', loc=0)
            used_once_file = False
            if (once_loc > 0):
                used_once_file = True

            # Find the internal timing data from the run output:
            openi_str, loc = find_shortest_str(
                log_str, 'Open Input Files: ', loc=loc)
            tot_timing_str, loc = find_shortest_str(
                log_str, 'Complete Conversion Process: ', loc=loc)
            openo_str, loc = find_shortest_str(
                log_str, 'Open Output Files: ', loc=loc)
            metaTI_str, loc = find_shortest_str(
                log_str, 'Write Time-Invariant Metadata: ', loc=loc)
            metaTV_str, loc = find_shortest_str(
                log_str, 'Write Time-Variant Metadata: ', loc=loc)
            TS_str, loc = find_shortest_str(
                log_str, 'Write Time-Series Variables: ', loc=loc)

            # Get read/write size amounts
            actual_str, nloc = find_shortest_str(
                log_str, 'Actual Data: ', loc=loc)
            requested_str, nloc = find_shortest_str(
                log_str, 'Requested Data: ', loc=loc)

            # Compute the job number from the log file timestamp
            timestamp = os.path.getmtime(log_name)
            dt = datetime.datetime.fromtimestamp(timestamp)
            job_num = dt.strftime('%y%m%d-%H%M%S')

            # Compute the elapsed time
            elapsed = float(tot_timing_str)

            # Display the internal timing info for comparison
            print '    Elapsed time (INTERNAL):', elapsed, 'sec'

            # Write the JSON data
            if (job_num not in timedb[common_name]['results'][method_name]):
                timedb[common_name]['results'][method_name][job_num] = {}

            timedb[common_name]['results'][
                method_name][job_num]['sys'] = "yellowstone"
            timedb[common_name]['results'][
                method_name][job_num]['cores'] = num_cores
            timedb[common_name]['results'][
                method_name][job_num]['nodes'] = num_nodes
            timedb[common_name]['results'][
                method_name][job_num]['real'] = elapsed
            timedb[common_name]['results'][
                method_name][job_num]['metadata'] = True
            timedb[common_name]['results'][
                method_name][job_num]['once'] = used_once_file
            timedb[common_name]['results'][
                method_name][job_num]['actual'] = float(actual_str)
            timedb[common_name]['results'][
                method_name][job_num]['request'] = float(requested_str)
            timedb[common_name]['results'][
                method_name][job_num]['openi'] = float(openi_str)
            timedb[common_name]['results'][
                method_name][job_num]['openo'] = float(openo_str)
            timedb[common_name]['results'][
                method_name][job_num]['metaTI'] = float(metaTI_str)
            timedb[common_name]['results'][
                method_name][job_num]['metaTV'] = float(metaTV_str)
            timedb[common_name]['results'][
                method_name][job_num]['TS'] = float(TS_str)

    # Write the JSON data file
    json_file = open(timefn, 'w')
    json_file.write(json.dumps(timedb, sort_keys=True,
                               indent=4, separators=(',', ': ')))
    json_file.close()
