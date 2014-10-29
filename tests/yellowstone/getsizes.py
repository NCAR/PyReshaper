#!/usr/bin/env python
#==============================================================================
#
# This script is designed to look through all of the datasets and determine
# there total data sizes.  In particular, it determines the
# total size, size of all time-series variables, and metadata sizes.
#
#==============================================================================

import os
import Nio
import json
import glob
import optparse
import numpy as np
np.set_printoptions(precision=2)

#==============================================================================
# Command-Line Interface Definition
#==============================================================================
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage=usage)
parser.add_option('-i', '--testinfo', default=None,
                  help='Location of the testinfo.json file '
                       '[Default: None]')

# Parse the CLI options and assemble the Reshaper inputs
(options, arguments) = parser.parse_args()

# Read the data file
data_file_name = 'testinfo.json'
if (options.testinfo != None and os.path.isfile(options.testinfo)):
    data_file_name = options.testinfo
data_file = open(data_file_name)
jsondata = dict(json.load(data_file))
data_file.close()

#==============================================================================
# Helper function to convert typecodes into dtype itemsizes
def get_itemsize(tc):
    dt = None
    if (tc == 'd'):
        dt = np.float64
    elif (tc == 'f'):
        dt = np.float32
    elif (tc == 'l'):
        dt = np.long
    elif (tc == 'i'):
        dt = np.int32
    elif (tc == 'h'):
        dt = np.int16
    elif (tc == 'b'):
        dt = np.int8
    elif (tc == 'S1'):
        dt = np.character
    else:
        dt = np.float
    return np.dtype(dt).itemsize

#==============================================================================
# Helper function to return the size of an array given its shape
def get_size(shp):
    if (len(shp) > 0):
        return reduce(lambda x, y: x * y, shp)
    else:
        return 1

#==============================================================================
# Helper function to represent an integer number of bytes as a string w/ units
def nbyte_str(n, exp=0):
    if (n > 1024.):
        return nbyte_str(n / 1024., exp=exp + 1)
    else:
        units = ''
        if (exp == 0):
            units = 'Bytes'
        elif (exp == 1):
            units = 'KB'
        elif (exp == 2):
            units = 'MB'
        elif (exp == 3):
            units = 'GB'
        elif (exp == 4):
            units = 'TB'
        elif (exp == 5):
            units = 'PB'
        else:
            n *= 1024.**(exp - 5)
            units = 'PB'
        return '%.3f %s' % (n, units)

#==============================================================================
# Loop over all tests
#==============================================================================
for test_name in jsondata:
    print test_name + ':'

    # Generate the list of input files
    input_dir = jsondata[test_name]['input_dir']
    input_files = []
    for glob_str in jsondata[test_name]['input_globs']:
        full_glob_str = os.path.join(input_dir, glob_str)
        glob_files = glob.glob(full_glob_str)
        input_files.extend(glob_files)

    # Sort by name
    input_files.sort()

    # Open the first file
    infile0 = Nio.open_file(input_files[0], 'r')

    # Get the name of the unlimited dimension (e.g., time)
    udim = None
    for dim in infile0.dimensions:
        if infile0.unlimited(dim):
            udim = dim
            continue

    # Get the data dimensions
    metadata_names = set(infile0.dimensions.keys())

    # Add the extra metadata variable names
    metadata_names.update(set(jsondata[test_name]['metadata']))

    # Determine sizes of variables
    timeseries_vars = {}
    ti_metadata_vars = {}
    tv_metadata_vars = {}
    for var_name in infile0.variables.keys():
        var_obj = infile0.variables[var_name]
        var_shape = var_obj.shape
        var_size = get_size(var_shape)
        var_bytesize = var_size * get_itemsize(var_obj.typecode())
        if (udim in var_obj.dimensions):
            if var_name in metadata_names:
                tv_metadata_vars[var_name] = var_bytesize
            else:
                timeseries_vars[var_name] = var_bytesize
        else:
            ti_metadata_vars[var_name] = var_bytesize

    tser_num = len(timeseries_vars)
    timd_num = len(ti_metadata_vars)
    tvmd_num = len(tv_metadata_vars)
    tot_num = tser_num + timd_num + tvmd_num

    # Close the file
    infile0.close()

    # Some output (numbers)
    print '  No. of TI Metadata Variables: ', timd_num
    print '  No. of TV Metadata Variables: ', tvmd_num
    print '  No. of Time-Series Variables: ', tser_num
    print '  No. of Variables (TOTAL):     ', tot_num

    # Loop over all input files and compute time-variant variable sizes
    for file_name in input_files[1:]:
        # Open the first file
        infile = Nio.open_file(file_name, 'r')

        # Compute variable sizes
        for var_name in infile.variables.keys():
            var_obj = infile.variables[var_name]
            if (udim in var_obj.dimensions):
                var_shape = var_obj.shape
                var_size = get_size(var_shape)
                var_bytesize = var_size * get_itemsize(var_obj.typecode())
                if var_name in metadata_names:
                    tv_metadata_vars[var_name] += var_bytesize
                else:
                    timeseries_vars[var_name] += var_bytesize

        # Close the file
        infile.close()

    # Compute totals
    tser_size = np.sum(timeseries_vars.values())
    timd_size = np.sum(ti_metadata_vars.values())
    tvmd_size = np.sum(tv_metadata_vars.values())
    tot_isize = tser_size + timd_size + tvmd_size

    print '  Size of TI Metadata Variables:', nbyte_str(timd_size)
    print '  Size of TV Metadata Variables:', nbyte_str(tvmd_size)
    print '  Size of Time-Series Variables:', nbyte_str(tser_size)
    print '  TOTAL Size of Variables:      ', nbyte_str(tot_isize)
    print
