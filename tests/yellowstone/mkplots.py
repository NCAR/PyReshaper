#!/usr/bin/env python
#==============================================================================
#
# This script is based on the script used to create the plots used in the
# poster and presentation at the 2014 CESM Workshop in Breckenridge, CO.  It
# is specialized to the data generated for the Workshop presentations, but
# it has been expanded and generalized in some respects.
#
# This script is meant to be run from the command-line with no options. It will
# search all of the tests found in the 'timings.json' data file.  It will pull
# out the timing and throughput data for each test in the data file, but it
# will only plot the timing and throughput data for tests and runs that
# satisfy desired conditions.  Namely, these conditions are that the plotted
# results are for test runs that:
#
#   1. Must include all metadata
#   2. Must use only a set number of cores
#   3. Method must be used in all other datasets
#   4. Pick the most recent test run, if there are multiple test runs
#
#==============================================================================

import os
import json
import numpy as np
import optparse
from matplotlib import pyplot as plt
plt.switch_backend('agg')

#==============================================================================
# Command-Line Interface Definition
#==============================================================================
usage = 'usage: %prog [options] [method1] [method2]'
parser = optparse.OptionParser(usage=usage)
parser.add_option('-f', '--data_file', default=None,
                  help='Location of the timings.json file '
                       '[Default: None]')
parser.add_option('-n', '--num_cores', default=16, type='int',
                  help='Default number of cores to assume in parallel plots '
                       '[Default: 16]')

# Parse the CLI options and assemble the Reshaper inputs
(options, arguments) = parser.parse_args()

# Read the data file
data_file_name = 'timings.json'
if (options.data_file != None and os.path.isfile(options.data_file)):
    data_file_name = options.data_file
data_file = open(data_file_name)
jsondata = dict(json.load(data_file))
data_file.close()

#==============================================================================
# Some common data for plotting
#==============================================================================
dataset_labels = {'POP-1.0': 'POP (1 deg)', 'POP-0.1': 'POP (0.1 deg)',
                  'CLM-1.0': 'CLM (1 deg)', 'CLM-0.25': 'CLM (1/4 deg)',
                  'CICE-1.0': 'CICE (1 deg)', 'CICE-0.1': 'CICE (0.1 deg)',
                  'CAMSE-1.0': 'CAM-SE (1 deg)',
                  'CAMSE-0.25': 'CAM-SE (1/4 deg)',
                  'CAMFV-1.0': 'CAM-FV (1 deg)'}
dataset_order = ['CAMFV-1.0', 'CAMSE-1.0', 'CICE-1.0', 'CLM-1.0', 'POP-1.0',
                 'CAMSE-0.25', 'CICE-0.1', 'CLM-0.25', 'POP-0.1']

method_labels = {'ncl': 'NCL 6.1.2',
                 'nco': 'NCO',
                 'ncr': 'ncReshaper',
                 'pynio': 'PyNIO (NetCDF3)',
                 'pynio4_0': 'PyNIO (NetCDF4)',
                 'pynio4_1': 'PyNIO (NetCDF4-CL1)',
                 'pyniompi': 'PyNIO+mpi4py (NetCDF3)',
                 'pyniompi4_0': 'PyNIO+mpi4py (NetCDF4)',
                 'pyniompi4_1': 'PyNIO+mpi4py (NetCDF4-CL1)',
                 'pagoda': 'Pagoda',
                 'cdo': 'CDO',
                 'pyreshaper': 'PyReshaper (NetCDF3)',
                 'pyreshaper4': 'PyReshaper (NetCDF4)',
                 'pyreshaper4c': 'PyReshaper (NetCDF4-CL1)'}
method_colors = {'ncl': 'magenta',
                 'nco': 'red',
                 'ncr': 'orange',
                 'pynio': 'purple',
                 'pynio4_0': 'green',
                 'pynio4_1': 'blue',
                 'pyniompi': 'purple',
                 'pyniompi4_0': 'green',
                 'pyniompi4_1': 'blue',
                 'pagoda': 'yellow',
                 'cdo': 'cyan',
                 'pyreshaper': 'purple',
                 'pyreshaper4': 'green',
                 'pyreshaper4c': 'blue'}
method_order = ['ncl', 'nco', 'pagoda', 'cdo',
                'pynio', 'pynio4_0', 'pynio4_1',
                'ncr',
                'pyniompi', 'pyniompi4_0', 'pyniompi4_1',
                'pyreshaper', 'pyreshaper4', 'pyreshaper4c']

#==============================================================================
# METHODS TO PLOT (From CL Arguments)
#==============================================================================
methods_to_plot = None
if (len(arguments) > 0):
    methods_to_plot = []
    for arg in arguments:
        if arg in method_labels:
            methods_to_plot.append(arg)


#==============================================================================
# Find all data according to the criteria:
# 1. Must include all metadata
# 2. Must use only 'ncores' cores
# 3. Method must be used in all datasets
# 4. Pick the most recent run, if there are multiple
# 5. Method must be in the methods_to_plot list
#==============================================================================
def find_subset(ncores=1):
    greater_subset = {}
    for dataset in jsondata:
        greater_subset[dataset] = {}
        results = jsondata[dataset]['results']
        for method in results:
            if (methods_to_plot != None and method not in methods_to_plot):
                continue
            runs = results[method]
            good_runs = []
            for run, rundata in runs.items():
                if ('cores' in rundata and 'metadata' in rundata):
                    if (rundata['metadata']):
                        if (rundata['cores'] == ncores):
                            good_runs.append(run)
                        elif (ncores > 1 and method == 'nco'):
                            good_runs.append(run)
            if (len(good_runs) > 0):
                most_recent_run = max(good_runs)
                greater_subset[dataset][method] = runs[most_recent_run]
                greater_subset[dataset][method]['isize'] = \
                    jsondata[dataset]['isize']
    dataset0 = greater_subset.keys()[0]
    common_methods = set(greater_subset[dataset0])
    for dataset in greater_subset:
        methods = greater_subset[dataset].keys()
        common_methods.intersection_update(methods)
    subset = {}
    for dataset in greater_subset:
        subset[dataset] = {}
        for method in greater_subset[dataset]:
            if method in common_methods:
                subset[dataset][method] = greater_subset[dataset][method]
    return subset, common_methods


#==============================================================================
# Extract the duration data in minutes from a data subset and method
#==============================================================================
def get_duration_data(data, method):
    yvalues = []
    for dataset in dataset_order:
        if dataset in data:
            yvalues.append(data[dataset][method]['real'] / 60.0)
    return np.array(yvalues)


#==============================================================================
# Extract the throughput data in MB/s from a data subset and method
#==============================================================================
def get_throughput_data(data, method):
    yvalues = []
    for dataset in dataset_order:
        if dataset in data:
            isize = data[dataset][method]['isize']
            time = data[dataset][method]['real']
            yvalues.append(isize / time)
    return np.array(yvalues)


#==============================================================================
# Get a function for retrieving the data with a name of type
#=============================================================================
def get_yvalues(name):
    if (name == 'throughput'):
        return get_throughput_data
    elif (name == 'duration'):
        return get_duration_data


#==============================================================================
# Get the ylabels from the name
#=============================================================================
def get_ylabel(name):
    if (name == 'throughput'):
        return 'Throughput [MB/sec]'
    elif (name == 'duration'):
        return 'Duration [min]'


#==============================================================================
# Make a plot from a data subset
#==============================================================================
def make_bar_plot(ncores, name, title, filename):

    # Get the data
    data, methods = find_subset(ncores)

    # Plot-sizing data
    width = 0.75 / len(methods)
    offset = -0.375
    xbase = np.arange(len(data)) + 1
    plt.figure(figsize=(14, 8))
    plt.subplots_adjust(left=0.11, right=0.975, top=0.915, bottom=0.275)

    # Plot every method and dataset
    for method in method_order:
        if method in methods:
            yvalues = get_yvalues(name)(data, method)
            xvalues = xbase + offset
            clr = method_colors[str(method)]
            lab = method_labels[str(method)]
            plt.bar(xvalues, yvalues, width, color=clr, label=lab)
            offset += width

    # Label the x-axis with
    xlabels = []
    for dataset in dataset_order:
        if dataset in data:
            xlabels.append(dataset_labels[str(dataset)])

    plt.title(title, fontsize=40, ha='center', va='bottom')
    plt.ylabel(get_ylabel(name), fontsize=32, ha='center', va='bottom')
    #plt.xlabel('Sample Dataset', fontsize=32, ha='center', va='top')
    plt.xticks(xbase, xlabels, rotation=35, ha='right', fontsize=24)
    plt.legend(loc=2, fontsize=24)
    plt.savefig(filename, format='png')

#==============================================================================
# SERIAL PLOTS
#==============================================================================

make_bar_plot(1, 'throughput',
              'Serial Throughput',
              'serial_throughput.png')

make_bar_plot(1, 'duration',
              'Serial Duration',
              'serial_duration.png')

#==============================================================================
# PARALLEL PLOTS
#==============================================================================

par_title = 'Parallel Throughput (' + str(options.num_cores) + ' cores)'
make_bar_plot(options.num_cores, 'throughput',
              par_title, 'parallel_throughput.png')

par_title = 'Parallel Duration (' + str(options.num_cores) + ' cores)'
make_bar_plot(options.num_cores, 'duration',
              par_title, 'parallel_duration.png')
