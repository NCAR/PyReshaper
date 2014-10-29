'''
This is the main reshaper module.  This is where the specific operations
are defined.  Currently, only one operation has been implemented (i.e.,
the time-slice to time-series operation).

________________________
Created on Apr 30, 2014

@author Kevin Paul <kpaul@ucar.edu>
'''

from specification import Specifier, Slice2SeriesSpecifier
from pyreshaper.messenger import create_messenger
from timekeeper import TimeKeeper

import os
import Nio
import numpy
import itertools


#==============================================================================
# create_reshaper factory function
#==============================================================================
def create_reshaper(specifier, serial=False, verbosity=1):
    '''
    Factory function for Reshaper class instantiations.

    @param specifier  An instantiation of a Specifier class that defines
                      the type of operation to be performed.  (That is, a
                      Slice2SeriesSpecifier will invoke the creation of a
                      matching Slice2SeriesReshaper object.)

                      Alternatively, this can be a list of Specifier objects,
                      or dictionary of named Specifier objects.
                      In this case, a reshaper will be created for each
                      specifier in the list, and each reshaper will be
                      created and run in sequence.

    @param serial     True or False, indicating whether the Reshaper object
                      should perform its operation in serial (True) or
                      parallel (False).

    @param verbosity  Level of printed output (stdout).  A value of 0 means
                      no output, and a higher value means more output.  The
                      default value is 1.

    @return  An instantiation of a Reshaper class that matches the type of
             the Specifier object given as input
    '''
    # Type checking
    if (type(serial) is not bool):
        err_msg = "Serial indicator must be True or False."
        raise TypeError(err_msg)

    # Determine the type of Reshaper object to instantiate
    if (type(specifier) is Slice2SeriesSpecifier):
        return Slice2SeriesReshaper(specifier,
                                    serial=serial,
                                    verbosity=verbosity)
    elif (isinstance(specifier, list)):
        spec_dict = {}
        for i in range(len(specifier)):
            spec_dict[i] = specifier[i]
        return MultiSpecReshaper(spec_dict,
                                 serial=serial,
                                 verbosity=verbosity)
    elif (isinstance(specifier, dict)):
        return MultiSpecReshaper(specifier,
                                 serial=serial,
                                 verbosity=verbosity)
    else:
        err_msg = 'Specifier of type ' + str(type(specifier)) + ' is not a ' \
                + 'valid Specifier object.'
        raise TypeError(err_msg)


#==============================================================================
# _pprint_dictionary - Helper method for printing diagnostic data
#==============================================================================
def _pprint_dictionary(title, dictionary, order=None):
    '''
    Hidden method for pretty-printing a dictionary of numeric values,
    with a given title.

    @param title The title to give to the printed table

    @param dictionary Dictionary of numeric values

    @param order The print order for the keys in the dictionary (only items
                 that are in both the order list and the dictionary will be
                 printed)

    @return A string with the pretty-printed dictionary data
    '''
    # Type checking
    if (type(title) is not str):
        err_msg = 'Title must be a str type'
        raise TypeError(err_msg)
    if (not isinstance(dictionary, dict)):
        err_msg = 'Input dictionary needs to be a dictionary type'
        raise TypeError(err_msg)
    if (order != None and not isinstance(order, list)):
        err_msg = 'Order list needs to be a list type'
        raise TypeError(err_msg)

    # Determine the print order, if present
    print_order = dictionary.keys()
    if (order != None):
        print_order = []
        for item in order:
            if (item in dictionary):
                print_order.append(item)

    # Header line with Title
    hline = '-' * 50 + os.linesep
    ostr = hline + ' ' + title.upper() + ':' + os.linesep + hline

    # Determine the longest timer name
    # and thus computer the column to line up values
    valcol = 0
    for name in print_order:
        if (len(str(name)) > valcol):
            valcol = len(str(name))
    valcol += 2

    # Print out the timer names and the accumulated times
    for name in print_order:
        spacer = ' ' * (valcol - len(str(name)))
        ostr += str(name) + ':' + spacer
        ostr += str(dictionary[name]) + os.linesep
    ostr += hline

    return ostr


#==============================================================================
# Reshaper Base Class
#==============================================================================
class Reshaper(object):
    '''
    Base Class for the PyReshaper general operations.  Specific operations
    should be defined as derived types.
    '''

    def __init__(self, specifier, serial=False, verbosity=1):
        '''
        Constructor

        @param specifier  An instance of the Specifier class, defining the
                          input specification for this reshaper operation.

        @param serial     True or False, indicating whether the operation
                          should be performed in serial (True) or parallel
                          (False).  The default is to assume parallel operation
                          (but serial will be chosen if the mpi4py cannot be
                          found when trying to initialize decomposition.

        @param verbosity  Level of printed output (stdout).  A value of 0 means
                          no output, and a higher value means more output.  The
                          default value is 1.
        '''

        # Check types
        if (not isinstance(specifier, Specifier)):
            err_msg = "Input must be given in the form of a Specifier object"
            raise TypeError(err_msg)
        if (type(serial) is not bool):
            err_msg = "Serial indicator must be True or False."
            raise TypeError(err_msg)
        if (type(verbosity) is not int):
            err_msg = "Verbosity level must be an integer."
            raise TypeError(err_msg)

        ## Internal timer data
        self._timer = TimeKeeper()

        ## Dictionary storing read/write data amounts
        self.assumed_block_size = float(4 * 1024 * 1024)
        self._byte_counts = {}

        self._timer.start('Initializing Messenger')
        ## Pointer to the messenger
        self._messenger = create_messenger(serial=serial)
        self._timer.stop('Initializing Messenger')

        # Set the verbosity in the messenger
        self._messenger.verbosity = verbosity

        # Debug output starting
        self._messenger.print_once('Initializing Reshaper', vlevel=1)

        ## Input specification
        self._specifier = specifier

        # Validate the user input data
        self._timer.start('Specifier Validation')
        self._specifier.validate()
        self._timer.stop('Specifier Validation')
        self._messenger.print_once('Specifier validated', vlevel=1)

    def convert(self, output_limit=0):
        '''
        Method to perform the Reshaper's designated operation.

        @param output_limit Limit on the number of output files to write
                            during the convert() operation.  If set to 0,
                            no limit is placed.  This limits the number of
                            output files produced by each processor in a
                            parallel run.
        '''
        pass

    def print_diagnostics(self):
        '''
        Print out timing and I/O information collected up to this point
        in the Reshaper's operation.
        '''
        # Get all totals and maxima
        max_times = self._messenger.max(self._timer.get_all_times())
        total_bytes = self._messenger.sum(self._byte_counts)

        # Synchronize
        self._messenger.sync()

        # Print timing maxima
        o = self._timer.get_order()
        time_table_str = _pprint_dictionary('TIMING DATA', max_times, order=o)
        self._messenger.print_once(time_table_str, vlevel=0)

        # Convert byte count to MB
        for name in total_bytes:
            total_bytes[name] = total_bytes[name] / float(1024 * 1024)

        # Print byte count totals
        byte_count_str = _pprint_dictionary('BYTE COUNTS (MB)', total_bytes)
        self._messenger.print_once(byte_count_str, vlevel=0)


#==============================================================================
# Slice2SeriesReshaper Class
#==============================================================================
class Slice2SeriesReshaper(Reshaper):
    '''
    This is the main reshaper class.  This is the class that defines how
    the time-slice to time-series reshaping operation is to be performed.
    '''

    def __init__(self, spec, serial=False, verbosity=1):
        '''
        Constructor

        @param spec  An instance of the Specifier class, defining the
                          input specification for this reshaper operation.

        @param serial     True or False, indicating whether the operation
                          should be performed in serial (True) or parallel
                          (False).  The default is to assume parallel operation
                          (but serial will be chosen if the mpi4py cannot be
                          found when trying to initialize decomposition.

        @param verbosity  Level of printed output (stdout).  A value of 0 means
                          no output, and a higher value means more output.  The
                          default value is 1.
        '''
        # Type checking (or double-checking)
        if (not isinstance(spec, Slice2SeriesSpecifier)):
            err_msg = "Slice2SeriesReshaper requires a Slice2SeriesSpecifier" \
                    + " as input."
            raise TypeError(err_msg)

        # Call the base-class constructor
        super(Slice2SeriesReshaper, self).__init__(spec,
                                                   serial=serial,
                                                   verbosity=verbosity)

        # Setup PyNIO options (including disabling the default PreFill option)
        opt = Nio.options()
        opt.PreFill = False

        # Determine the Format and CompressionLevel options
        # from the NetCDF format string in the Specifier
        if (self._specifier.netcdf_format == 'netcdf'):
            opt.Format = 'Classic'
        elif (self._specifier.netcdf_format == 'netcdf4'):
            opt.Format = 'NetCDF4Classic'
            opt.CompressionLevel = 0
        elif (self._specifier.netcdf_format == 'netcdf4c'):
            opt.Format = 'NetCDF4Classic'
            opt.CompressionLevel = 1
        self._nio_options = opt
        self._messenger.print_once('PyNIO options set', vlevel=2)

        # Open all of the input files
        self._timer.start('Open Input Files')
        self._input_files = []
        for filename in self._specifier.input_file_list:
            self._input_files.append(Nio.open_file(filename, "r"))
        self._timer.stop('Open Input Files')
        self._messenger.print_once('Input files opened', vlevel=2)

        # Validate the input files themselves
        self._timer.start('Input File Validation')
        self._validate_input_files()
        self._timer.stop('Input File Validation')
        self._messenger.print_once('Input files validated', vlevel=2)

        # Sort the input files by time
        self._timer.start('Sort Input Files')
        self._sort_input_files_by_time()
        self._timer.stop('Sort Input Files')
        self._messenger.print_once('Input files sorted', vlevel=2)

        # Retrieve and sort the variables in each time-slice file
        # (To determine if it is time-invariant metadata, time-variant
        # metadata, or if it is a time-series variable)
        self._timer.start('Sort Variables')
        self._sort_variables()
        self._timer.stop('Sort Variables')
        self._messenger.print_once('Variables sorted', vlevel=2)

        # Helpful debugging message
        self._messenger.print_once('Reshaper initialized.', vlevel=1)

    def _validate_input_files(self):
        '''
        Perform validation of input data files themselves.  We check
        the file contents here, assuming that the files are already open.
        '''

        # Helpful debugging message
        self._messenger.print_once('Validating input files', vlevel=1)

        # Make a pass through each file to make sure there is a 'time'
        # dimension and a 'time' variable
        for i in range(len(self._input_files)):
            ifile = self._input_files[i]
            if ('time' not in ifile.dimensions):
                err_msg = 'Time dimension not found in file (' \
                        + self._specifier.input_file_list[i] + ')'
                raise LookupError(err_msg)
            if ('time' not in ifile.variables):
                err_msg = 'Time variable not found in file (' \
                        + self._specifier.input_file_list[i] + ')'
                raise LookupError(err_msg)

        # Make sure that the list of variables in each file is the same
        variables = self._input_files[0].variables
        var_names = set(variables.keys())
        missing_vars = set()
        for ifile in self._input_files[1:]:
            var_names_next = set(ifile.variables.keys())
            missing_vars.update(var_names - var_names_next)
        if (len(missing_vars) != 0):
            print "WARNING: The first input file has variables that are " \
                + "not in all input files:"
            warning = ' '
            for var in missing_vars:
                warning += ' ' + str(var)
            print warning

    def _sort_input_files_by_time(self):
        '''
        Internal method for sorting the input files by time, and check to
        make sure that all of the times spanning across each file do not
        overlap with each other (i.e., that the times across all files are
        monotonicly increasing).  We also store the array of times across
        all input files for future writing.

        @note Currently, this method assumes that all of the input files
              have the same 'time:units' attribute, such that all time variable
              values are measured from the same date-time.  When this is true,
              we do not need to consider the value of the 'time:units'
              attribute itself.  If this assumption is not true, then we need
              to consider the 'time:units" attribute of each file, together
              with that file's time variable values.  This is what the CFTime
              class is intended to do.
        '''

        # Helpful debugging message
        self._messenger.print_once('Sorting input files', vlevel=1)

        # Get the time attributes (for convenience) and, for each file,
        # add the times to a list.  (Each file will have an array of times
        # associated with it.  Each array will be added to a list, such
        # that the outer-most list contains an array for each input file)
        time_values = []
        for ifile in self._input_files:
            time_values.append(ifile.variables['time'].get_value())

        # Determine the sort order based on the first time in the time values
        order = range(len(self._input_files))
        new_order = sorted(order, key=lambda i: time_values[i][0])

        # Re-order the list of input files and filenames
        new_file_list = [None] * len(new_order)
        new_filenames = [None] * len(new_order)
        new_values = [None] * len(new_order)
        for i in order:
            new_file_list[i] = self._input_files[new_order[i]]
            new_filenames[i] = self._specifier.input_file_list[new_order[i]]
            new_values[i] = time_values[new_order[i]]

        # Save this data in the new orders
        self._input_files = new_file_list
        self._input_filenames = new_filenames

        # Now, check that the largest time in each file is less than the
        # smallest time in the next file (so that the time spans of each file
        # do not overlap)
        for i in order[:-1]:
            if (new_values[i][-1] >= new_values[i + 1][0]):
                err_msg = 'Times in input files ' + str(new_filenames[i]) \
                        + ' and ' + str(new_filenames[i + 1]) + ' appear to ' \
                        + 'overlap.'
                raise ValueError(err_msg)

        # Now that this is validated, let's string together the numpy array
        # of all times (using the new_values array)
        self._all_time_values = \
            numpy.fromiter(itertools.chain.from_iterable(new_values),
                           dtype='float')

    def _sort_variables(self):
        '''
        Internal method for sorting the variables that exist in each
        time-slice file.  This method determines if each variable is to be
        treated as time-invariant metadata, time-variant metadata (user
        defined), or time-series variables.  All metadata is written to
        every time-series file, and any time-series variable is written to
        its own file.  The time-variant metadata variables are determined
        by user input, and are contained in the Specifier data member,
        Specifier.time_variant_metadata.
        '''

        # Helpful debugging message
        self._messenger.print_once('Sorting variables', vlevel=1)

        # Initialize the list of variable names for each category
        self._time_variant_metadata = []
        self._time_invariant_metadata = []
        self._time_series_variables = []

        # Categorize each variable (only looking at first file)
        variables = self._input_files[0].variables
        for var_name in self._input_files[0].variables.keys():
            if ('time' not in variables[var_name].dimensions):
                self._time_invariant_metadata.append(var_name)
            elif (var_name in self._specifier.time_variant_metadata):
                self._time_variant_metadata.append(var_name)
            else:
                self._time_series_variables.append(var_name)

        # Debug output
        self._messenger.print_once('Time-Invariant Metadata: ' + \
                                  str(self._time_invariant_metadata), vlevel=2)
        self._messenger.print_once('Time-Variant Metadata: ' + \
                                  str(self._time_variant_metadata), vlevel=2)
        self._messenger.print_once('Time-Series Variables: ' + \
                                  str(self._time_series_variables), vlevel=2)

    def convert(self, output_limit=0):
        '''
        Method to perform the Reshaper's designated operation.  In this case,
        convert a list of time-slice files to time-series.

        @param output_limit Limit on the number of output (time-series) files
                            to write during the convert() operation.  If set
                            to 0, no limit is placed.  This limits the number
                            of output files produced by each processor in a
                            parallel run.
        '''
        # Type checking input
        if (type(output_limit) is not int):
            err_msg = 'Output limit must be an integer'
            raise TypeError(err_msg)

        # Start the total convert process timer
        self._timer.start('Complete Conversion Process')

        # Debugging output
        self._messenger.sync()
        self._messenger.print_once('Converting time-slices to time-series',
                                   vlevel=1)

        # For data common to all input files, we reference only the first
        ref_infile = self._input_files[0]

        # Store the common dimensions and attributes for each file
        # (taken from the first input file in the list)
        common_dims = ref_infile.dimensions
        common_atts = ref_infile.attributes

        # Partition the time-series variables across all processors
        tsv_names_loc = self._messenger.partition(self._time_series_variables)
        if (output_limit > 0):
            tsv_names_loc = tsv_names_loc[0:output_limit]

        # Print partitions for all ranks
        dbg_msg = 'Local time-series variables are ' \
                + str(tsv_names_loc)
        self._messenger.print_all(dbg_msg, vlevel=2)

        # Initialize the byte count dictionary
        self._byte_counts['Requested Data'] = 0
        self._byte_counts['Actual Data'] = 0

        # NOTE: In the prototype, we check for the existance of the output
        # directory at this point.  If it does not exist, we create it (but
        # only from the master rank).  This requires synchronization with
        # the decomp utility.  Instead, we assume the output directory
        # already exists (and is checked by the Specifier's validation).  No
        # synchronization is needed.

        # For each time-series variable, create the corresponding output file
        # (Also defines the header info for each output file)
        out_files = {}
        out_tvm_vars = {}
        for out_name in tsv_names_loc:
            out_filename = self._specifier.output_file_prefix \
                         + out_name + self._specifier.output_file_suffix
            self._messenger.print_all('Creating output file for variable: ' + \
                                      out_name, vlevel=1)

            # Open each output file and create the dimensions and attributes
            # NOTE: If the output file already exists, abort!
            self._timer.start('Open Output Files')
            if (os.path.exists(out_filename)):
                err_msg = 'Found existing output file: ' + out_filename
                raise OSError(err_msg)
            out_file = Nio.open_file(out_filename, 'w',
                                     options=self._nio_options)
            for att_name, att_val in common_atts.iteritems():
                setattr(out_file, att_name, att_val)
            for dim_name, dim_val in common_dims.iteritems():
                if (dim_name == 'time'):
                    out_file.create_dimension(dim_name, None)
                else:
                    out_file.create_dimension(dim_name, dim_val)
            self._timer.stop('Open Output Files')

            # Create the time-invariant metadata variables
            self._timer.start('Create Time-Invariant Metadata')
            for name in self._time_invariant_metadata:
                in_var = ref_infile.variables[name]
                out_var = out_file.create_variable(name,
                    in_var.typecode(), in_var.dimensions)
                for att_name, att_val in in_var.attributes.iteritems():
                    setattr(out_var, att_name, att_val)
            self._timer.stop('Create Time-Invariant Metadata')

            # Create the time-variant metadata variables
            self._timer.start('Create Time-Variant Metadata')
            for name in self._time_variant_metadata:
                in_var = ref_infile.variables[name]
                out_tvm_vars[name] = out_file.create_variable(name,
                    in_var.typecode(), in_var.dimensions)
                for att_name, att_val in in_var.attributes.iteritems():
                    setattr(out_tvm_vars[name], att_name, att_val)
            self._timer.stop('Create Time-Variant Metadata')

            # Create the time-series variable itself
            self._timer.start('Create Time-Series Variables')
            in_var = ref_infile.variables[out_name]
            out_var = out_file.create_variable(out_name,
                in_var.typecode(), in_var.dimensions)
            self._timer.stop('Create Time-Series Variables')

            # Append the output file to list
            out_files[out_name] = out_file

        # Now that each output file has been created, start writing the data
        # (Looping over output file index, which is common in name lists)
        for out_name, out_file in out_files.iteritems():

            dbg_msg = 'Writing output file for variable: ' + out_name
            self._messenger.print_all(dbg_msg, vlevel=1)

            # Create the attributes of the time-series variable
            in_var = ref_infile.variables[out_name]
            out_var = out_file.variables[out_name]
            for att_name, att_val in in_var.attributes.iteritems():
                setattr(out_var, att_name, att_val)

            # Write the time-invariant metadata
            self._timer.start('Write Time-Invariant Metadata')
            for name in self._time_invariant_metadata:
                in_meta = ref_infile.variables[name]
                out_meta = out_file.variables[name]
                if in_meta.rank > 0:
                    out_meta[:] = in_meta[:]
                else:
                    out_meta.assign_value(in_meta.get_value())
            self._timer.stop('Write Time-Invariant Metadata')

            # Write each time-variant variable
            t_index = 0
            for in_file in self._input_files:

                # Write the time-varient metadata
                self._timer.start('Write Time-Variant Metadata')
                for name in self._time_variant_metadata:
                    in_meta = in_file.variables[name]
                    out_meta = out_file.variables[name]
                    out_meta[t_index] = in_meta[:]

                    requested_nbytes = in_meta[:].nbytes
                    self._byte_counts['Requested Data'] += requested_nbytes
                    actual_nbytes = self.assumed_block_size \
                      * numpy.ceil(requested_nbytes / self.assumed_block_size)
                    self._byte_counts['Actual Data'] += actual_nbytes

                self._timer.stop('Write Time-Variant Metadata')

                # Write the time-series variables
                self._timer.start('Write Time-Series Variables')
                out_var[t_index] = in_file.variables[out_name][:]
                requested_nbytes = in_file.variables[out_name][:].nbytes
                self._byte_counts['Requested Data'] += requested_nbytes
                actual_nbytes = self.assumed_block_size \
                    * numpy.ceil(requested_nbytes / self.assumed_block_size)
                self._byte_counts['Actual Data'] += actual_nbytes
                self._timer.stop('Write Time-Series Variables')

                # Increment the time index
                t_index += 1

            # Close the output file
            self._timer.start('Close Output Files')
            out_file.close()
            self._timer.stop('Close Output Files')
            self._messenger.print_all('Closed output file for variable: ' + \
                                      out_name, vlevel=1)

        # Information
        self._messenger.sync()
        self._messenger.print_once(
            'Finished converting time-slices to time-series.', vlevel=1)

        # Finish clocking the entire convert procedure
        self._timer.stop('Complete Conversion Process')


#==============================================================================
# MultiSpecReshaper Class
#==============================================================================
class MultiSpecReshaper(object):
    '''
    This class is designed to deal with lists of multiple Specifiers at a time.
    Instead of being instantiated (or initialized) with a single Specifier,
    it takes a list of Specifier objects.

    @note This class does not inherit from the base Reshaper class.  This is
          really more like a container + wrapper class around the base Reshaper
          class.
    '''

    def __init__(self, specifiers, serial=False, verbosity=1):
        '''
        Constructor

        @param specifiers A dictionary of instances of the Specifier class,
                          defining the input specifications for each reshaper
                          operation.  The keys unquely label each Specifier
                          with string name, and the values of the dictionary
                          are the Specifier objects themselves.

        @param serial     True or False, indicating whether the operations
                          should be performed in serial (True) or parallel
                          (False).  The default is to assume parallel operation
                          (but serial will be chosen if the mpi4py cannot be
                          found when trying to initialize decomposition.

        @param verbosity  Level of printed output (stdout).  A value of 0 means
                          no output, and a higher value means more output.  The
                          default value is 1.
        '''

        # Check types
        if (not isinstance(specifiers, dict)):
            err_msg = "Input must be given in a dictionary of Specifiers"
            raise TypeError(err_msg)
        if (type(serial) is not bool):
            err_msg = "Serial indicator must be True or False."
            raise TypeError(err_msg)
        if (type(verbosity) is not int):
            err_msg = "Verbosity level must be an integer."
            raise TypeError(err_msg)

        ## Store the list of specifiers
        self._specifiers = specifiers

        ## Store the serial specifier
        self._serial = serial

        ## Pointer to its own messenger
        self._messenger = create_messenger(serial=serial)

        ## Store the verbosity
        self._verbosity = verbosity

        # Set the messenger's verbosity
        self._messenger.verbosity = verbosity

        ## Storage for timing data
        self._times = {}

        ## Orders for printing timing data
        self._time_orders = {}

        ## Storage for all byte counters
        self._byte_counts = {}

    def convert(self, output_limit=0):
        '''
        Method to perform each Reshaper's designated operation.  Loops through
        and creates each Reshaper, calls each Reshaper's convert() method,
        and pulls the timing data out for each convert operation.

        @param output_limit Limit on the number of output (time-series) files
                            to write during the convert() operation.  If set
                            to 0, no limit is placed.  This limits the number
                            of output files produced by each processor in a
                            parallel run.
        '''
        # Type checking input
        if (type(output_limit) is not int):
            err_msg = 'Output limit must be an integer'
            raise TypeError(err_msg)

        # Loop over all specifiers
        for spec_name in self._specifiers:
            self._messenger.print_once('--- Converting Specifier: ' \
                + spec_name, vlevel=0)

            rshpr = create_reshaper(self._specifiers[spec_name],
                                    serial=self._serial,
                                    verbosity=self._verbosity)
            rshpr.convert(output_limit=output_limit)

            this_times = rshpr._timer.get_all_times()
            self._times[spec_name] = rshpr._messenger.max(this_times)
            self._time_orders[spec_name] = rshpr._timer.get_order()
            this_count = rshpr._byte_counts
            self._byte_counts[spec_name] = rshpr._messenger.sum(this_count)

            self._messenger.print_once('--- Finished converting Specifier: ' \
                + spec_name + os.linesep, vlevel=0)
            self._messenger.sync()

    def print_diagnostics(self):
        '''
        Print out timing and I/O information collected up to this point
        in all of the contained timers.
        '''
        # Loop through all timers
        for name in self._specifiers:
            self._messenger.print_once('Specifier: ' + str(name), vlevel=0)

            times = self._times[name]
            o = self._time_orders[name]
            times_str = _pprint_dictionary('TIMING DATA', times, order=o)
            self._messenger.print_once(times_str, vlevel=0)

            counts = self._byte_counts[name]
            for name in counts:
                counts[name] = counts[name] / float(1024 * 1024)
            counts_str = _pprint_dictionary('BYTE COUNTS (MB)', counts)
            self._messenger.print_once(counts_str, vlevel=0)
