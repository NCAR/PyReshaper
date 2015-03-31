'''
The module containing the Reshaper class

This is the main reshaper module.  This is where the specific operations
are defined.  Currently, only one operation has been implemented (i.e.,
the time-slice to time-series operation).

'''

from specification import Specifier, Slice2SeriesSpecifier
from asaptools.simplecomm import create_comm, SimpleComm
from asaptools.timekeeper import TimeKeeper
from asaptools.partition import WeightBalanced
from asaptools.vprinter import VPrinter

import abc
import os
import Nio
import numpy
import itertools


#==============================================================================
# create_reshaper factory function
#==============================================================================
def create_reshaper(specifier, serial=False, verbosity=1,
                    once=False, simplecomm=None):
    '''
    Factory function for Reshaper class instantiations.

    Parameters:
        specifier (Specifier): An instantiation of a Specifier class that 
            defines the type of operation to be performed.  (That is, a
            Slice2SeriesSpecifier will invoke the creation of a
            matching Slice2SeriesReshaper object.)
            Alternatively, this can be a list of Specifier objects,
            or dictionary of named Specifier objects.
            In this case, a reshaper will be created for each
            specifier in the list, and each reshaper will be
            created and run in sequence.

    Keyword Arguments:
        serial (bool): True or False, indicating whether the Reshaper object
            should perform its operation in serial (True) or
            parallel (False).
        verbosity (int): Level of printed output (stdout).  A value of 0 means
            no output, and a higher value means more output.  The
            default value is 1.
        once (bool): True or False, indicating whether the Reshaper should
            write all metadata to a 'once' file (separately).
        simplecomm (SimpleComm): A SimpleComm object to handle the parallel
            communication, if necessary

    Returns:
        Reshaper: An instance of the Reshaper object requested
    '''
    # Determine the type of Reshaper object to instantiate
    if type(specifier) is Slice2SeriesSpecifier:
        return Slice2SeriesReshaper(specifier,
                                    serial=serial,
                                    verbosity=verbosity,
                                    once=once,
                                    simplecomm=simplecomm)
    elif isinstance(specifier, list):
        spec_dict = dict([(str(i), s) for (i, s) in enumerate(specifier)])
        return create_reshaper(spec_dict,
                               serial=serial,
                               verbosity=verbosity,
                               once=once,
                               simplecomm=simplecomm)
    elif isinstance(specifier, dict):
        spec_types = set([type(s) for s in specifier.values()])
        if len(spec_types) > 1:
            err_msg = 'Multiple specifiers must all have the same type'
            raise TypeError(err_msg)
        spec_type = spec_types.pop()
        if spec_type is Slice2SeriesSpecifier:
            return MultiSpecS2SReshaper(specifier,
                                        serial=serial,
                                        verbosity=verbosity,
                                        once=once,
                                        simplecomm=simplecomm)
        else:
            err_msg = 'Multiple specifiers of type ' + str(spec_type) \
                + ' are not valid.'
            raise TypeError(err_msg)
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

    Parameters:
        title (str): The title to give to the printed table
        dictionary (dict): A dictionary of numeric values

    Keyword Arguments:
        order (list): The print order for the keys in the dictionary (only 
            items that are in both the order list and the dictionary will be
            printed)

    Return:
        str: A string with the pretty-printed dictionary data
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
    The Reshaper abstract base class  

    Specific operations should be defined as derived types.
    '''
    __metaclass__ = abc.ABCMeta

    def __init__(self, specifier, serial=False,
                 verbosity=1, once=False, simplecomm=None):
        '''
        Constructor

        Parameters:
            specifier (Specifier): An instance of the Specifier class, 
                defining the input specification for this reshaper operation.

        Keyword Arguments:
            serial (bool): True or False, indicating whether the operation
                should be performed in serial (True) or parallel
                (False).  The default is to assume parallel operation
                (but serial will be chosen if the mpi4py cannot be
                found when trying to initialize decomposition.
            verbosity(int): Level of printed output (stdout).  A value of 0 
                means no output, and a higher value means more output.  The
                default value is 1.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel 
                communication, if necessary
        '''

        # Check types
        if not isinstance(specifier, Specifier):
            err_msg = "Input must be given in the form of a Specifier object"
            raise TypeError(err_msg)
        if type(serial) is not bool:
            err_msg = "Serial indicator must be True or False."
            raise TypeError(err_msg)
        if type(verbosity) is not int:
            err_msg = "Verbosity level must be an integer."
            raise TypeError(err_msg)
        if type(once) is not bool:
            err_msg = "Once-file indicator must be True or False."
            raise TypeError(err_msg)
        if simplecomm is not None:
            if simplecomm is not isinstance(simplecomm, SimpleComm):
                err_msg = "Simple communicator object is not a SimpleComm"
                raise TypeError(err_msg)

        # Whether to write a once file
        self._use_once_file = once

        # Internal timer data
        self._timer = TimeKeeper()

        # Dictionary storing read/write data amounts
        self.assumed_block_size = float(4 * 1024 * 1024)
        self._byte_counts = {}

        self._timer.start('Initializing Simple Communicator')
        if simplecomm is None:
            simplecomm = create_comm(serial=serial)
        # Reference to the simple communicator
        self._simplecomm = simplecomm
        self._timer.stop('Initializing Simple Communicator')

        # Contruct the print header
        header = ''.join(['[', str(self._simplecomm.get_rank()),
                          '/', str(self._simplecomm.get_size()), '] '])

        # Reference to the verbose printer tool
        self._vprint = VPrinter(header=header, verbosity=verbosity)

        # Debug output starting
        if self._simplecomm.is_manager():
            self._vprint('Initializing Reshaper', verbosity=1)

        # Input specification
        self._specifier = specifier

        # Validate the user input data
        self._timer.start('Specifier Validation')
        self._specifier.validate()
        self._timer.stop('Specifier Validation')
        if self._simplecomm.is_manager():
            self._vprint('Specifier validated', verbosity=1)

    @abc.abstractmethod
    def convert(self, output_limit=0):
        '''
        Method to perform the Reshaper's designated operation.

        Keyword Arguments:
            output_limit (int): Limit on the number of output (time-series) 
                files to write during the convert() operation.  If set
                to 0, no limit is placed.  This limits the number
                of output files produced by each processor in a
                parallel run.
        '''
        return

    def print_diagnostics(self):
        '''
        Print out timing and I/O information collected up to this point
        '''

        # Get all totals and maxima
        my_times = self._timer.get_all_times()
        # print str(self._simplecomm.get_rank()) + ': all_times =', my_times
        max_times = self._simplecomm.allreduce(my_times, op='max')
        # print str(self._simplecomm.get_rank()) + ': max_times =', max_times
        my_bytes = self._byte_counts
        # print str(self._simplecomm.get_rank()) + ': byte_counts =', my_bytes
        total_bytes = self._simplecomm.allreduce(my_bytes, op='sum')
        # print str(self._simplecomm.get_rank()) + ': total_bytes =',
        # total_bytes

        # Synchronize
        self._simplecomm.sync()

        # Print timing maxima
        o = self._timer.get_names()
        time_table_str = _pprint_dictionary('TIMING DATA', max_times, order=o)
        if self._simplecomm.is_manager():
            self._vprint(time_table_str, verbosity=0)

        # Convert byte count to MB
        for name in total_bytes:
            total_bytes[name] = total_bytes[name] / float(1024 * 1024)

        # Print byte count totals
        byte_count_str = _pprint_dictionary('BYTE COUNTS (MB)', total_bytes)
        if self._simplecomm.is_manager():
            self._vprint(byte_count_str, verbosity=0)


#==============================================================================
# Slice2SeriesReshaper Class
#==============================================================================
class Slice2SeriesReshaper(Reshaper):

    '''
    The time-slice to time-series Reshaper class

    This is the class that defines how the time-slice to time-series 
    reshaping operation is to be performed.
    '''

    def __init__(self, specifier, serial=False,
                 verbosity=1, once=False, simplecomm=None):
        '''
        Constructor

        Parameters:
            specifier (Specifier): An instance of the Specifier class, 
                defining the input specification for this reshaper operation.

        Keyword Arguments:
            serial (bool): True or False, indicating whether the operation
                should be performed in serial (True) or parallel
                (False).  The default is to assume parallel operation
                (but serial will be chosen if the mpi4py cannot be
                found when trying to initialize decomposition.
            verbosity(int): Level of printed output (stdout).  A value of 0 
                means no output, and a higher value means more output.  The
                default value is 1.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel 
                communication, if necessary
        '''

        # Type checking (or double-checking)
        if not isinstance(specifier, Slice2SeriesSpecifier):
            err_msg = "Slice2SeriesReshaper requires a Slice2SeriesSpecifier" \
                + " as input."
            raise TypeError(err_msg)

        # Call the base-class constructor
        super(Slice2SeriesReshaper, self).__init__(specifier,
                                                   serial=serial,
                                                   verbosity=verbosity,
                                                   once=once,
                                                   simplecomm=simplecomm)

        # Setup PyNIO options (including disabling the default PreFill option)
        opt = Nio.options()
        opt.PreFill = False

        # Determine the Format and CompressionLevel options
        # from the NetCDF format string in the Specifier
        if self._specifier.netcdf_format == 'netcdf':
            opt.Format = 'Classic'
        elif self._specifier.netcdf_format == 'netcdf4':
            opt.Format = 'NetCDF4Classic'
            opt.CompressionLevel = 0
        elif self._specifier.netcdf_format == 'netcdf4c':
            opt.Format = 'NetCDF4Classic'
            opt.CompressionLevel = 1
        self._nio_options = opt
        if self._simplecomm.is_manager():
            self._vprint('PyNIO options set', verbosity=2)

        # Open all of the input files
        self._timer.start('Open Input Files')
        self._input_files = []
        for filename in self._specifier.input_file_list:
            self._input_files.append(Nio.open_file(filename, "r"))
        self._timer.stop('Open Input Files')
        if self._simplecomm.is_manager():
            self._vprint('Input files opened', verbosity=2)

        # Validate the input files themselves
        self._timer.start('Input File Validation')
        self._validate_input_files()
        self._timer.stop('Input File Validation')
        if self._simplecomm.is_manager():
            self._vprint('Input files validated', verbosity=2)

        # Sort the input files by time
        self._timer.start('Sort Input Files')
        self._sort_input_files_by_time()
        self._timer.stop('Sort Input Files')
        if self._simplecomm.is_manager():
            self._vprint('Input files sorted', verbosity=2)

        # Retrieve and sort the variables in each time-slice file
        # (To determine if it is time-invariant metadata, time-variant
        # metadata, or if it is a time-series variable)
        self._timer.start('Sort Variables')
        self._sort_variables()
        self._timer.stop('Sort Variables')
        if self._simplecomm.is_manager():
            self._vprint('Variables sorted', verbosity=2)

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Reshaper initialized.', verbosity=1)

        # Sync before continuing..
        self._simplecomm.sync()

    def _validate_input_files(self):
        '''
        Perform validation of input data files themselves.  

        We check the file contents here, assuming that the files are already 
        open.
        '''

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Validating input files', verbosity=1)

        # In the first file, look for the 'unlimited' dimension
        ifile = self._input_files[0]
        self._unlimited_dim = None
        for dim in ifile.dimensions:
            if ifile.unlimited(dim):
                self._unlimited_dim = dim
                continue  # There can only be 1!
        if self._unlimited_dim == None:
            err_msg = 'Unlimited dimension not identified.'
            raise LookupError(err_msg)

        # Make a pass through each file and:
        # (1) Make sure it has the 'unlimited' dimension
        # (2) Make sure this dimension is truely 'unlimited'
        # (3) Check that this dimension has a corresponding variable
        for i in range(len(self._input_files)):
            ifile = self._input_files[i]
            if self._unlimited_dim not in ifile.dimensions:
                err_msg = 'Unlimited dimension not found in file (' \
                    + self._specifier.input_file_list[i] + ')'
                raise LookupError(err_msg)
            if not ifile.unlimited(self._unlimited_dim):
                err_msg = 'Unlimited dimension not unlimited in file (' \
                    + self._specifier.input_file_list[i] + ')'
                raise LookupError(err_msg)
            if self._unlimited_dim not in ifile.variables:
                err_msg = 'Unlimited dimension variable not found in file (' \
                    + self._specifier.input_file_list[i] + ')'
                raise LookupError(err_msg)

        # Make sure that the list of variables in each file is the same
        variables = self._input_files[0].variables
        var_names = set(variables.keys())
        missing_vars = set()
        for ifile in self._input_files[1:]:
            var_names_next = set(ifile.variables.keys())
            missing_vars.update(var_names - var_names_next)
        if len(missing_vars) != 0:
            warning = "WARNING: The first input file has variables that are " \
                + "not in all input files:" + os.linesep + '   '
            for var in missing_vars:
                warning += ' ' + str(var)
            self._vprint(warning, header=True, verbosity=1)

    def _sort_input_files_by_time(self):
        '''
        Internal method for sorting the input files by time

        This assumes that 'time' is the unlimited dimension, and it checks
        to make sure that all of the times spanning across each file do not 
        overlap with each other (i.e., that the times across all files are 
        monotonicly increasing).

        Currently, this method assumes that all of the input files
        have the same 'time:units' attribute, such that all time variable
        values are measured from the same date-time.  When this is true,
        we do not need to consider the value of the 'time:units'
        attribute itself.  If this assumption is not true, then we need
        to consider the 'time:units" attribute of each file, together
        with that file's time variable values.  To do that properly,
        however, one should use UDUNITS to do the comparisons.
        '''

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Sorting input files', verbosity=1)

        # Get the time attributes (for convenience) and, for each file,
        # add the times to a list.  (Each file will have an array of times
        # associated with it.  Each array will be added to a list, such
        # that the outer-most list contains an array for each input file)
        time_values = []
        for ifile in self._input_files:
            time_values.append(
                ifile.variables[self._unlimited_dim].get_value())

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
            if new_values[i][-1] >= new_values[i + 1][0]:
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
        Internal method for sorting the variables in each time-slice file

        This method determines if each variable is to be treated as 
        time-invariant metadata, time-variant metadata (user defined), or 
        time-series variables.  All metadata is written to every time-series 
        file, and any time-series variable is written to its own file.  
        The time-variant metadata variables are determined by user input, 
        and are contained in the Specifier data member:

            Specifier.time_variant_metadata.
        '''

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Sorting variables', verbosity=1)

        # Initialize the dictionary of variable names for each category
        # (Keys are variable names, Values are variable sizes)
        self._time_variant_metadata = {}
        self._time_invariant_metadata = {}
        self._time_series_variables = {}

        # Categorize each variable (only looking at first file)
        variables = self._input_files[0].variables
        for var_name in variables.keys():
            var = variables[var_name]
            size = numpy.dtype(var.typecode()).itemsize
            size = size * numpy.prod(var.shape)
            if self._unlimited_dim not in var.dimensions:
                self._time_invariant_metadata[var_name] = size
            elif var_name in self._specifier.time_variant_metadata:
                self._time_variant_metadata[var_name] = size
            else:
                self._time_series_variables[var_name] = size

        # Debug output
        if self._simplecomm.is_manager():
            self._vprint('Time-Invariant Metadata: ' +
                         str(self._time_invariant_metadata.keys()), verbosity=2)
            self._vprint('Time-Variant Metadata: ' +
                         str(self._time_variant_metadata.keys()), verbosity=2)
            self._vprint('Time-Series Variables: ' +
                         str(self._time_series_variables.keys()), verbosity=2)

        # Add 'once' variable if writing to a once file
        # NOTE: This is a "cheat"!  There is no 'once' variable.  It's just
        #       a catch for all metadata IFF the 'once-file' is enabled.
        if self._use_once_file:
            self._time_series_variables['once'] = 1

    def convert(self, output_limit=0):
        '''
        Method to perform the Reshaper's designated operation.

        In this case, convert a list of time-slice files to time-series files.

        Keyword Arguments:
            output_limit (int): Limit on the number of output (time-series) 
                files to write during the convert() operation.  If set
                to 0, no limit is placed.  This limits the number
                of output files produced by each processor in a
                parallel run.
        '''
        # Type checking input
        if type(output_limit) is not int:
            err_msg = 'Output limit must be an integer'
            raise TypeError(err_msg)

        # Start the total convert process timer
        self._simplecomm.sync()
        self._timer.start('Complete Conversion Process')

        # Debugging output
        if self._simplecomm.is_manager():
            self._vprint('Converting time-slices to time-series', verbosity=1)

        # For data common to all input files, we reference only the first
        ref_infile = self._input_files[0]

        # Store the common dimensions and attributes for each file
        # (taken from the first input file in the list)
        common_dims = ref_infile.dimensions
        common_atts = ref_infile.attributes

        # Partition the time-series variables across all processors
        tsv_names_loc = self._simplecomm.partition(self._time_series_variables.items(),
                                                   func=WeightBalanced(),
                                                   involved=True)
        if output_limit > 0:
            tsv_names_loc = tsv_names_loc[0:output_limit]

        # Print partitions for all ranks
        dbg_msg = 'Local time-series variables are ' + str(tsv_names_loc)
        self._vprint(dbg_msg, header=True, verbosity=2)

        # Reset all of the timer values (as it is possible that there are no
        # time-series variables in the local list procuded above)
        self._timer.reset('Open Output Files')
        self._timer.reset('Create Time-Invariant Metadata')
        self._timer.reset('Create Time-Variant Metadata')
        self._timer.reset('Create Time-Series Variables')
        self._timer.reset('Write Time-Invariant Metadata')
        self._timer.reset('Write Time-Variant Metadata')
        self._timer.reset('Write Time-Series Variables')
        self._timer.reset('Close Output Files')

        # Initialize the byte count dictionary
        self._byte_counts['Requested Data'] = 0
        self._byte_counts['Actual Data'] = 0

        # Defining a simple helper function to determine whether to
        # write time-series data and/or write metadata.  This is useful
        # for adding the ability to write a "once" file
        def _get_once_info(vname):
            is_once_file = (vname == 'once')
            write_meta = True
            write_tser = True
            if self._use_once_file:
                write_meta = is_once_file
                write_tser = not is_once_file
            return is_once_file, write_meta, write_tser

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
            is_once_file, write_meta, write_tser = _get_once_info(out_name)

            # Determine the output file name for this variable
            out_filename = self._specifier.output_file_prefix \
                + out_name + self._specifier.output_file_suffix
            dbg_msg = 'Creating output file for variable: ' + out_name
            if is_once_file:
                dbg_msg = 'Creating "once" file.'
            self._vprint(dbg_msg, header=True, verbosity=1)

            # Open each output file and create the dimensions and attributes
            # NOTE: If the output file already exists, abort!
            self._timer.start('Open Output Files')
            if os.path.exists(out_filename):
                err_msg = 'Found existing output file: ' + out_filename
                raise OSError(err_msg)
            out_file = Nio.open_file(out_filename, 'w',
                                     options=self._nio_options)
            for att_name, att_val in common_atts.iteritems():
                setattr(out_file, att_name, att_val)
            for dim_name, dim_val in common_dims.iteritems():
                if dim_name == self._unlimited_dim:
                    out_file.create_dimension(dim_name, None)
                else:
                    out_file.create_dimension(dim_name, dim_val)
            self._timer.stop('Open Output Files')

            # Create the time-invariant metadata variables
            if (write_meta):
                self._timer.start('Create Time-Invariant Metadata')
                for name in self._time_invariant_metadata:
                    in_var = ref_infile.variables[name]
                    out_var = out_file.create_variable(name,
                                                       in_var.typecode(),
                                                       in_var.dimensions)
                    for att_name, att_val in in_var.attributes.iteritems():
                        setattr(out_var, att_name, att_val)
                self._timer.stop('Create Time-Invariant Metadata')

            # Create the time-variant metadata variables
            if write_meta:
                self._timer.start('Create Time-Variant Metadata')
                for name in self._time_variant_metadata:
                    in_var = ref_infile.variables[name]
                    out_tvm_vars[name] = out_file.create_variable(name,
                                                                  in_var.typecode(), in_var.dimensions)
                    for att_name, att_val in in_var.attributes.iteritems():
                        setattr(out_tvm_vars[name], att_name, att_val)
                self._timer.stop('Create Time-Variant Metadata')

            # Create the time-series variable itself
            if write_tser:
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
            is_once_file, write_meta, write_tser = _get_once_info(out_name)

            dbg_msg = 'Writing output file for variable: ' + out_name
            if is_once_file:
                dbg_msg = 'Writing "once" file.'
            self._vprint(dbg_msg, header=True, verbosity=1)

            # Create the attributes of the time-series variable
            if write_tser:
                in_var = ref_infile.variables[out_name]
                out_var = out_file.variables[out_name]
                for att_name, att_val in in_var.attributes.iteritems():
                    setattr(out_var, att_name, att_val)

            # Write the time-invariant metadata
            if write_meta:
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
            series_step_index = 0
            for in_file in self._input_files:

                # Get the number of time steps in this slice file
                num_steps = in_file.dimensions[self._unlimited_dim]

                # Loop over the time steps in this slice file
                for slice_step_index in range(num_steps):

                    # Write the time-varient metadata
                    if write_meta:
                        self._timer.start('Write Time-Variant Metadata')
                        for name in self._time_variant_metadata:
                            in_meta = in_file.variables[name]
                            out_meta = out_file.variables[name]
                            ndims = len(in_meta.dimensions)
                            udidx = in_meta.dimensions.index(
                                self._unlimited_dim)
                            in_slice = [slice(None)] * ndims
                            in_slice[udidx] = slice_step_index
                            out_slice = [slice(None)] * ndims
                            out_slice[udidx] = series_step_index
                            out_meta[tuple(out_slice)] = in_meta[
                                tuple(in_slice)]

                            requested_nbytes = in_meta[:].nbytes
                            self._byte_counts[
                                'Requested Data'] += requested_nbytes
                            actual_nbytes = self.assumed_block_size \
                                * numpy.ceil(requested_nbytes / self.assumed_block_size)
                            self._byte_counts['Actual Data'] += actual_nbytes
                        self._timer.stop('Write Time-Variant Metadata')

                    # Write the time-series variables
                    if write_tser:
                        self._timer.start('Write Time-Series Variables')
                        in_var = in_file.variables[out_name]
                        ndims = len(in_var.dimensions)
                        udidx = in_var.dimensions.index(self._unlimited_dim)
                        in_slice = [slice(None)] * ndims
                        in_slice[udidx] = slice_step_index
                        out_slice = [slice(None)] * ndims
                        out_slice[udidx] = series_step_index
                        out_var[tuple(out_slice)] = in_var[tuple(in_slice)]

                        requested_nbytes = in_file.variables[
                            out_name][:].nbytes
                        self._byte_counts['Requested Data'] += requested_nbytes
                        actual_nbytes = self.assumed_block_size \
                            * numpy.ceil(requested_nbytes / self.assumed_block_size)
                        self._byte_counts['Actual Data'] += actual_nbytes
                        self._timer.stop('Write Time-Series Variables')

                    # Increment the time-series step index
                    series_step_index += 1

            # Close the output file
            self._timer.start('Close Output Files')
            out_file.close()
            self._timer.stop('Close Output Files')
            dbg_msg = 'Closed output file for variable: ' + out_name
            if is_once_file:
                dbg_msg = 'Closed "once" file.'
            self._vprint(dbg_msg, header=True, verbosity=1)

        # Information
        self._simplecomm.sync()
        if self._simplecomm.is_manager():
            self._vprint(
                'Finished converting time-slices to time-series.', verbosity=1)

        # Finish clocking the entire convert procedure
        self._timer.stop('Complete Conversion Process')


#==============================================================================
# MultiSpecReshaper Base Class
#==============================================================================
class MultiSpecReshaper(object):

    '''
    Multiple-Specification Reshaper abstract base class

    This class is designed to deal with dictionaries of multiple named 
    Specifiers at a time.  Instead of being instantiated (or initialized) with
    a single Specifier, it takes a dictionary of named Specifier objects.

    This class does not inherit from the base Reshaper class.  This is really
    more like a container + wrapper class around the base Reshaper class.
    '''
    __metaclass__ = abc.ABCMeta

    def __init__(self, specifiers, serial=False,
                 verbosity=1, once=False, simplecomm=None):
        '''
        Constructor

        Parameters:
            specifier (Specifier): An instance of the Specifier class, 
                defining the input specification for this reshaper operation.

        Keyword Arguments:
            serial (bool): True or False, indicating whether the operation
                should be performed in serial (True) or parallel
                (False).  The default is to assume parallel operation
                (but serial will be chosen if the mpi4py cannot be
                found when trying to initialize decomposition.
            verbosity(int): Level of printed output (stdout).  A value of 0 
                means no output, and a higher value means more output.  The
                default value is 1.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel 
                communication, if necessary
        '''

        # Check types
        if not isinstance(specifiers, dict):
            err_msg = "Input must be given in a dictionary of Specifiers"
            raise TypeError(err_msg)
        if type(serial) is not bool:
            err_msg = "Serial indicator must be True or False."
            raise TypeError(err_msg)
        if type(verbosity) is not int:
            err_msg = "Verbosity level must be an integer."
            raise TypeError(err_msg)
        if simplecomm is not None:
            if simplecomm is not isinstance(simplecomm, SimpleComm):
                err_msg = "Simple communicator object is not a SimpleComm"
                raise TypeError(err_msg)

        # Whether to write to a once file
        self._use_once_file = once

        # Store the list of specifiers
        self._specifiers = specifiers

        # Store the serial specifier
        self._serial = serial

        # Check for a SimpleComm, and if none create it
        if simplecomm is None:
            simplecomm = create_comm(serial=serial)

        # Pointer to its own messenger
        self._simplecomm = simplecomm

        # Store the verbosity
        self._verbosity = verbosity

        # Set the verbose printer
        self._vprint = VPrinter(verbosity=verbosity)

        # Storage for timing data
        self._times = {}

        # Orders for printing timing data
        self._time_orders = {}

        # Storage for all byte counters
        self._byte_counts = {}

    @abc.abstractmethod
    def convert(self, output_limit):
        '''
        Method to perform the Reshaper's designated operation.

        Keyword Arguments:
            output_limit (int): Limit on the number of output (time-series) 
                files to write during the convert() operation.  If set
                to 0, no limit is placed.  This limits the number
                of output files produced by each processor in a
                parallel run.
        '''
        return

    def print_diagnostics(self):
        '''
        Print out timing and I/O information collected up to this point
        '''
        # Loop through all timers
        for name in self._specifiers:
            if self._simplecomm.is_manager():
                self._vprint('Specifier: ' + str(name), verbosity=0)

            times = self._times[name]
            o = self._time_orders[name]
            times_str = _pprint_dictionary('TIMING DATA', times, order=o)
            if self._simplecomm.is_manager():
                self._vprint(times_str, verbosity=0)

            counts = self._byte_counts[name]
            for name in counts:
                counts[name] = counts[name] / float(1024 * 1024)
            counts_str = _pprint_dictionary('BYTE COUNTS (MB)', counts)
            if self._simplecomm.is_manager():
                self._vprint(counts_str, verbosity=0)


#==============================================================================
# MultiSpecS2SReshaper Class
#==============================================================================
class MultiSpecS2SReshaper(MultiSpecReshaper):

    '''
    Multiple Slice-to-Series Reshaper class

    This class is designed to deal with lists of multiple 
    Slice2SeriesSpecifiers at a time.  Instead of being instantiated 
    (or initialized) with a single Slice2SeriesSpecifier,
    it takes a dictionary of Slice2SeriesSpecifier objects.
    '''

    def __init__(self, specifiers, serial=False,
                 verbosity=1, once=False, simplecomm=None):
        '''
        Constructor

        Parameters:
            specifier (Specifier): An instance of the Specifier class, 
                defining the input specification for this reshaper operation.

        Keyword Arguments:
            serial (bool): True or False, indicating whether the operation
                should be performed in serial (True) or parallel
                (False).  The default is to assume parallel operation
                (but serial will be chosen if the mpi4py cannot be
                found when trying to initialize decomposition.
            verbosity(int): Level of printed output (stdout).  A value of 0 
                means no output, and a higher value means more output.  The
                default value is 1.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel 
                communication, if necessary
        '''

        # Call the base class
        super(MultiSpecS2SReshaper, self).__init__(specifiers,
                                                   serial=serial,
                                                   verbosity=verbosity,
                                                   once=once,
                                                   simplecomm=simplecomm)

    def convert(self, output_limit=0):
        '''
        Method to perform each Reshaper's designated operation.

        Loops through and creates each Reshaper, calls each Reshaper's 
        convert() method, and pulls the timing data out for each convert 
        operation.

        Keyword Arguments:
            output_limit (int): Limit on the number of output (time-series) 
                files to write during the convert() operation.  If set
                to 0, no limit is placed.  This limits the number
                of output files produced by each processor in a
                parallel run.
        '''
        # Type checking input
        if type(output_limit) is not int:
            err_msg = 'Output limit must be an integer'
            raise TypeError(err_msg)

        # Loop over all specifiers
        for spec_name in self._specifiers:
            if self._simplecomm.is_manager():
                self._vprint('--- Converting Specifier: '
                             + str(spec_name), verbosity=0)

            rshpr = create_reshaper(self._specifiers[spec_name],
                                    serial=self._serial,
                                    verbosity=self._verbosity,
                                    once=self._use_once_file)
            rshpr.convert(output_limit=output_limit)

            this_times = rshpr._timer.get_all_times()
            self._times[spec_name] = rshpr._simplecomm.allreduce(
                this_times, op='max')
            self._time_orders[spec_name] = rshpr._timer.get_names()
            this_count = rshpr._byte_counts
            self._byte_counts[spec_name] = rshpr._simplecomm.allreduce(
                this_count, op='sum')

            if self._simplecomm.is_manager():
                self._vprint('--- Finished converting Specifier: '
                             + str(spec_name) + os.linesep, verbosity=0)
            self._simplecomm.sync()
