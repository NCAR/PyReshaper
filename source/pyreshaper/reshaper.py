"""
The module containing the Reshaper class

This is the main reshaper module.  This is where the specific operations
are defined.  Currently, only one operation has been implemented (i.e.,
the time-slice to time-series operation).

Copyright 2015, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

# Built-in imports
import abc
from os import linesep, remove
from os.path import exists, isfile

# Third-party imports
import numpy
from Nio import open_file as nio_open_file
from Nio import options as nio_options
from asaptools.simplecomm import create_comm, SimpleComm
from asaptools.timekeeper import TimeKeeper
from asaptools.partition import WeightBalanced
from asaptools.vprinter import VPrinter

# PyReshaper imports
from specification import Specifier


#==============================================================================
# create_reshaper factory function
#==============================================================================
def create_reshaper(specifier, serial=False, verbosity=1,
                    skip_existing=False, overwrite=False,
                    once=False, simplecomm=None):
    """
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
        skip_existing (bool): Flag specifying whether to skip the generation
            of time-series for variables with time-series files that already
            exist.  Default is False.
        overwrite (bool): Flag specifying whether to forcefully overwrite
            output files if they already exist.  Default is False.
        once (bool): True or False, indicating whether the Reshaper should
            write all metadata to a 'once' file (separately).
        simplecomm (SimpleComm): A SimpleComm object to handle the parallel
            communication, if necessary
        compressionlevel (int): An integer specifying the level of
            compression to be used for the output files, if compression
            is allowed by the NetCDF format

    Returns:
        Reshaper: An instance of the Reshaper object requested
    """
    # Determine the type of Reshaper object to instantiate
    if isinstance(specifier, Specifier):
        return Slice2SeriesReshaper(specifier,
                                    serial=serial,
                                    verbosity=verbosity,
                                    skip_existing=skip_existing,
                                    overwrite=overwrite,
                                    once=once,
                                    simplecomm=simplecomm)
    elif isinstance(specifier, (list, tuple)):
        spec_dict = dict([(str(i), s) for (i, s) in enumerate(specifier)])
        return create_reshaper(spec_dict,
                               serial=serial,
                               verbosity=verbosity,
                               skip_existing=skip_existing,
                               overwrite=overwrite,
                               once=once,
                               simplecomm=simplecomm)
    elif isinstance(specifier, dict):
        if not all([isinstance(s, Specifier) for s in specifier.values()]):
            err_msg = 'Multiple specifiers must all be of Specifier type'
            raise TypeError(err_msg)
        return MultiSpecReshaper(specifier,
                                 serial=serial,
                                 verbosity=verbosity,
                                 skip_existing=skip_existing,
                                 overwrite=overwrite,
                                 once=once,
                                 simplecomm=simplecomm)
    else:
        err_msg = 'Specifier of type ' + str(type(specifier)) + ' is not a ' \
            + 'valid Specifier object.'
        raise TypeError(err_msg)


#==============================================================================
# _pprint_dictionary - Helper method for printing diagnostic data
#==============================================================================
def _pprint_dictionary(title, dictionary, order=None):
    """
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
    """
    # Type checking
    if (type(title) is not str):
        err_msg = 'Title must be a str type'
        raise TypeError(err_msg)
    if (not isinstance(dictionary, dict)):
        err_msg = 'Input dictionary needs to be a dictionary type'
        raise TypeError(err_msg)
    if (order is not None and not isinstance(order, list)):
        err_msg = 'Order list needs to be a list type'
        raise TypeError(err_msg)

    # Determine the print order, if present
    print_order = dictionary.keys()
    if (order is not None):
        print_order = []
        for item in order:
            if (item in dictionary):
                print_order.append(item)

    # Header line with Title
    hline = '-' * 50 + linesep
    ostr = hline + ' ' + title.upper() + ':' + linesep + hline

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
        ostr += str(dictionary[name]) + linesep
    ostr += hline

    return ostr


#==============================================================================
# Reshaper Abstract Base Class
#==============================================================================
class Reshaper(object):

    """
    Abstract base class for Reshaper objects
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def convert(self):
        """
        Method to perform the Reshaper's designated operation.
        """
        return

    @abc.abstractmethod
    def print_diagnostics(self):
        """
        Print out timing and I/O information collected up to this point
        """
        return


#==============================================================================
# Reshaper Class
#==============================================================================
class Slice2SeriesReshaper(Reshaper):

    """
    The time-slice to time-series Reshaper class

    This is the class that defines how the time-slice to time-series
    reshaping operation is to be performed.
    """

    def __init__(self, specifier, serial=False, verbosity=1,
                 skip_existing=False, overwrite=False,
                 once=False, simplecomm=None):
        """
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
            skip_existing (bool): Flag specifying whether to skip the generation
                of time-series for variables with time-series files that already
                exist.  Default is False.
            overwrite (bool): Flag specifying whether to forcefully overwrite
                output files if they already exist.  Default is False.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel
                communication, if necessary
        """

        # Type checking (or double-checking)
        if not isinstance(specifier, Specifier):
            err_msg = "Input must be given in the form of a Specifier object"
            raise TypeError(err_msg)
        if type(serial) is not bool:
            err_msg = "Serial indicator must be True or False."
            raise TypeError(err_msg)
        if type(verbosity) is not int:
            err_msg = "Verbosity level must be an integer."
            raise TypeError(err_msg)
        if type(skip_existing) is not bool:
            err_msg = "Skip_existing flag must be True or False."
            raise TypeError(err_msg)
        if type(once) is not bool:
            err_msg = "Once-file indicator must be True or False."
            raise TypeError(err_msg)
        if simplecomm is not None:
            if not isinstance(simplecomm, SimpleComm):
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

        # Validate the user input data
        self._timer.start('Specifier Validation')
        specifier.validate()
        self._timer.stop('Specifier Validation')
        if self._simplecomm.is_manager():
            self._vprint('Specifier validated', verbosity=1)

        # Store the input file names
        self._input_filenames = specifier.input_file_list

        # Initialize the dictionary of variable names for each category
        # (Keys are variable names, Values are variable sizes)
        self._time_variant_metadata = {}
        self._time_invariant_metadata = {}
        self._time_series_variables = {}

        # Setup PyNIO options (including disabling the default PreFill option)
        opt = nio_options()
        opt.PreFill = False

        # Determine the Format and CompressionLevel options
        # from the NetCDF format string in the Specifier
        if specifier.netcdf_format == 'netcdf':
            opt.Format = 'Classic'
        elif specifier.netcdf_format in ['netcdf4', 'netcdf4c']:
            opt.Format = 'NetCDF4Classic'
            opt.CompressionLevel = specifier.compression_level
        self._nio_options = opt
        if self._simplecomm.is_manager():
            self._vprint('PyNIO options set', verbosity=2)

        # Validate the input files themselves
        if self._simplecomm.is_manager():
            self._vprint('Inspecting input files...', verbosity=1)
        self._timer.start('Inspect Input Files')
        self._inspect_input_files(specifier.time_variant_metadata)
        self._timer.stop('Inspect Input Files')
        if self._simplecomm.is_manager():
            self._vprint('Input files inspected.', verbosity=1)

        # Validate the output files
        if self._simplecomm.is_manager():
            self._vprint('Checking output file status...', verbosity=1)
        self._timer.start('Check Output File Status')
        self._check_output_file_status(specifier, skip_existing, overwrite)
        self._timer.stop('Check Output File Status')
        if self._simplecomm.is_manager():
            self._vprint('Output file status checked.', verbosity=1)

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Reshaper initialized.', verbosity=1)

        # Sync before continuing..
        self._simplecomm.sync()

    def _inspect_input_files(self, metadata_names=[]):
        """
        Inspect the input data files themselves.

        We check the file contents here.

        Parameters:
            metadata_names (list): List of string names.  If a variable with
                one of these names is found, assume it is metadata.
        """

        #===== INSPECT FIRST INPUT FILE =====

        # Open first file
        ifile = nio_open_file(self._input_filenames[0])

        # Look for the 'unlimited' dimension
        try:
            self._unlimited_dim = next(dim for dim in ifile.dimensions
                                       if ifile.unlimited(dim))
        except StopIteration:
            err_msg = 'Unlimited dimension not found.'
            raise LookupError(err_msg)

        # Get the time values
        time_values = [ifile.variables[self._unlimited_dim].get_value()]

        # Get the list of variable names and missing variables
        var_names = set(ifile.variables.keys())
        missing_vars = set()

        # Categorize each variable (only looking at first file)
        for var_name, var in ifile.variables.iteritems():
            size = numpy.dtype(var.typecode()).itemsize
            size = size * numpy.prod(var.shape)
            if self._unlimited_dim not in var.dimensions:
                self._time_invariant_metadata[var_name] = size
            elif var_name in metadata_names:
                self._time_variant_metadata[var_name] = size
            else:
                self._time_series_variables[var_name] = size

        # Close the first file
        ifile.close()

        #===== INSPECT REMAINING INPUT FILES =====

        # Make a pass through remaining files and:
        # (1) Make sure it has the 'unlimited' dimension
        # (2) Make sure this dimension is truely 'unlimited'
        # (3) Check that this dimension has a corresponding variable
        # (4) Check if there are any missing variables
        # (5) Get the time values from the files
        for ifilename in self._input_filenames[1:]:
            ifile = nio_open_file(ifilename)

            # Determine the unlimited dimension
            if self._unlimited_dim not in ifile.dimensions:
                err_msg = ('Unlimited dimension not found '
                           'in file "{}"').format(ifilename)
                raise LookupError(err_msg)
            if not ifile.unlimited(self._unlimited_dim):
                err_msg = ('Dimension "{}" not unlimited in file '
                           '"{}"').format(self._unlimited_dim, ifilename)
                raise LookupError(err_msg)
            if self._unlimited_dim not in ifile.variables:
                err_msg = ('Unlimited dimension variable not found in file '
                           '"{}"').format(ifilename)
                raise LookupError(err_msg)

            # Get the time values (list of NDArrays)
            time_values.append(
                ifile.variables[self._unlimited_dim].get_value())

            # Get the missing variables
            var_names_next = set(ifile.variables.keys())
            missing_vars.update(var_names - var_names_next)

            # Close the file
            ifile.close()

        #===== CHECK FOR MISSING VARIABLES =====

        # Make sure that the list of variables in each file is the same
        if len(missing_vars) != 0:
            warning = ("WARNING: The first input file has variables that are "
                       "not in all input files:{}{}").format(linesep, '   ')
            for var in missing_vars:
                warning += ' {}'.format(var)
            self._vprint(warning, header=True, verbosity=1)

        #===== SORT INPUT FILES BY TIME =====

        # Determine the sort order based on the first time in the time values
        old_order = range(len(self._input_filenames))
        new_order = sorted(old_order, key=lambda i: time_values[i][0])

        # Re-order the list of input filenames and time values
        new_filenames = [self._input_filenames[i] for i in new_order]
        new_values = [time_values[i] for i in new_order]

        # Now, check that the largest time in each file is less than the
        # smallest time in the next file (so that the time spans of each file
        # do not overlap)
        for i in xrange(1, len(new_values)):
            if new_values[i - 1][-1] >= new_values[i][0]:
                err_msg = ('Times in input files {} and {} appear '
                           'to overlap').format(new_filenames[i - 1],
                                                new_filenames[i])
                raise ValueError(err_msg)

        # Now that this is validated, save the time values and filename in
        # the new order
        self._input_filenames = new_filenames

        #===== FINALIZING OUTPUT =====

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

    def _check_output_file_status(self, specifier,
                                  skip_existing=False, overwrite=False):
        """
        Perform validation of output data files themselves.

        We compute the output file name from the prefix and suffix, and then
        we check whether the output files exist.  By default, if the output
        file

        Parameters:
            specifier (Specifier): The reshaper specifier object

        Keyword Arguments:
            skip_existing (bool): Flag specifying whether to skip the generation
                of time-series for variables with time-series files that already
                exist.  Default is False.
            overwrite (bool): Flag specifying whether to forcefully overwrite
                output files if they already exist.  Default is False.
        """

        # Loop through the time-series variables and generate output filenames
        prefix = specifier.output_file_prefix
        suffix = specifier.output_file_suffix
        self._time_series_filenames = \
            dict([(variable, prefix + variable + suffix)
                  for variable in self._time_series_variables])

        # Find which files already exist
        existing = []
        for variable, filename in self._time_series_filenames.items():
            if isfile(filename):
                existing.append(variable)

        # Return if no existing files found
        if len(existing) == 0:
            return

        # If overwrite is enabled, delete all existing files first
        if overwrite:
            if self._simplecomm.is_manager():
                self._vprint('WARNING: Deleting existing output files for '
                             'time-series variables: {}'.format(existing),
                             verbosity=1)
            for variable in existing:
                remove(self._time_series_filenames[variable])

        # Or, if skip_existing is set, remove the existing time-series
        # variables from the list of time-series variables to convert
        elif skip_existing:
            if self._simplecomm.is_manager():
                self._vprint('WARNING: Skipping time-series variables with '
                             'existing output files: {}'.format(existing),
                             verbosity=1)
            for variable in existing:
                self._time_series_variables.pop(variable)

        # Otherwise, throw an exception if any existing output files are found
        else:
            err_msg = ("Found existing output files for time-series "
                       "variables: {}").format(existing)
            raise RuntimeError(err_msg)

    def convert(self, output_limit=0):
        """
        Method to perform the Reshaper's designated operation.

        In this case, convert a list of time-slice files to time-series files.

        Keyword Arguments:
            output_limit (int): Limit on the number of output (time-series)
                files to write during the convert() operation.  If set
                to 0, no limit is placed.  This limits the number
                of output files produced by each processor in a
                parallel run.
        """
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

        # Partition the time-series variables across all processors
        tsv_names_loc = self._simplecomm.partition(
            self._time_series_variables.items(), func=WeightBalanced(),
            involved=True)
        if output_limit > 0:
            tsv_names_loc = tsv_names_loc[0:output_limit]

        # Print partitions for all ranks
        dbg_msg = 'Local time-series variables are {}'.format(tsv_names_loc)
        self._vprint(dbg_msg, header=True, verbosity=2)

        # Reset all of the timer values (as it is possible that there are no
        # time-series variables in the local list procuded above)
        self._timer.reset('Open Output Files')
        self._timer.reset('Close Output Files')
        self._timer.reset('Open Input Files')
        self._timer.reset('Close Input Files')
        self._timer.reset('Create Time-Invariant Metadata')
        self._timer.reset('Create Time-Variant Metadata')
        self._timer.reset('Create Time-Series Variables')
        self._timer.reset('Read Time-Invariant Metadata')
        self._timer.reset('Read Time-Variant Metadata')
        self._timer.reset('Read Time-Series Variables')
        self._timer.reset('Write Time-Invariant Metadata')
        self._timer.reset('Write Time-Variant Metadata')
        self._timer.reset('Write Time-Series Variables')

        # Initialize the byte count dictionary
        self._byte_counts['Requested Data'] = 0
        self._byte_counts['Actual Data'] = 0

        # Defining a simple helper function to determine the bytes size of
        # a variable given to it, whether an NDArray or not
        def _get_bytesize(data):
            if hasattr(data, 'nbytes'):
                return data.nbytes
            else:
                return 0

        #===== LOOP OVER TIME_SERIES VARIABLES =====

        # Loop over all time-series variables
        for out_name in tsv_names_loc:

            # Determine the output file name for this variable
            out_filename = self._time_series_filenames[out_name]
            dbg_msg = 'Creating output file for variable: {}'.format(out_name)
            if out_name == 'once':
                dbg_msg = 'Creating "once" file.'
            self._vprint(dbg_msg, header=True, verbosity=1)

            # Open the output file
            # NOTE: If the output file already exists, abort!
            self._timer.start('Open Output Files')
            if exists(out_filename):
                err_msg = 'Found existing output file: {}'.format(out_filename)
                raise OSError(err_msg)
            out_file = nio_open_file(out_filename, 'w',
                                     options=self._nio_options)
            self._timer.stop('Open Output Files')

            # Start the loop over input files (i.e., time-steps)
            series_step_index = 0
            for in_filename in self._input_filenames:

                # Open the input file
                self._timer.start('Open Input Files')
                in_file = nio_open_file(in_filename, 'r')
                self._timer.stop('Open Input Files')

                # Create header info, if this is the first input file
                if in_filename == self._input_filenames[0]:

                    # Copy file attributes and dimensions to output file
                    for name, val in in_file.attributes.iteritems():
                        setattr(out_file, name, val)
                    for name, val in in_file.dimensions.iteritems():
                        if name == self._unlimited_dim:
                            out_file.create_dimension(name, None)
                        else:
                            out_file.create_dimension(name, val)

                    # Create the metadata variables
                    if not (self._use_once_file and out_name != 'once'):

                        # Time-invariant metadata variables
                        self._timer.start('Create Time-Invariant Metadata')
                        for name in self._time_invariant_metadata:
                            in_var = in_file.variables[name]
                            out_var = out_file.create_variable(
                                name, in_var.typecode(), in_var.dimensions)
                            for att_name, att_val in in_var.attributes.iteritems():
                                setattr(out_var, att_name, att_val)
                        self._timer.stop('Create Time-Invariant Metadata')

                        # Time-variant metadata variables
                        self._timer.start('Create Time-Variant Metadata')
                        for name in self._time_variant_metadata:
                            in_var = in_file.variables[name]
                            out_var = out_file.create_variable(
                                name, in_var.typecode(), in_var.dimensions)
                            for att_name, att_val in in_var.attributes.iteritems():
                                setattr(out_var, att_name, att_val)
                        self._timer.stop('Create Time-Variant Metadata')

                    # Create the time-series variable
                    if not (self._use_once_file and out_name == 'once'):

                        # Time-series variable
                        self._timer.start('Create Time-Series Variables')
                        in_var = in_file.variables[out_name]
                        out_var = out_file.create_variable(
                            out_name, in_var.typecode(), in_var.dimensions)
                        for att_name, att_val in in_var.attributes.iteritems():
                            setattr(out_var, att_name, att_val)
                        self._timer.stop('Create Time-Series Variables')

                    dbg_msg = ('Writing output file for variable: '
                               '{}').format(out_name)
                    if out_name == 'once':
                        dbg_msg = 'Writing "once" file.'
                    self._vprint(dbg_msg, header=True, verbosity=1)

                    # Copy the time-invariant metadata
                    if not (self._use_once_file and out_name != 'once'):

                        for name in self._time_invariant_metadata:
                            in_var = in_file.variables[name]
                            out_var = out_file.variables[name]
                            self._timer.start('Read Time-Invariant Metadata')
                            tmp_data = in_var.get_value()
                            self._timer.stop('Read Time-Invariant Metadata')
                            self._timer.start('Write Time-Invariant Metadata')
                            out_var.assign_value(tmp_data)
                            self._timer.stop('Write Time-Invariant Metadata')

                            requested_nbytes = _get_bytesize(tmp_data)
                            self._byte_counts[
                                'Requested Data'] += requested_nbytes
                            actual_nbytes = self.assumed_block_size \
                                * numpy.ceil(requested_nbytes / self.assumed_block_size)
                            self._byte_counts['Actual Data'] += actual_nbytes

                # Get the number of time steps in this slice file
                num_steps = in_file.dimensions[self._unlimited_dim]

                # Explicitly loop over time steps (to control memory use)
                for slice_step_index in xrange(num_steps):

                    # Copy the time-varient metadata
                    if not (self._use_once_file and out_name != 'once'):

                        for name in self._time_variant_metadata:
                            in_var = in_file.variables[name]
                            out_var = out_file.variables[name]
                            ndims = len(in_var.dimensions)
                            udidx = in_var.dimensions.index(self._unlimited_dim)
                            in_slice = [slice(None)] * ndims
                            in_slice[udidx] = slice_step_index
                            out_slice = [slice(None)] * ndims
                            out_slice[udidx] = series_step_index
                            self._timer.start('Read Time-Variant Metadata')
                            tmp_data = in_var[tuple(in_slice)]
                            self._timer.stop('Read Time-Variant Metadata')
                            self._timer.start('Write Time-Variant Metadata')
                            out_var[tuple(out_slice)] = tmp_data
                            self._timer.stop('Write Time-Variant Metadata')

                            requested_nbytes = _get_bytesize(tmp_data)
                            self._byte_counts[
                                'Requested Data'] += requested_nbytes
                            actual_nbytes = self.assumed_block_size \
                                * numpy.ceil(requested_nbytes / self.assumed_block_size)
                            self._byte_counts['Actual Data'] += actual_nbytes

                    # Copy the time-series variables
                    if not (self._use_once_file and out_name == 'once'):

                        in_var = in_file.variables[out_name]
                        out_var = out_file.variables[out_name]
                        ndims = len(in_var.dimensions)
                        udidx = in_var.dimensions.index(self._unlimited_dim)
                        in_slice = [slice(None)] * ndims
                        in_slice[udidx] = slice_step_index
                        out_slice = [slice(None)] * ndims
                        out_slice[udidx] = series_step_index
                        self._timer.start('Read Time-Series Variables')
                        tmp_data = in_var[tuple(in_slice)]
                        self._timer.stop('Read Time-Series Variables')
                        self._timer.start('Write Time-Series Variables')
                        out_var[tuple(out_slice)] = tmp_data
                        self._timer.stop('Write Time-Series Variables')

                        requested_nbytes = _get_bytesize(tmp_data)
                        self._byte_counts['Requested Data'] += requested_nbytes
                        actual_nbytes = self.assumed_block_size \
                            * numpy.ceil(requested_nbytes / self.assumed_block_size)
                        self._byte_counts['Actual Data'] += actual_nbytes

                    # Increment the time-series step index
                    series_step_index += 1

                # Close the input file
                self._timer.start('Close Input Files')
                in_file.close()
                self._timer.stop('Close Input Files')

            # Close the output file
            self._timer.start('Close Output Files')
            out_file.close()
            self._timer.stop('Close Output Files')

            # Output message to user
            dbg_msg = 'Closed output file for variable: {}'.format(out_name)
            if out_name == 'once':
                dbg_msg = 'Closed "once" file.'
            self._vprint(dbg_msg, header=True, verbosity=1)

        # Information
        self._simplecomm.sync()
        if self._simplecomm.is_manager():
            self._vprint(('Finished converting time-slices '
                          'to time-series.'), verbosity=1)

        # Finish clocking the entire convert procedure
        self._timer.stop('Complete Conversion Process')

    def print_diagnostics(self):
        """
        Print out timing and I/O information collected up to this point
        """

        # Get all totals and maxima
        my_times = self._timer.get_all_times()
        max_times = self._simplecomm.allreduce(my_times, op='max')
        my_bytes = self._byte_counts
        total_bytes = self._simplecomm.allreduce(my_bytes, op='sum')

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
# MultiSpecReshaper Class
#==============================================================================
class MultiSpecReshaper(Reshaper):

    """
    Multiple Slice-to-Series Reshaper class

    This class is designed to deal with lists of multiple
    Slice2SeriesSpecifiers at a time.  Instead of being instantiated
    (or initialized) with a single Slice2SeriesSpecifier,
    it takes a dictionary of Slice2SeriesSpecifier objects.
    """

    def __init__(self, specifiers, serial=False, verbosity=1,
                 skip_existing=False, overwrite=False,
                 once=False, simplecomm=None):
        """
        Constructor

        Parameters:
            specifiers (dict): A dict of named Specifier instances, each
                defining an input specification for this reshaper operation.

        Keyword Arguments:
            serial (bool): True or False, indicating whether the operation
                should be performed in serial (True) or parallel
                (False).  The default is to assume parallel operation
                (but serial will be chosen if the mpi4py cannot be
                found when trying to initialize decomposition.
            verbosity(int): Level of printed output (stdout).  A value of 0
                means no output, and a higher value means more output.  The
                default value is 1.
            skip_existing (bool): Flag specifying whether to skip the generation
                of time-series for variables with time-series files that already
                exist.  Default is False.
            overwrite (bool): Flag specifying whether to forcefully overwrite
                output files if they already exist.  Default is False.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel
                communication, if necessary
        """

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
        if type(skip_existing) is not bool:
            err_msg = "Skip_existing flag must be True or False."
            raise TypeError(err_msg)
        if type(once) is not bool:
            err_msg = "Once-file indicator must be True or False."
            raise TypeError(err_msg)
        if simplecomm is not None:
            if not isinstance(simplecomm, SimpleComm):
                err_msg = "Simple communicator object is not a SimpleComm"
                raise TypeError(err_msg)

        # Whether to write to a once file
        self._use_once_file = once

        # Whether to write to a once file
        self._skip_existing = skip_existing

        # Whether to write to overwrite output files
        self._overwrite = overwrite

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

    def convert(self, output_limit=0):
        """
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
        """
        # Type checking input
        if type(output_limit) is not int:
            err_msg = 'Output limit must be an integer'
            raise TypeError(err_msg)

        # Loop over all specifiers
        for spec_name in self._specifiers:
            if self._simplecomm.is_manager():
                self._vprint('--- Converting Specifier: ' +
                             str(spec_name), verbosity=0)

            rshpr = create_reshaper(self._specifiers[spec_name],
                                    serial=self._serial,
                                    verbosity=self._verbosity,
                                    skip_existing=self._skip_existing,
                                    overwrite=self._overwrite,
                                    once=self._use_once_file,
                                    simplecomm=self._simplecomm)
            rshpr.convert(output_limit=output_limit)

            this_times = rshpr._timer.get_all_times()
            self._times[spec_name] = rshpr._simplecomm.allreduce(
                this_times, op='max')
            self._time_orders[spec_name] = rshpr._timer.get_names()
            this_count = rshpr._byte_counts
            self._byte_counts[spec_name] = rshpr._simplecomm.allreduce(
                this_count, op='sum')

            if self._simplecomm.is_manager():
                self._vprint('--- Finished converting Specifier: ' +
                             str(spec_name) + linesep, verbosity=0)
            self._simplecomm.sync()

    def print_diagnostics(self):
        """
        Print out timing and I/O information collected up to this point
        """
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
