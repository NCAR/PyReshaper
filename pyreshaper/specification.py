'''
This is a configuration specification class, through which the input to
the PyReshaper code is specified.  Currently all types of supported
operations for the PyReshaper are specified with derived dypes of the
Specification class.

_______________________
Created on Apr 30, 2014

@author Kevin Paul <kpaul@ucar.edu>
'''

from os import path as ospath


#==============================================================================
# create_specifier factory function
#==============================================================================
def create_specifier(spec_type='slice-to-series', **kwargs):
    '''
    Factory function for Specifier class objects.  Defined for convenience.

    @param spec_type A string specifying the type of Specifier class to
                     instantiate.  Currently, this only accepts the value
                     'slice-to-series', which is the default.

    @param kwargs    Optional arguments to be passed to the newly created
                     Specifier object's constructor.

    @return An instantiation of the type of Specifier class desired.
    '''
    # Type checking
    if (not isinstance(spec_type, str)):
        err_msg = 'Specification type must be given as a string.'
        raise TypeError(err_msg)

    # Select the right Specifier object to return
    if (spec_type == 'slice-to-series'):
        return Slice2SeriesSpecifier(**kwargs)
    else:
        err_msg = 'Specifier of type ' + str(spec_type) + ' is not defined.'
        raise ValueError(err_msg)


#==============================================================================
# Specifier Base Class
#==============================================================================
class Specifier(object):
    '''
    This is the base class for the PyReshaper input specification.  The
    base class has no functionality, but it defines the type of specifier
    that is created.
    '''

    def __init__(self):
        '''
        Constructor
        '''

        ## String specifier type
        self.specifier_type = 'undetermined'

    def validate(self):
        '''
        Perform self-validation of internal data
        '''

        # Validate types
        self.validate_types()

        # Validate values
        self.validate_values()

    def validate_types(self):
        '''
        Perform type-checking on internal data
        '''
        pass

    def validate_values(self):
        '''
        Perform value validation of internal data.

        @note We impose the (somewhat arbitrary) rule that the Specifier
              should not validate values what require "cracking" open the
              input files themselves.  Hence, we validate values that can
              be checked without any PyNIO file I/O (including reading the
              header information).
        '''
        pass


#==============================================================================
# Input Specification Class for Time-Slice to Time-Series
#==============================================================================
class Slice2SeriesSpecifier(Specifier):
    '''
    This class acts as a container for the various input data needed
    by the Reshaper to perform the time-slice to time-series operation.

    In the future, this should be made into a base class for input
    specifications, with a derived class for each type of Reshaper
    operation.  Currently, since only the time slice-to-series operation
    is defined, there is no need for more than one class.
    '''

    def __init__(self, infiles=[],
                 ncfmt='netcdf4c',
                 prefix='tseries.',
                 suffix='.nc',
                 metadata=[]):
        '''
        Initializes the internal data with optional arguments.

        @param infiles  List of full-path input filenames

        @param ncfmt    String specifying the NetCDF data format
                        ('netcdf','netcdf4','netcdf4c')

        @param prefix   String specifying the full-path prefix common
                        to all time-series output files

        @param suffix   String specifying the suffix common
                        to all time-series output files

        @param metadata List of variable names specifying the
                        variables that should be included in every
                        time-series output file

        @return An instance of the Specification class

        @note The time-series output files are named according to the
              convention:

                  output_file_name = prefix + variable_name + suffix

              The output_file_name should be a full-path filename.
        '''
        ## Type string
        self.specifier_type = 'slice-to-series'

        ## The list of input (time-slice) NetCDF files (absolute paths)
        self.input_file_list = infiles

        ## The string specifying the NetCDF file format for output
        self.netcdf_format = ncfmt

        ## The common prefix to all output files (following the rule:
        #  prefix + variable_name + suffix)
        self.output_file_prefix = prefix

        ## The common suffix to all output files (following the rule:
        #  prefix + variable_name + suffix)
        self.output_file_suffix = suffix

        ## List of time-variant variables that should be included in all
        #  output files.
        self.time_variant_metadata = metadata

    def validate_types(self):
        '''
        Method for checking the types of the Specifier data.
        Called by the validate() method.
        '''

        # Validate the type of the input file list
        if (not isinstance(self.input_file_list, list)):
            err_msg = "Input file list must be a list"
            raise TypeError(err_msg)

        # Validate that each input file name is a string
        for ifile_name in self.input_file_list:
            if (not isinstance(ifile_name, str)):
                err_msg = "Input file names must be given as strings"
                raise TypeError(err_msg)

        # Validate the netcdf format string
        if (not isinstance(self.netcdf_format, str)):
            err_msg = "NetCDF format must be given as a string"
            raise TypeError(err_msg)

        # Validate the output file prefix
        if (not isinstance(self.output_file_prefix, str)):
            err_msg = "Output file prefix must be given as a string"
            raise TypeError(err_msg)

        # Validate the output file suffix
        if (not isinstance(self.output_file_suffix, str)):
            err_msg = "Output file suffix must be given as a string"
            raise TypeError(err_msg)

        # Validate the type of the time-variant metadata list
        if (not isinstance(self.time_variant_metadata, list)):
            err_msg = "Input file list must be a list"
            raise TypeError(err_msg)

        # Validate the type of each time-variant metadata variable name
        for var_name in self.time_variant_metadata:
            if (not isinstance(var_name, str)):
                err_msg = "Time-variant metadata variable names must be " + \
                          "given as strings"
                raise TypeError(err_msg)

    def validate_values(self):
        '''
        Method to validate the values of the Specifier data.
        Called by the validate() method.

        @note We impose the (somewhat arbitrary) rule that the Specifier
              should not validate values what require "cracking" open the
              input files themselves.  Hence, we validate values that can
              be checked without any PyNIO file I/O (including reading the
              header information).

        @note This method will correct some input if it is safe to do so.
        '''

        # Make sure there is at least 1 input file given
        if (len(self.input_file_list) == 0):
            err_msg = "There must be at least one input file given."
            raise ValueError(err_msg)

        # Validate that each input file exists and is a regular file
        for ifile_name in self.input_file_list:
            if (not ospath.isfile(ifile_name)):
                err_msg = "Input file " + str(ifile_name) + \
                          " is not a regular file"
                raise ValueError(err_msg)

        # Validate the value of the netcdf format string
        valid_formats = ['netcdf', 'netcdf4', 'netcdf4c']
        if (self.netcdf_format not in valid_formats):
            err_msg = "Output NetCDF file format " \
                    + str(self.netcdf_format) \
                    + " is not valid"
            raise ValueError(err_msg)

        # Validate the output file directory
        abs_output_prefix = ospath.abspath(self.output_file_prefix)
        abs_output_dir = ospath.dirname(abs_output_prefix)
        if (not ospath.isdir(abs_output_dir)):
            err_msg = "Output directory " + str(abs_output_dir) + \
                    " implied in output prefix " + \
                    str(self.output_file_prefix) + " is not valid"
            raise ValueError(err_msg)
        self.output_file_prefix = abs_output_prefix

        # Validate the output file suffix string (should end in .nc)
        if (self.output_file_suffix[-3:] != '.nc'):
            self.output_file_suffix += '.nc'
