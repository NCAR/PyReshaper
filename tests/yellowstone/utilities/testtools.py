#==============================================================================
#
#  TestTools
#
#  This is a collection of functions that are useful for running the PyReshaper
#  tests on the Yellowstone compute system.
#
#==============================================================================

# Builtin Modules
import os
import glob
import json

# Third-Party Modules
import numpy as np
import Nio

# Package Modules
from pyreshaper import specification


#==============================================================================
# Private Bytesize from Typecode Calculator
#==============================================================================
def _bytesize(tc):
    DTYPE_MAP = {'d': np.float64, 'f': np.float32, 'l': np.long, 'i': np.int32,
                 'h': np.int16, 'b': np.int8, 'S1': np.character}
    return np.dtype(DTYPE_MAP.get(tc, np.float)).itemsize


#==============================================================================
# Private Size from Shape Calculator
#==============================================================================
def _shape2size(shp):
    return 1 if len(shp) < 1 else reduce(lambda x, y: x * y, shp)


#==============================================================================
# Private Bytesize to Unit-string Converter
#==============================================================================
def _nbyte_str(n, exp=0):
    BYTE_UNITS = ['Bytes', 'KB', 'MB', 'GB', 'PB']
    if (n > 1024.):
        return _nbyte_str(n / 1024., exp=exp + 1)
    else:
        if exp < len(BYTE_UNITS):
            units = BYTE_UNITS[exp]
        else:
            n *= 1024.**(exp + 1 - len(BYTE_UNITS))
            units = BYTE_UNITS[-1]
        return '%.4f %s' % (n, units)


#==============================================================================
# Test Database Class
#==============================================================================
class TestDB(object):

    def __init__(self, file_name=None):
        """
        Constructor

        Parameters:
            file_name (str): The name of the test database file.  Defaults
                to 'testinfo.json'.

        Raises:
            ValueError: If the test database file cannot be opened and/or
                read.
        """
        # See if there is a user-defined testinfo file,
        # otherwise look for default
        abs_path = ''
        if file_name:
            abs_path = os.path.abspath(file_name)
        else:
            this_dir = os.path.dirname(__file__)
            abs_path = os.path.join(this_dir, 'testinfo.json')

        # Try opening and reading the testinfo file
        self._database = {}
        try:
            dbfile = open(file_name, 'r')
            self._database = dict(json.load(dbfile))
            dbfile.close()
        except:
            err_msg = 'Problem reading and parsing test info file: ' \
                + str(abs_path)
            raise ValueError(err_msg)

        # Reset the analysis state to False
        self._analyzed = False

    def list(self):
        """
        List the tests in the test database.
        """
        print 'Tests found in the Test Database are:'
        print
        for test_name in self._database:
            print '   ' + str(test_name)
        return

    def create_specifier(self, test_name, ncfmt='netcdf4c', **kwargs):
        """
        Create a Specifier object for the given named test.

        Parameters:
            test_name (str): The string name of the test in the database
                for which to construct the Specifier.
            ncfmt (str): The NetCDF format string to be passed to the
                Specifier.
            kwargs (dict): A dictionary of additional options to be
                sent to the Specifier.

        Returns:
            Specifier: A Specifier instance with the information to run the
                named test.
        """

        # Check types
        if type(test_name) is not str:
            err_msg = "Test name must be a string"
            raise TypeError(err_msg)
        if type(ncfmt) is not str:
            err_msg = "NetCDF format must be a string"
            raise TypeError(err_msg)

        # Check for the given test name
        if test_name not in self._database:
            err_msg = "Test '" + test_name + "' not found in database"
            raise ValueError(err_msg)

        # Define the Specifier input
        input_dir = str(self._database[test_name]['input_dir'])
        infiles = []
        for input_glob in self._database[test_name]['input_globs']:
            input_glob_str = str(input_glob)
            full_input_glob = str(os.path.join(input_dir, input_glob))
            infiles.extend(glob.glob(full_input_glob))
        prefix = str(self._database[test_name]['output_prefix'])
        suffix = str(self._database[test_name]['output_suffix'])
        metadata = map(str, self._database[test_name]['metadata'])

        # Remove duplicate parameters from the kwargs dictionary
        kwargs.pop('infiles', None)
        kwargs.pop('prefix', None)
        kwargs.pop('suffix', None)
        kwargs.pop('metadata', None)

        # Return the Specifier
        return specification.Specifier(infiles=infiles, ncfmt=ncfmt,
                                       prefix=prefix, suffix=suffix,
                                       metadata=metadata, **kwargs)

    def create_specifiers(self, test_list=[], ncfmt='netcdf4c', **kwargs):
        """
        Create a dictionary of named Specifier objects for a list of tests.

        Parameters:
            test_list (list): A list of string names of tests in the database
                for which to construct Specifiers.
            ncfmt (str): The NetCDF format string to be passed to the
                Specifier.
            kwargs (dict): A dictionary of additional options to be
                sent to the Specifier.

        Returns:
            Specifier: A Specifier instance with the information to run the
                named test.
        """

        # Check types
        if type(test_list) is not list:
            err_msg = "Test list must be a list of string test names"
            raise TypeError(err_msg)
        test_types = [type(test_name) is str for test_name in test_list]
        if not all(test_types):
            err_msg = "Test list must be a list of string test names"
            raise TypeError(err_msg)

        # If test list is empty, assume all tests
        if len(test_list) == 0:
            test_list = self._database.keys()

        # Construct the list of specifiers
        spec_list = [self.create_specifier(test_name, ncfmt, **kwargs)
                     for test_name in test_list]

        # Return a dictionary of named specifiers
        return dict(zip(test_list, spec_list))

    def analyze(self):
        """
        Analyze the test database to determine test statistics
        """

        # If analysis has already been done, skip
        if self._analyzed:
            return

        # Initialize test statistics database
        self._statistics = {}

        # Generate statistics for each test
        for test_name in self._database:

            # Create a specifier for this test
            spec = self.create_specifier(str(test_name), ncfmt='netcdf')

            # Validate the test information
            spec.validate()

            # Sort the input files by name
            spec.input_file_list.sort()

            # Open the first input file
            infile = Nio.open_file(spec.input_file_list[0], 'r')

            # Get the name of the unlimited dimension (e.g., time)
            tdim = None
            for dim in infile.dimensions:
                if infile.unlimited(dim):
                    tdim = dim
                    continue

            # Get the data dimensions
            metadata_names = set(infile.dimensions.keys())

            # Add the extra metadata variable names
            metadata_names.update(set(spec.time_variant_metadata))

            # Gather statistics for variables in dataset
            self._statistics[test_name] = {}
            self._statistics[test_name]['length'] = infile.dimensions[tdim]
            self._statistics[test_name]['variables'] = {}
            for var_name in infile.variables.keys():
                self._statistics[test_name]['variables'][var_name] = {}
                var_obj = infile.variables[var_name]

                tvariant = False
                xshape = var_obj.shape
                if tdim in var_obj.dimensions:
                    xshape = list(xshape)
                    tindex = var_obj.dimensions.index(tdim)
                    xshape.pop(tindex)
                    xshape = tuple(xshape)
                    tvariant = True
                self._statistics[test_name]['variables'][
                    var_name]['tvariant'] = tvariant
                self._statistics[test_name]['variables'][
                    var_name]['xshape'] = xshape

                xlen = _shape2size(xshape)
                xsize = xlen * _bytesize(var_obj.typecode())
                self._statistics[test_name]['variables'][
                    var_name]['xsize'] = xsize

                if var_name in metadata_names:
                    self._statistics[test_name][
                        'variables'][var_name]['meta'] = True
                else:
                    self._statistics[test_name]['variables'][
                        var_name]['meta'] = False

            # Close the first file
            infile.close()

            # Loop over all input files and compute time-variant variable sizes
            for file_name in spec.input_file_list[1:]:

                # Open the file
                infile = Nio.open_file(file_name, 'r')

                # And number of time steps to the test data
                self._statistics[test_name][
                    'length'] += infile.dimensions[tdim]

                # Close the file
                infile.close()

            # Compute self-analysis parameters
            num_steps = self._statistics[test_name]['length']
            var_stats = self._statistics[test_name]['variables']
            num_vars = len(var_stats)
            meta_mask = [var_stats[v]['meta'] for v in var_stats]
            tvar_mask = [var_stats[v]['tvariant'] for v in var_stats]

            timd_mask = [meta_mask[i] and not tvar_mask[i]
                         for i in range(num_vars)]
            tvmd_mask = [meta_mask[i] and tvar_mask[i]
                         for i in range(num_vars)]
            tser_mask = [not meta_mask[i] and tvar_mask[i]
                         for i in range(num_vars)]
            lost_mask = [not meta_mask[i] and not tvar_mask[i]
                         for i in range(num_vars)]

            # Compute numbers/counts
            self._statistics[test_name]['counts'] = {}
            self._statistics[test_name]['counts'][
                'tseries'] = sum(tser_mask)
            self._statistics[test_name]['counts'][
                'tvariant'] = sum(tvmd_mask)
            self._statistics[test_name]['counts'][
                'tinvariant'] = sum(timd_mask)
            self._statistics[test_name]['counts'][
                'other'] = sum(lost_mask)

            # Compute shapes
            self._statistics[test_name]['xshapes'] = {}
            self._statistics[test_name]['xshapes'][
                'tseries'] = set([d['xshape'] for d, m in
                                  zip(var_stats.values(), tser_mask) if m])
            self._statistics[test_name]['xshapes'][
                'tvariant'] = set([d['xshape'] for d, m in
                                   zip(var_stats.values(), tvmd_mask) if m])
            self._statistics[test_name]['xshapes'][
                'tinvariant'] = set([d['xshape'] for d, m in
                                     zip(var_stats.values(), timd_mask) if m])

            # Compute bytesizes
            self._statistics[test_name]['totalsizes'] = {}
            self._statistics[test_name]['totalsizes'][
                'tseries'] = sum([d['xsize'] for d, m in
                                  zip(var_stats.values(), tser_mask) if m]) * num_steps
            self._statistics[test_name]['totalsizes'][
                'tvariant'] = sum([d['xsize'] for d, m in
                                   zip(var_stats.values(), tvmd_mask) if m]) * num_steps
            self._statistics[test_name]['totalsizes'][
                'tinvariant'] = sum([d['xsize'] for d, m in
                                     zip(var_stats.values(), timd_mask) if m])

            # Compute maxima
            self._statistics[test_name]['maxsizes'] = {}
            self._statistics[test_name]['maxsizes'][
                'tseries'] = max([d['xsize'] for d, m in
                                  zip(var_stats.values(), tser_mask) if m]) * num_steps
            self._statistics[test_name]['maxsizes'][
                'tvariant'] = max([d['xsize'] for d, m in
                                   zip(var_stats.values(), tvmd_mask) if m]) * num_steps
            self._statistics[test_name]['maxsizes'][
                'tinvariant'] = max([d['xsize'] for d, m in
                                     zip(var_stats.values(), timd_mask) if m])

        # Set the analyzed flag to True
        self._analyzed = True

    def print_statistics(self):
        """
        Print the statistics information determined from self analysis
        """

        # Perform self analysis, if needed
        if not self._analyzed:
            self.analyze()

        # Print the statistics information
        for test_name in self._statistics:
            print "Statistics for Test:", test_name
            print

            test_stats = self._statistics[test_name]
            var_stats = test_stats['variables']
            num_steps = test_stats['length']

            print "   Number of Time Steps:", num_steps
            print

            num_tser = test_stats['counts']['tseries']
            print "   Number of Time-Series Variables:            ", num_tser
            num_tvmd = test_stats['counts']['tvariant']
            print "   Number of Time-Variant Metadata Variables:  ", num_tvmd
            num_timd = test_stats['counts']['tinvariant']
            print "   Number of Time-Invariant Metadata Variables:", num_timd
            num_lost = test_stats['counts']['other']
            if num_lost > 0:
                print "   WARNING:", num_lost, " unclassified variables"
            print

            print "   Time-Series Variable Transverse Shapes:"
            print "      ", " ".join(test_stats['xshapes']['tseries'])
            print "   Time-Variant Metadata Transverse Shapes:"
            print "      ", " ".join(test_stats['xshapes']['tvariant'])
            print "   Time-Invariant Metadata Transverse Shapes:"
            print "      ", " ".join(test_stats['xshapes']['tinvariant'])
            print

            tser_totsize = test_stats['totalsizes']['tseries']
            tvmd_totsize = test_stats['totalsizes']['tvariant']
            timd_totsize = test_stats['totalsizes']['tinvariant']

            print "   Time-Series Variable Total Size:   ", _nbyte_str(tser_totsize)
            print "   Time-Variant Metadata Total Size:  ", _nbyte_str(tvmd_totsize)
            print "   Time-Invariant Metadata Total Size:", _nbyte_str(timd_totsize)
            print

            tser_maxsize = test_stats['maxsizes']['tseries']
            tvmd_maxsize = test_stats['maxsizes']['tvariant']
            timd_maxsize = test_stats['maxsizes']['tinvariant']

            print "   Time-Series Variable Max Size:   ", _nbyte_str(tser_maxsize)
            print "   Time-Variant Metadata Max Size:  ", _nbyte_str(tvmd_maxsize)
            print "   Time-Invariant Metadata Max Size:", _nbyte_str(timd_maxsize)
            print

    def save_statistics(self, filename="teststats.json"):
        """
        Save the statistics information to a JSON data file

        Parameters:
            filename (str): The name of the JSON data file to write
        """

        # Check type
        if type(filename) is not str:
            err_msg = "File name must be a string"
            raise TypeError(err_msg)

        # Check if statistics are available
        if self._analyzed:
            json.dump(self._statistics, open(filename, 'w'))

    def load_statistics(self, filename="teststats.json"):
        """
        Load the statistics information from a JSON data file

        Parameters:
            filename (str): The name of the JSON data file to read
        """

        # Check type
        if type(filename) is not str:
            err_msg = "File name must be a string"
            raise TypeError(err_msg)

        # Check if statistics are available
        try:
            statfile = open(filename, 'r')
            self._statistics = dict(json.load(statfile))
            statfile.close()
        except:
            err_msg = "Failed to open and read statistics file '" + \
                str(filename) + "'"
            raise RuntimeError(err_msg)

        # Set the analyzed flag
        self._analyzed = True
