#==============================================================================
#
#  TestTools
#
#  This is a collection of functions that are useful for running the PyReshaper
#  tests on the Yellowstone compute system.
#
#==============================================================================

# Builtin Modules
import glob

# Third-Party Modules
import numpy as np
import Nio

# Package Modules
from pyreshaper import specification


#==============================================================================
# Data Type Mapping & Sizes
#==============================================================================
__DTYPE_MAP = {'d': np.float64, 'f': np.float32, 'l': np.long, 'i': np.int32,
               'h': np.int16, 'b': np.int8, 'S1': np.character}


def __dtype_itemsize(tc):
    dt = np.float
    if tc in __DTYPE_MAP:
        dt = __DTYPE_MAP(tc)
    return np.dtype(dt).itemsize


def __shape_size(shp):
    if (len(shp) > 0):
        return reduce(lambda x, y: x * y, shp)
    else:
        return 1


#==============================================================================
# Byte Size Conversion & Units
#==============================================================================
__BYTE_UNITS = ['Bytes', 'KB', 'MB', 'GB', 'PB']


def __nbyte_str(n, exp=0):
    if (n > 1024.):
        return __nbyte_str(n / 1024., exp=exp + 1)
    else:
        if exp < len(__BYTE_UNITS):
            units = __BYTE_UNITS[exp]
        else:
            n *= 1024.**(exp + 1 - len(__BYTE_UNITS))
            units = __BYTE_UNITS[-1]
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
            abs_path = os.path.abs_path(file_name)
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
        prefix = str(
            os.path.join(output_dir,
                         str(self._database[test_name]['output_prefix'])))
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
            spec = self.create_specifier(test_name, ncfmt='netcdf')

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
            for var_name in infile.variables.keys():
                self._statistics[test_name][var_name] = {}
                var_obj = infile.variables[var_name]

                xshape = var_obj.shape
                tlength = 1
                tvariant = False
                if tdim in var_obj.dimensions:
                    xshape = list(xshape)
                    tindex = var_obj.dimensions.index(tdim)
                    tlength = xshape.pop(tindex)
                    xshape = tuple(xshape)
                    tvariant = True
                self._statistics[test_name][var_name]['tvariant'] = tvariant
                self._statistics[test_name][var_name]['xshape'] = xshape
                self._statistics[test_name][var_name]['length'] = tlength

                xlen = __shape_size(xshape)
                typecode = var_obj.typecode()
                xsize = xlen * __dtype_itemsize(typecode)
                self._statistics[test_name][var_name]['xsize'] = xsize

                if var_name in metadata_names:
                    self._statistics[test_name][var_name]['meta'] = True
                else:
                    self._statistics[test_name][var_name]['meta'] = False

            # Close the first file
            infile.close()

            # Loop over all input files and compute time-variant variable sizes
            for file_name in spec.input_file_list[1:]:

                # Open the file
                infile = Nio.open_file(file_name, 'r')

                # Get the number of time steps in the file
                tlength = infile.dimensions[tdim]

                # Add the time length to each time-variant variable
                for var_name in self._statistics[test_name]:
                    if self._statistics[test_name][var_name]['tvariant']:
                        self._statistics[test_name][
                            var_name]['length'] += tlength

                # Close the file
                infile.close()

        # Set the analyzed flag to Tru
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
            num_vars = len(test_stats)
            meta_mask = [test_stats[v]['meta'] for v in test_stats]
            tvar_mask = [test_stats[v]['tvariant'] for v in test_stats]

            tser_mask = [not meta_mask[i] and tvar_mask[i]
                         for i in range(num_vars)]
            imet_mask = [meta_mask[i] and not tvar_mask[i]
                         for i in range(num_vars)]
            tmet_mask = [meta_mask[i] and tvar_mask[i]
                         for i in range(num_vars)]
            lost_mask = [not meta_mask[i] and not tvar_mask[i]
                         for i in range(num_vars)]

            num_tser = sum(tser_mask)
            print "   Number of Time-Series Variables:            ", num_tser
            num_tmet = sum(tmet_mask)
            print "   Number of Time-Variant Metadata Variables:  ", num_tmet
            num_imet = sum(imet_mask)
            print "   Number of Time-Invariant Metadata Variables:", num_imet
            num_lost = sum(lost_mask)
