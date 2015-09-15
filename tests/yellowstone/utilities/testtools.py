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
import textwrap

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

    def __init__(self, dbname=None, stname=None):
        """
        Initializer

        Parameters:
            dbname (str): The name of the test database file.  Defaults
                to 'testinfo.json'.
            stname (str): The name of the test statistics file.  Defaults
                to 'teststats.json'.

        Raises:
            ValueError: If the test database file cannot be opened and/or
                read.
        """
        # See if there is a user-defined testinfo file,
        # otherwise look for default
        abs_path = ''
        if dbname:
            abs_path = os.path.abspath(dbname)
        else:
            abs_path = os.path.join(os.getcwd(), 'testinfo.json')

        # Try opening and reading the testinfo file
        self._database = {}
        try:
            dbfile = open(abs_path, 'r')
            self._database = dict(json.load(dbfile))
            dbfile.close()
        except:
            err_msg = 'Problem reading and parsing test info file: ' \
                + str(abs_path)
            raise ValueError(err_msg)

        # Initialize the statistics database
        if stname:
            abs_path = os.path.abspath(stname)
        else:
            abs_path = os.path.join(os.getcwd(), 'teststats.json')

        # Try opening and reading the testinfo file
        self._statistics = {}
        try:
            stfile = open(abs_path, 'r')
            self._statistics = dict(json.load(stfile))
            stfile.close()
        except:
            print "Warning: Statistics file (" + stname + ") not found."

    def get_database(self):
        """
        Return the testing database as a dictionary

        Returns:
            dict: The testing database
        """
        return self._database

    def get_statistics(self):
        """
        Return the test analysis statistics as a dictionary

        Returns:
            dict: The statistics database
        """
        return self._statistics

    def print_tests(self):
        """
        List the tests in the test database.
        """
        print 'Tests found in the Test Database are:'
        print
        for test_name in self._database:
            print '   ' + str(test_name)
        return

    def create_specifier(self, test_name, ncfmt='netcdf4c',
                         outdir='', **kwargs):
        """
        Create a Specifier object for the given named test.

        Parameters:
            test_name (str): The string name of the test in the database
                for which to construct the Specifier.
            ncfmt (str): The NetCDF format string to be passed to the
                Specifier.
            outdir (str): An optional path string to be prepended to the
                "output_prefix" argument of the Specifier.  To be used to
                direct output to a different location.  Leave empty if using
                absolute paths in the test's "output_prefix".
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
        if type(outdir) is not str:
            err_msg = "Output directory must be a string"
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
            full_input_glob = str(os.path.join(input_dir, input_glob_str))
            infiles.extend(glob.glob(full_input_glob))
        prefix = str(os.path.join(
            outdir, self._database[test_name]['output_prefix']))
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

    def create_specifiers(self, tests=[], ncfmt='netcdf4c',
                          outdir='', **kwargs):
        """
        Create a dictionary of named Specifier objects for a list of tests.

        Parameters:
            tests (list): A list of string names of tests in the database
                for which to construct Specifiers.
            ncfmt (str): The NetCDF format string to be passed to the
                Specifier.
            outdir (str): An optional path string to be prepended to the
                "output_prefix" argument of the Specifier.  To be used to
                direct output to a different location.  Leave empty if using
                absolute paths in the test's "output_prefix".
            kwargs (dict): A dictionary of additional options to be
                sent to the Specifier.

        Returns:
            list: A list of Specifier instances with the information to run the
                named tests.
        """

        # Check types
        if not isinstance(tests, (list, tuple)):
            err_msg = "Test list must be a list or tuple of string test names"
            raise TypeError(err_msg)
        test_types = [type(test_name) is str for test_name in tests]
        if not all(test_types):
            err_msg = "Test list must be a list of string test names"
            raise TypeError(err_msg)

        # If test list is empty, assume all tests
        if len(tests) == 0:
            test_list = self._database.keys()
        else:
            test_list = tests

        # Construct the list of specifiers
        specs = [self.create_specifier(test_name, ncfmt, outdir, **kwargs)
                 for test_name in test_list]

        # Return a dictionary of named specifiers
        return dict(zip(test_list, specs))

    def analyze(self, tests=[], force=False):
        """
        Analyze the test database to determine test statistics

        Parameters:
            tests (list): A list of string names of tests in the database
                to analyze.  If empty, assume all tests.
            force (bool): Whether to force reanalysis of tests that have
                already been analyzed
        """

        # Check type
        if not isinstance(tests, (list, tuple)):
            err_msg = "Test name list must be of list or tuple type"
            raise TypeError(err_msg)

        # Error if tests not in database
        bad_names = [t for t in tests if t not in self._database]
        if len(bad_names) > 0:
            err_msg = "Tests not found in database: " + ", ".join(bad_names)
            raise ValueError(err_msg)

        # Assume all tests to be analyzed if empty input
        if len(tests) == 0:
            tests = self._database.keys()

        # If analysis has already been done, remove those tests
        if not force:
            for test_name in [t for t in tests if t in self._statistics]:
                print "Not Analyzing Test:", str(test_name)
            tests = [t for t in tests if t not in self._statistics]

        # Generate statistics for each test
        for test_name in tests:
            print "Analyzing Test:", str(test_name)

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

            # Add the transverse (i.e., non-time) dimensions and sizes
            xcoords = dict([(v, s) for (v, s) in infile.dimensions.items()
                            if v != tdim])

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
                    xshape.pop(var_obj.dimensions.index(tdim))
                    xshape = tuple(xshape)
                    tvariant = True
                if not tvariant:
                    metadata_names.add(var_name)
                self._statistics[test_name]['variables'][
                    var_name]['tvariant'] = tvariant
                self._statistics[test_name]['variables'][
                    var_name]['xshape'] = xshape

                xsize = _shape2size(xshape) * _bytesize(var_obj.typecode())
                self._statistics[test_name]['variables'][
                    var_name]['xsize'] = xsize

                if var_name in metadata_names:
                    self._statistics[test_name][
                        'variables'][var_name]['meta'] = True
                else:
                    self._statistics[test_name][
                        'variables'][var_name]['meta'] = False

            # Close the first file
            infile.close()

            # Loop over all input files and compute time-variant variable sizes
            for filename in spec.input_file_list[1:]:

                # Open the file
                infile = Nio.open_file(filename, 'r')

                # And number of time steps to the test data
                self._statistics[test_name]['length'] += \
                    infile.dimensions[tdim]

                # Close the file
                infile.close()

            # Compute self-analysis parameters
            num_steps = self._statistics[test_name]['length']
            var_stats = self._statistics[test_name]['variables']
            tser_vars = [str(v) for (v, s) in var_stats.items()
                         if not s['meta'] and s['tvariant']]
            tvmd_vars = [str(v) for (v, s) in var_stats.items()
                         if s['meta'] and s['tvariant']]
            timd_vars = [str(v) for (v, s) in var_stats.items()
                         if s['meta'] and not s['tvariant']]
            lost_vars = [str(v) for (v, s) in var_stats.items()
                         if not s['meta'] and not s['tvariant']]

            # Store the transverse (to time) coordinate sizes
            self._statistics[test_name]['xcoords'] = xcoords

            # Record the variables names
            self._statistics[test_name]['names'] = {}
            self._statistics[test_name]['names']['tseries'] = tser_vars
            self._statistics[test_name]['names']['tvariant'] = tvmd_vars
            self._statistics[test_name]['names']['tinvariant'] = timd_vars
            self._statistics[test_name]['names']['other'] = lost_vars

            # Compute numbers/counts
            num_tser = len(tser_vars)
            num_tvmd = len(tvmd_vars)
            num_timd = len(timd_vars)
            num_lost = len(lost_vars)
            self._statistics[test_name]['counts'] = {}
            self._statistics[test_name]['counts']['tseries'] = num_tser
            self._statistics[test_name]['counts']['tvariant'] = num_tvmd
            self._statistics[test_name]['counts']['tinvariant'] = num_timd
            self._statistics[test_name]['counts']['other'] = num_lost

            # Compute shapes
            self._statistics[test_name]['xshapes'] = {}
            self._statistics[test_name]['xshapes']['tseries'] = list(
                set([var_stats[v]['xshape'] for v in tser_vars]))
            self._statistics[test_name]['xshapes']['tvariant'] = list(
                set([var_stats[v]['xshape'] for v in tvmd_vars]))
            self._statistics[test_name]['xshapes']['tinvariant'] = list(
                set([var_stats[v]['xshape'] for v in timd_vars]))

            # Compute bytesizes
            self._statistics[test_name]['totalsizes'] = {}
            self._statistics[test_name]['totalsizes']['tseries'] = \
                sum([var_stats[v]['xsize'] for v in tser_vars]) * num_steps
            self._statistics[test_name]['totalsizes']['tvariant'] = \
                sum([var_stats[v]['xsize'] for v in tvmd_vars]) * num_steps
            self._statistics[test_name]['totalsizes']['tinvariant'] = \
                sum([var_stats[v]['xsize'] for v in timd_vars])

            # Compute maxima
            self._statistics[test_name]['maxsizes'] = {}
            maxsize = 0 if num_tser == 0 else \
                max([var_stats[v]['xsize'] for v in tser_vars]) * num_steps
            self._statistics[test_name]['maxsizes']['tseries'] = maxsize
            maxsize = 0 if num_tvmd == 0 else \
                max([var_stats[v]['xsize'] for v in tvmd_vars]) * num_steps
            self._statistics[test_name]['maxsizes']['tvariant'] = maxsize
            maxsize = 0 if num_timd == 0 else \
                max([var_stats[v]['xsize'] for v in timd_vars])
            self._statistics[test_name]['maxsizes']['tinvariant'] = maxsize

    def print_statistics(self, tests=[]):
        """
        Print the statistics information determined from self analysis

        Parameters:
            tests (list): A list of string names of tests in the database
                to print.  If empty, assume all tests.
        """

        # Check type
        if not isinstance(tests, (list, tuple)):
            err_msg = "Test name list must be of list or tuple type"
            raise TypeError(err_msg)

        # Error if tests not in database
        bad_names = [t for t in tests if t not in self._statistics]
        if len(bad_names) > 0:
            err_msg = "Tests not found in statistics: " + ", ".join(bad_names)
            raise ValueError(err_msg)

        # Assume all tests to be analyzed if empty input
        if len(tests) == 0:
            tests = self._statistics.keys()

        # Print the statistics information
        for test_name in self._statistics:
            print "Statistics for Test:", test_name
            print

            test_stats = self._statistics[test_name]
            num_steps = test_stats['length']

            print "   Number of Time Steps:", num_steps
            print

            # Print counts
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

            # Print the coordinate data
            print "   Transverse Coordinate Shapes:"
            maxlenxcoord = max([len(xc) for xc in test_stats['xcoords']])
            for xcoord, cxsize in test_stats['xcoords'].items():
                spcr = ' ' * (maxlenxcoord - len(xcoord))
                print "      " + xcoord + ": " + spcr + str(cxsize)
            print

            # Print names
            print "   Time-Series Variables:"
            vlist = ", ".join([str(v) for v in
                               test_stats['names']['tseries']])
            print "      " + "\n      ".join(textwrap.wrap(vlist))
            print "   Time-Variant Metadata Variables:"
            vlist = ", ".join([str(v) for v in
                               test_stats['names']['tvariant']])
            print "      " + "\n      ".join(textwrap.wrap(vlist))
            print "   Time-Invariant Metadata Variables:"
            vlist = ", ".join([str(v) for v in
                               test_stats['names']['tinvariant']])
            print "      " + "\n      ".join(textwrap.wrap(vlist))
            if num_lost > 0:
                print "   Unclassified Variables (neither meta nor time-variant):"
                vlist = ", ".join([str(v) for v in
                                   test_stats['names']['other']])
                print "      " + "\n      ".join(textwrap.wrap(vlist))
            print

            # Print Transverse Shapes
            print "   Time-Series Variable Transverse Shapes:"
            print "      " + " ".join([str(s) for s in
                                       test_stats['xshapes']['tseries']])
            print "   Time-Variant Metadata Transverse Shapes:"
            print "      " + " ".join([str(s) for s in
                                       test_stats['xshapes']['tvariant']])
            print "   Time-Invariant Metadata Transverse Shapes:"
            print "      " + " ".join([str(s) for s in
                                       test_stats['xshapes']['tinvariant']])
            print

            # Print total bytesizes
            tser_totsize = test_stats['totalsizes']['tseries']
            print "   Time-Series Variable Total Size:   ", _nbyte_str(tser_totsize)
            tvmd_totsize = test_stats['totalsizes']['tvariant']
            print "   Time-Variant Metadata Total Size:  ", _nbyte_str(tvmd_totsize)
            timd_totsize = test_stats['totalsizes']['tinvariant']
            print "   Time-Invariant Metadata Total Size:", _nbyte_str(timd_totsize)
            print

            # Print maximum bytesizes
            tser_maxsize = test_stats['maxsizes']['tseries']
            print "   Time-Series Variable Max Size:   ", _nbyte_str(tser_maxsize)
            tvmd_maxsize = test_stats['maxsizes']['tvariant']
            print "   Time-Variant Metadata Max Size:  ", _nbyte_str(tvmd_maxsize)
            timd_maxsize = test_stats['maxsizes']['tinvariant']
            print "   Time-Invariant Metadata Max Size:", _nbyte_str(timd_maxsize)
            print

    def save_statistics(self, filename="teststats.json"):
        """
        Save the statistics information to a JSON data file

        Parameters:
            filename (str): The name of the JSON data file to write
        """

        # Check types
        if isinstance(filename, str):
            fp = open(filename, 'w')
        else:
            err_msg = "Statistics filename must be a string"
            raise TypeError(err_msg)

        # Dump JSON data to file
        try:
            json.dump(self._statistics, fp)
        except:
            err_msg = "Failed to write statistics file"
            raise RuntimeError(err_msg)

        # Close the file
        fp.close()

    def load_statistics(self, filename="teststats.json"):
        """
        Load the statistics information from a JSON data file

        Parameters:
            filename (str): The name of the JSON data file to read
        """

        # Check types
        if isinstance(filename, str):
            fp = open(filename, 'r')
        else:
            err_msg = "Statistics filename must be a string"
            raise TypeError(err_msg)

        # Try reading the statistics
        try:
            self._statistics = dict(json.load(fp))
        except:
            err_msg = "Failed to read statistics file"
            raise RuntimeError(err_msg)

        # Close the file
        fp.close()
