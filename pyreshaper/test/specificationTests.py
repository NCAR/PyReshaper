'''
Specifier unit tests

------------------------
Created on Apr 30, 2014

@author: kpaul
'''

import unittest
from pyreshaper import specification


class SpecifierTests(unittest.TestCase):
    '''
    SpecifierTests Class

    This class defines all of the unit tests for the specification module.
    '''

    def test_init(self):
        spec = specification.Slice2SeriesSpecifier()
        self.assertListEqual(spec.input_file_list, [],
                             'Input file list not initialized to empty')
        self.assertEqual(spec.netcdf_format, 'netcdf4c',
                         'NetCDF format not initialized to netcdf4c')
        self.assertEqual(spec.output_file_prefix, 'tseries.',
                         'Output file prefix not initialized to tseries.')
        self.assertEqual(spec.output_file_suffix, '.nc',
                         'Output file prefix not initialized to .nc')
        self.assertListEqual(spec.time_variant_metadata, [],
            'Time variant metadata list not initialized to empty')

    def test_init_full(self):
        in_list = ['a', 'b', 'c']
        fmt = 'netcdf'
        prefix = 'pre.'
        suffix = '.suf.nc'
        metadata = ['x', 'y', 'z']
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        self.assertListEqual(spec.input_file_list, in_list,
                             'Input file list not initialized properly')
        self.assertEqual(spec.netcdf_format, fmt,
                         'NetCDF format not initialized properly')
        self.assertEqual(spec.output_file_prefix, prefix,
                         'Output file prefix not initialized properly')
        self.assertEqual(spec.output_file_suffix, suffix,
                         'Output file prefix not initialized properly')
        self.assertListEqual(spec.time_variant_metadata, metadata,
            'Time variant metadata list not initialized properly')

    def test_validate_types(self):
        in_list = ['a', 'b', 'c']
        fmt = 'netcdf'
        prefix = 'pre.'
        suffix = '.suf.nc'
        metadata = ['x', 'y', 'z']
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        spec.validate_types()

    def test_validate_types_fail_1(self):
        in_list = ['a', 2, 'c']
        fmt = 'netcdf'
        prefix = 'pre.'
        suffix = '.suf.nc'
        metadata = ['x', 'y', 'z']
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        self.assertRaises(TypeError, spec.validate_types)

    def test_validate_types_fail_2(self):
        in_list = ['a', 'b', 'c']
        fmt = 2342
        prefix = 'pre.'
        suffix = '.suf.nc'
        metadata = ['x', 'y', 'z']
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        self.assertRaises(TypeError, spec.validate_types)

    def test_validate_types_fail_3(self):
        in_list = ['a', 'b', 'c']
        fmt = 'netcdf'
        prefix = dict()
        suffix = '.suf.nc'
        metadata = ['x', 'y', 'z']
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        self.assertRaises(TypeError, spec.validate_types)

    def test_validate_types_fail_4(self):
        in_list = ['a', 'b', 'c']
        fmt = 'netcdf'
        prefix = 'pre.'
        suffix = list()
        metadata = ['x', 'y', 'z']
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        self.assertRaises(TypeError, spec.validate_types)

    def test_validate_types_fail_5(self):
        in_list = ['a', 'b', 'c']
        fmt = 'netcdf'
        prefix = 'pre.'
        suffix = '.suf.nc'
        metadata = ['x', 'y', 2]
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        self.assertRaises(TypeError, spec.validate_types)

    def test_validate_values_fail_1(self):
        in_list = ['a', 'b', 'c']
        fmt = 'netcdf'
        prefix = 'pre.'
        suffix = '.suf.nc'
        metadata = []
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        spec.validate_types()
        self.assertRaises(ValueError, spec.validate_values)

    def test_validate_values_fail_2(self):
        in_list = ['timekeeperTests.py', 'messengerTests.py']
        fmt = 'netcdf9'
        prefix = 'pre.'
        suffix = '.suf.nc'
        metadata = []
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        spec.validate_types()
        self.assertRaises(ValueError, spec.validate_values)

    def test_validate_values_fail_3(self):
        in_list = ['timekeeperTests.py', 'messengerTests.py']
        fmt = 'netcdf4'
        prefix = '/sfcsrytsdfv/pre.'
        suffix = '.suf.nc'
        metadata = []
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        spec.validate_types()
        self.assertRaises(ValueError, spec.validate_values)

    def test_validate_values_mod(self):
        in_list = ['timekeeperTests.py', 'messengerTests.py']
        fmt = 'netcdf4'
        prefix = 'pre.'
        suffix = '.suf'
        metadata = []
        spec = specification.Slice2SeriesSpecifier(infiles=in_list,
            ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata)
        spec.validate_types()
        spec.validate_values()
        self.assertEqual(spec.output_file_suffix, suffix + '.nc',
                         'Suffix was not changed to .nc extension')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
