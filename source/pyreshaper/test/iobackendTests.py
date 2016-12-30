"""
Unit tests for the iobackend module

Copyright 2016, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

from __future__ import print_function, absolute_import

import unittest
import numpy as np
from numpy import testing as npt
from functools import reduce

from pyreshaper import iobackend
from os import linesep, remove
from os.path import exists

try:
    import netCDF4
except ImportError:
    netCDF4 = None

try:
    import Nio
except ImportError:
    Nio = None


#===============================================================================
# print_test_msg
#===============================================================================
def print_test_msg(testname, indata=None, actual=None, expected=None):
    msg = '{0}:{1}'.format(testname, linesep)
    if indata is not None:
        msg += '   - input:    {0!r}{1}'.format(indata, linesep)
    if actual is not None:
        msg += '   - actual:   {0!r}{1}'.format(actual, linesep)
    if expected is not None:
        msg += '   - expected: {0!r}{1}'.format(expected, linesep)
    print(msg)


#===============================================================================
# ReadTests
#===============================================================================
class ReadTests(unittest.TestCase):
    """
    Read Tests TestCase Class
    """
    
    def setUp(self):
        self.fname = 'readtest.nc'
        self.cleanUp()
        
        self.fatts = {'a1': 'attribute 1', 'a2': 'attribute 2'}
        self.uldim = 't'
        self.fdims = {'t': 10, 'x': 5, 'c': 16}
        self.vtyps = {'T': 'd', 'X': 'd', 'V': 'f'}
        self.vdims = {'T': ('t',), 'X': ('x',), 'V': ('t', 'x')}
        self.vdata = {'T': np.arange(0, self.fdims['t'], dtype=self.vtyps['T']),
                      'X': np.random.ranf(self.fdims['x']).astype(self.vtyps['X']),
                      'V': np.random.ranf(self.fdims['t']*self.fdims['x']).reshape(
                          self.fdims['t'], self.fdims['x']).astype(self.vtyps['V'])}
        self.vatts = {'T': {'long_name': 'time', 'units': 'sec'},
                      'X': {'long_name': 'space', 'units': 'meters'},
                      'V': {'long_name': 'variable', 'units': 'vunits'}}
        self.backend = None
    
    def cleanUp(self):
        if exists(self.fname):
            remove(self.fname)
        
    def tearDown(self):
        self.cleanUp()
        
    def test_avail(self):
        actual = iobackend._AVAILABLE_
        print_test_msg('{0}: _AVAILABLE_'.format(self.backend), actual=actual)
        self.assertTrue(self.backend in iobackend._AVAILABLE_,
                        '{0} not available'.format(self.backend))

    def test_set_backend(self):
        indata = self.backend
        iobackend.set_backend(indata)
        actual = iobackend._BACKEND_
        expected = indata
        print_test_msg('{0}: set_backend({0})'.format(indata), indata, actual, expected)
        self.assertEqual(iobackend._BACKEND_, indata,
                        '{0} backend name not set'.format(self.backend))

    def test_set_backend_x(self):
        indata = 'x'
        actual = iobackend._BACKEND_
        print_test_msg('{0}: set_backend({1})'.format(self.backend, indata), indata, actual)
        self.assertRaises(KeyError, iobackend.set_backend, indata)

    def test_NCFile_init_mode_x(self):
        expected = ValueError
        print_test_msg('{0}: NCFile.__init__(mode=x)'.format(self.backend), expected=expected)
        self.assertRaises(expected, iobackend.NCFile, self.fname, mode='x')

    def test_NCFile_init_read(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        actual = type(ncf)
        expected = iobackend.NCFile
        ncf.close()
        print_test_msg('{0}: NCFile.__init__()'.format(self.backend),
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected, 'NCFile not created with correct type')
        
    def test_NCFile_dimensions(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        actual = ncf.dimensions
        expected = self.fdims
        ncf.close()
        print_test_msg('{0}: NCFile.dimensions'.format(self.backend),
                       actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected), 'NCFile dimensions not correct length')
        for dn in expected:
            self.assertTrue(dn in actual, 'NCFile dimension {0!r} not present'.format(dn))
            self.assertEqual(actual[dn], expected[dn], 'NCFile dimension {0!r} not correct'.format(dn))
            
    def test_NCFile_unlimited(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        for dn in ncf.dimensions:
            actual = ncf.unlimited(dn)
            expected = (dn == self.uldim)
            print_test_msg('{0}: NCFile.unlimited({1})'.format(self.backend, dn),
                           actual=actual, expected=expected)
            self.assertEqual(actual, expected, 'NCFile dimension has incorrect unlimited status')
        ncf.close()

    def test_NCFile_ncattrs(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        actual = ncf.ncattrs
        expected = list(self.fatts)
        print_test_msg('{0}: NCFile.ncattrs'.format(self.backend),
                       actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected), 'NCFile ncattrs not correct length')
        for dn in expected:
            self.assertTrue(dn in actual,
                            'NCFile ncattrs {0!r} not present'.format(dn))
            self.assertEqual(ncf.getncattr(dn), self.fatts[dn],
                            'NCFile ncattrs {0!r} not correct'.format(dn))
        ncf.close()

    def test_NCFile_variables(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        actual = ncf.variables
        print_test_msg('Nio: NCFile.variables', actual=list(actual))
        for vn in self.vtyps:
            self.assertTrue(vn in actual, '{0} variable not in NCFile'.format(vn))
        for vn in actual:
            self.assertTrue(isinstance(actual[vn], iobackend.NCVariable),
                            'Variable {0!r} has wrong type'.format(vn))
        ncf.close()
            
    def test_NCFile_close(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        actual = ncf.close()
        print_test_msg('{0}: NCFile.close'.format(self.backend), actual=actual)

    def test_NCVariable_ncattrs(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        for vn in ncf.variables:
            actual = ncf.variables[vn].ncattrs
            expected = list(self.vatts[vn])
            print_test_msg('{0}: NCVariable {1} ncattrs'.format(self.backend, vn),
                           actual=actual, expected=expected)
            self.assertEqual(len(actual), len(expected),
                             'NCVariable {0} ncattrs not correct length'.format(vn))
            for a in expected:
                self.assertTrue(a in actual,
                                'Attribute {0!r} not found in variable {1}'.format(a, vn))
                self.assertEqual(ncf.variables[vn].getncattr(a), self.vatts[vn][a],
                                 'Attribute {0!r} not correct in variable {1}'.format(a, vn))
        ncf.close()

    def test_NCVariable_dimensions(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        for vn in ncf.variables:
            actual = ncf.variables[vn].dimensions
            expected = self.vdims[vn]
            print_test_msg('{0}: NCVariable {1} dimensions'.format(self.backend, vn),
                           actual=actual, expected=expected)
            self.assertEqual(actual, expected, 'NCVariable {0} dimensions not correct'.format(vn))
        ncf.close()

    def test_NCVariable_datatype(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        for vn in ncf.variables:
            actual = ncf.variables[vn].datatype
            expected = np.dtype(self.vtyps[vn])
            print_test_msg('{0}: NCVariable {1} datatype'.format(self.backend, vn),
                           actual=actual, expected=expected)
            self.assertEqual(actual, expected, 'NCVariable {0} datatype not correct'.format(vn))
        ncf.close()
        
    def test_NCVariable_shape(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        for vn in ncf.variables:
            actual = ncf.variables[vn].shape
            expected = tuple(self.fdims[dn] for dn in self.vdims[vn])
            print_test_msg('{0}: NCVariable {1} shape'.format(self.backend, vn),
                           actual=actual, expected=expected)
            self.assertEqual(actual, expected, 'NCVariable {0} shape not correct'.format(vn))
        ncf.close()
        
    def test_NCVariable_size(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        for vn in ncf.variables:
            actual = ncf.variables[vn].size
            expected = reduce(lambda x,y: x*y, (self.fdims[dn] for dn in self.vdims[vn]))
            print_test_msg('{0}: NCVariable {1} size'.format(self.backend, vn),
                           actual=actual, expected=expected)
            self.assertEqual(actual, expected, 'NCVariable {0} size not correct'.format(vn))
        ncf.close()

    def test_NCVariable_getitem(self):
        iobackend.set_backend(self.backend)
        ncf = iobackend.NCFile(self.fname)
        for vn in ncf.variables:
            actual = ncf.variables[vn][:]
            expected = self.vdata[vn][:]
            print_test_msg('{0}: NCVariable[{1}].__getitem__'.format(self.backend, vn),
                           actual=actual, expected=expected)
            npt.assert_array_equal(actual, expected)
        ncf.close()


#===================================================================================================
# NioReadTests
#===================================================================================================
class NioReadTests(ReadTests):
    
    def setUp(self):
        super(NioReadTests, self).setUp()
        self.backend = 'Nio'

        ncfile = Nio.open_file(self.fname, 'c')
        for k in self.fatts:
            setattr(ncfile, k, self.fatts[k])
        for d in self.fdims:
            if d == self.uldim:
                ncfile.create_dimension(d, None)
            else:
                ncfile.create_dimension(d, self.fdims[d])
        for n in self.vtyps:
            ncvar = ncfile.create_variable(n, self.vtyps[n], self.vdims[n])
            for k in self.vatts[n]:
                setattr(ncvar, k, self.vatts[n][k])
            ncvar.assign_value(self.vdata[n])
        ncfile.close()


#===================================================================================================
# NC4ReadTests
#===================================================================================================
class NC4ReadTests(ReadTests):
    
    def setUp(self):
        super(NC4ReadTests, self).setUp()
        self.backend = 'netCDF4'

        ncfile = netCDF4.Dataset(self.fname, 'w')
        for k in self.fatts:
            ncfile.setncattr(k, self.fatts[k])
        for d in self.fdims:
            if d == self.uldim:
                ncfile.createDimension(d)
            else:
                ncfile.createDimension(d, self.fdims[d])
        for n in self.vtyps:
            ncvar = ncfile.createVariable(n, self.vtyps[n], self.vdims[n])
            for k in self.vatts[n]:
                ncvar.setncattr(k, self.vatts[n][k])
            ncvar[:] = self.vdata[n]
        ncfile.close()    


# #===============================================================================
# # IOBackendWriteTests
# #===============================================================================
# class IOBackendWriteTests(unittest.TestCase):
# 
#     """
#     IOBackendWriteTests Class
# 
#     This class defines all of the unit tests for the iobackend module.
#     """
#     
#     def setUp(self):
#         self.ncfwname = 'writetest.nc'
#         self.ncattrs = {'a1': 'attribute 1',
#                         'a2': 'attribute 2'}
#         self.ncdims = {'t': 10, 'x': 5}
#         self.t = np.arange(0, self.ncdims['t'], dtype='d')
#         self.x = np.random.ranf(self.ncdims['x']).astype('d')
#         self.v = np.random.ranf(self.ncdims['t']*self.ncdims['x']).reshape(10,5).astype('f')
#         self.vattrs = {'long_name': 'variable',
#                        'units': 'meters'}
#     
#     def tearDown(self):
#         if exists(self.ncfwname):
#             remove(self.ncfwname)
# 
#     def test_nio_NCFile_init_write(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         actual = type(ncf)
#         ncf.close()
#         expected = iobackend.NCFile
#         print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile not created with correct type')
# 
#     def test_nc4_NCFile_init_write(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         actual = type(ncf)
#         ncf.close()
#         expected = iobackend.NCFile
#         print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile not created with correct type')
# 
#     def test_nio_NCFile_setncattr(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         for a,v in self.ncattrs.iteritems():
#             ncf.setncattr(a, v)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.attributes
#         expected = self.ncattrs
#         ncfr.close()
#         print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
#         for a,v in expected.iteritems():
#             self.assertTrue(a in actual,
#                             'NCFile attribute {0!r} not found'.format(a))
#             self.assertEqual(actual[a], v,
#                              'NCFile attribute {0!r} incorrect'.format(a))
# 
#     def test_nc4_NCFile_setncattr(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         for a,v in self.ncattrs.iteritems():
#             ncf.setncattr(a, v)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.attributes
#         expected = self.ncattrs
#         ncfr.close()
#         print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
#         for a,v in expected.iteritems():
#             self.assertTrue(a in actual,
#                             'NCFile attribute {0!r} not found'.format(a))
#             self.assertEqual(actual[a], v,
#                              'NCFile attribute {0!r} incorrect'.format(a))
# 
#     def test_nio_NCFile_create_dimension(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('x', self.ncdims['x'])
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.dimensions['x']
#         expected = self.ncdims['x']
#         ncfr.close()
#         print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile x-dimension incorrect')
# 
#     def test_nc4_NCFile_create_dimension(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('x', self.ncdims['x'])
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.dimensions['x']
#         expected = self.ncdims['x']
#         ncfr.close()
#         print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile x-dimension incorrect')
# 
#     def test_nio_NCFile_create_dimension_unlimited(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.dimensions['t']
#         expected = 0
#         ncfr.close()
#         print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile t-dimension incorrect')
# 
#     def test_nc4_NCFile_create_dimension_unlimited(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.dimensions['t']
#         expected = 0
#         ncfr.close()
#         print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile t-dimension incorrect')
# 
#     def test_nio_NCFile_create_variable(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('x', self.ncdims['x'])
#         x = ncf.create_variable('x', np.dtype('d'), ('x',))
#         x[:] = self.x
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['x'][:]
#         expected = self.x
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile x-variable incorrect')
# 
#     def test_nc4_NCFile_create_variable(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('x', self.ncdims['x'])
#         x = ncf.create_variable('x', np.dtype('d'), ('x',))
#         x[:] = self.x
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['x'][:]
#         expected = self.x
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile x-variable incorrect')
# 
#     def test_nio_NCFile_create_variable_unlimited(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         t = ncf.create_variable('t', np.dtype('d'), ('t',))
#         t[:] = self.t
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['t'][:]
#         expected = self.t
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile t-variable incorrect')
# 
#     def test_nc4_NCFile_create_variable_unlimited(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         t = ncf.create_variable('t', np.dtype('d'), ('t',))
#         t[:] = self.t
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['t'][:]
#         expected = self.t
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile t-variable incorrect')
# 
#     def test_nio_NCFile_create_variable_ndim(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         ncf.create_dimension('x', self.ncdims['x'])
#         v = ncf.create_variable('v', np.dtype('f'), ('t', 'x'))
#         v[:] = self.v
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['v'][:]
#         expected = self.v
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile 2d-variable incorrect')
# 
#     def test_nc4_NCFile_create_variable_ndim(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         ncf.create_dimension('x', self.ncdims['x'])
#         v = ncf.create_variable('v', np.dtype('f'), ('t', 'x'))
#         v[:] = self.v
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['v'][:]
#         expected = self.v
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile 2d-variable incorrect')
# 
#     def test_nio_NCVariable_setncattr(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         ncf.create_dimension('x', self.ncdims['x'])
#         v = ncf.create_variable('v', np.dtype('f'), ('t', 'x'))
#         for attr,value in self.vattrs.iteritems():
#             v.setncattr(attr, value)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['v'].attributes
#         expected = self.vattrs
#         ncfr.close()
#         print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
#         for attr, value in expected.iteritems():
#             self.assertTrue(attr in actual,
#                             'Variable attribute {0!r} not found'.format(attr))
#             self.assertEqual(actual[attr], value,
#                              'Variable attribute {0!r} incorrect'.format(attr))
# 
#     def test_nc4_NCVariable_setncattr(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfwname, mode='w')
#         ncf.create_dimension('t')
#         ncf.create_dimension('x', self.ncdims['x'])
#         v = ncf.create_variable('v', np.dtype('f'), ('t', 'x'))
#         for attr,value in self.vattrs.iteritems():
#             v.setncattr(attr, value)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfwname)
#         actual = ncfr.variables['v'].attributes
#         expected = self.vattrs
#         ncfr.close()
#         print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
#         for attr, value in expected.iteritems():
#             self.assertTrue(attr in actual,
#                             'Variable attribute {0!r} not found'.format(attr))
#             self.assertEqual(actual[attr], value,
#                              'Variable attribute {0!r} incorrect'.format(attr))
# 
# 
# #===============================================================================
# # IOBackendAppendTests
# #===============================================================================
# class IOBackendAppendTests(unittest.TestCase):
# 
#     """
#     IOBackendAppendTests Class
# 
#     This class defines all of the unit tests for the iobackend module.
#     """
#     
#     def setUp(self):
#         self.ncfaname = 'appendtest.nc'
#         self.ncattrs = {'a1': 'attribute 1',
#                         'a2': 'attribute 2'}
#         self.ncdims = {'t': 10, 'x': 5}
#         self.t = np.arange(self.ncdims['t'], dtype='d')
#         self.x = np.arange(self.ncdims['x'], dtype='d')
#         self.v = np.arange(self.ncdims['t']*self.ncdims['x'],
#                            dtype='f').reshape(self.ncdims['t'], self.ncdims['x'])
#         self.vattrs = {'long_name': 'variable',
#                        'units': 'meters'}
#         
#         self.fattrs2 = {'a3': 'attribute 3',
#                         'a4': 'attribute 4'}
#         self.t2 = np.arange(self.ncdims['t'], 2*self.ncdims['t'], dtype='d')
#         self.v2 = np.arange(self.ncdims['t']*self.ncdims['x'],
#                             dtype='f').reshape(self.ncdims['t'], self.ncdims['x'])
#         self.vattrs2 = {'standard_name': 'variable'}
# 
#         ncfile = netCDF4.Dataset(self.ncfaname, 'w')
#         for a,v in self.ncattrs.iteritems():
#             setattr(ncfile, a, v)
#         ncfile.createDimension('t')
#         ncfile.createDimension('x', self.ncdims['x'])
#         t = ncfile.createVariable('t', 'd', ('t',))
#         t[:] = self.t
#         x = ncfile.createVariable('x', 'd', ('x',))
#         x[:] = self.x
#         v = ncfile.createVariable('v', 'f', ('t', 'x'))
#         for a,val in self.vattrs.iteritems():
#             v.setncattr(a, val)
#         v[:,:] = self.v
#         
#         ncfile.close()
#     
#     def tearDown(self):
#         if exists(self.ncfaname):
#             remove(self.ncfaname)
# 
#     def test_nio_NCFile_init_append(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         actual = type(ncf)
#         ncf.close()
#         expected = iobackend.NCFile
#         print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile not created with correct type')
# 
#     def test_nc4_NCFile_init_append(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         actual = type(ncf)
#         ncf.close()
#         expected = iobackend.NCFile
#         print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
#         self.assertEqual(actual, expected,
#                          'NCFile not created with correct type')
# 
#     def test_nio_NCFile_setncattr(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         for a,v in self.fattrs2.iteritems():
#             ncf.setncattr(a, v)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.attributes
#         expected = self.ncattrs
#         expected.update(self.fattrs2)
#         ncfr.close()
#         print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
#         for a,v in expected.iteritems():
#             self.assertTrue(a in actual,
#                             'NCFile attribute {0!r} not found'.format(a))
#             self.assertEqual(actual[a], v,
#                              'NCFile attribute {0!r} incorrect'.format(a))
# 
#     def test_nc4_NCFile_setncattr(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         for a,v in self.fattrs2.iteritems():
#             ncf.setncattr(a, v)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.attributes
#         expected = self.ncattrs
#         expected.update(self.fattrs2)
#         ncfr.close()
#         print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
#         for a,v in expected.iteritems():
#             self.assertTrue(a in actual,
#                             'NCFile attribute {0!r} not found'.format(a))
#             self.assertEqual(actual[a], v,
#                              'NCFile attribute {0!r} incorrect'.format(a))
# 
#     def test_nio_NCFile_create_variable_ndim(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         v2 = ncf.create_variable('v2', np.dtype('f'), ('t', 'x'))
#         v2[:] = self.v2
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.variables['v2'][:]
#         expected = self.v2
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile 2d-variable incorrect')
# 
#     def test_nc4_NCFile_create_variable_ndim(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         v2 = ncf.create_variable('v2', np.dtype('f'), ('t', 'x'))
#         v2[:] = self.v2
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.variables['v2'][:]
#         expected = self.v2
#         ncfr.close()
#         print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile 2d-variable incorrect')
# 
#     def test_nio_NCFile_variable_append(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         nt = self.ncdims['t']
#         t = ncf.variables['t']
#         t[nt:] = self.t2
#         v = ncf.variables['v']
#         v[nt:, :] = self.v2
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.variables['t'][:]
#         expected = np.concatenate((self.t, self.t2))
#         print_test_msg('NCVariable append', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile t-variable incorrect')
#         actual = ncfr.variables['v'][:]
#         expected = np.concatenate((self.v, self.v2))
#         print_test_msg('NCVariable append', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile 2d-variable incorrect')
#         ncfr.close()
# 
#     def test_nc4_NCFile_variable_append(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         nt = self.ncdims['t']
#         t = ncf.variables['t']
#         t[nt:] = self.t2
#         v = ncf.variables['v']
#         v[nt:, :] = self.v2
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.variables['t'][:]
#         expected = np.concatenate((self.t, self.t2))
#         print_test_msg('NCVariable append', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile t-variable incorrect')
#         actual = ncfr.variables['v'][:]
#         expected = np.concatenate((self.v, self.v2))
#         print_test_msg('NCVariable append', actual=actual, expected=expected)
#         npt.assert_array_equal(actual, expected,
#                                'NCFile 2d-variable incorrect')
#         ncfr.close()
# 
#     def test_nio_NCVariable_setncattr(self):
#         iobackend.set_backend('Nio')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         v = ncf.variables['v']
#         for attr,value in self.vattrs2.iteritems():
#             v.setncattr(attr, value)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.variables['v'].attributes
#         expected = self.vattrs
#         expected.update(self.vattrs2)
#         ncfr.close()
#         print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
#         for attr, value in expected.iteritems():
#             self.assertTrue(attr in actual,
#                             'Variable attribute {0!r} not found'.format(attr))
#             self.assertEqual(actual[attr], value,
#                              'Variable attribute {0!r} incorrect'.format(attr))
# 
#     def test_nc4_NCVariable_setncattr(self):
#         iobackend.set_backend('netCDF4')
#         ncf = iobackend.NCFile(self.ncfaname, mode='a')
#         v = ncf.variables['v']
#         for attr,value in self.vattrs2.iteritems():
#             v.setncattr(attr, value)
#         ncf.close()
#         ncfr = Nio.open_file(self.ncfaname)
#         actual = ncfr.variables['v'].attributes
#         expected = self.vattrs
#         expected.update(self.vattrs2)
#         ncfr.close()
#         print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
#         for attr, value in expected.iteritems():
#             self.assertTrue(attr in actual,
#                             'Variable attribute {0!r} not found'.format(attr))
#             self.assertEqual(actual[attr], value,
#                              'Variable attribute {0!r} incorrect'.format(attr))


#===============================================================================
# CLI
#===============================================================================
if __name__ == "__main__":
    import sys
    sys.argv = ['']
    if Nio is not None:
        print('RUNNING PYNIO TESTS')
        sys.argv.append('NioReadTests')
    if netCDF4 is not None:
        print('RUNNING NETCDF4-PYTHON TESTS')
        sys.argv.append('NC4ReadTests')
    unittest.main()
