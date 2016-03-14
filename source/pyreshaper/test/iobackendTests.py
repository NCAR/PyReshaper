"""
Unit tests for the iobackend module

Copyright 2016, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

import unittest
import numpy as np
import numpy.testing as npt
import netCDF4
import Nio

from pyreshaper import iobackend
from os import linesep, remove
from os.path import exists


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
    print msg
    

#===============================================================================
# IOBackendReadTests
#===============================================================================
class IOBackendReadTests(unittest.TestCase):

    """
    IOBackendReadTests Class

    This class defines all of the unit tests for the iobackend module.
    """
    
    def setUp(self):
        self.ncfrname = 'readtest.nc'
        self.ncattrs = {'a1': 'attribute 1',
                        'a2': 'attribute 2'}
        self.ncdims = {'t': 10, 'x': 5}
        self.t = np.arange(0, self.ncdims['t'], dtype='d')
        self.x = np.random.ranf(self.ncdims['x']).astype('d')
        self.v = np.random.ranf(self.ncdims['t']*self.ncdims['x']).reshape(10,5).astype('f')
        self.vattrs = {'long_name': 'variable',
                       'units': 'meters'}

        ncfile = netCDF4.Dataset(self.ncfrname, 'w')
        for a,v in self.ncattrs.iteritems():
            setattr(ncfile, a, v)
        ncfile.createDimension('t')
        ncfile.createDimension('x', self.ncdims['x'])
        t = ncfile.createVariable('t', 'd', ('t',))
        t[:] = self.t
        x = ncfile.createVariable('x', 'd', ('x',))
        x[:] = self.x
        v = ncfile.createVariable('v', 'f', ('t', 'x'))
        for a,val in self.vattrs.iteritems():
            v.setncattr(a, val)
        v[:,:] = self.v
        
        ncfile.close()
    
    def tearDown(self):
        if exists(self.ncfrname):
            remove(self.ncfrname)

    def test_avail(self):
        actual = iobackend._AVAIL_
        print_test_msg('_AVAIL_', actual=actual)
        self.assertTrue('Nio' in iobackend._AVAIL_,
                        'Nio importable but not available')
        self.assertTrue('netCDF4' in iobackend._AVAIL_,
                        'netCDF4 importable but not available')

    def test_set_backend_nio(self):
        indata = 'Nio'
        iobackend.set_backend(indata)
        actual = iobackend._BACKEND_
        expected = indata
        print_test_msg('set_backend()', indata, actual, expected)
        self.assertEqual(iobackend._BACKEND_, indata,
                        'PyNIO backend name not set')

    def test_set_backend_nc4(self):
        indata = 'netCDF4'
        iobackend.set_backend(indata)
        actual = iobackend._BACKEND_
        expected = indata
        print_test_msg('set_backend()', indata, actual, expected)
        self.assertEqual(iobackend._BACKEND_, indata,
                        'netCDF4 backend name not set')

    def test_set_backend_x(self):
        indata = 'x'
        actual = iobackend._BACKEND_
        print_test_msg('set_backend()', indata, actual)
        self.assertRaises(KeyError, iobackend.set_backend, indata)

    def test_NCFile_init_mode_x(self):
        expected = ValueError
        print_test_msg('NCFile.__init__(mode=x)', expected=expected)
        self.assertRaises(expected, iobackend.NCFile, self.ncfrname,
                          mode='x')

    def test_nio_NCFile_init_read(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = type(ncf)
        expected = iobackend.NCFile
        ncf.close()
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nc4_NCFile_init_read(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = type(ncf)
        expected = iobackend.NCFile
        ncf.close()
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nio_NCFile_dimensions(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.dimensions
        expected = self.ncdims
        ncf.close()
        print_test_msg('NCFile.dimensions', actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected),
                         'NCFile dimensions not correct length')
        for dn, dv in expected.iteritems():
            self.assertTrue(dn in actual,
                            'NCFile dimension {0!r} not present'.format(dn))
            self.assertEqual(actual[dn], dv,
                            'NCFile dimension {0!r} not correct'.format(dn))

    def test_nc4_NCFile_dimensions(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.dimensions
        expected = self.ncdims
        ncf.close()
        print_test_msg('NCFile.dimensions', actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected),
                         'NCFile dimensions not correct length')
        for dn, dv in expected.iteritems():
            self.assertTrue(dn in actual,
                            'NCFile dimension {0!r} not present'.format(dn))
            self.assertEqual(actual[dn], dv,
                            'NCFile dimension {0!r} not correct'.format(dn))

    def test_nio_NCFile_unlimited(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.unlimited('t')
        expected = True
        print_test_msg('NCFile.unlimited(t)', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile dimension t not unlimited')
        actual = ncf.unlimited('x')
        expected = False
        print_test_msg('NCFile.unlimited(x)', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile dimension x not limited')
        ncf.close()

    def test_nc4_NCFile_unlimited(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.unlimited('t')
        expected = True
        print_test_msg('NCFile.unlimited(t)', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile dimension t not unlimited')
        actual = ncf.unlimited('x')
        expected = False
        print_test_msg('NCFile.unlimited(x)', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile dimension x not limited')
        ncf.close()

    def test_nio_NCFile_ncattrs(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.ncattrs
        expected = self.ncattrs.keys()
        ncf.close()
        print_test_msg('NCFile.ncattrs', actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected),
                         'NCFile ncattrs not correct length')
        for dn in expected:
            self.assertTrue(dn in actual,
                            'NCFile ncattrs {0!r} not present'.format(dn))
            self.assertEqual(ncf.getncattr(dn), self.ncattrs[dn],
                            'NCFile ncattrs {0!r} not correct'.format(dn))

    def test_nc4_NCFile_ncattrs(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.ncattrs
        expected = self.ncattrs.keys()
        print_test_msg('NCFile.ncattrs', actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected),
                         'NCFile ncattrs not correct length')
        for xname in expected:
            self.assertTrue(xname in actual,
                            'NCFile ncattrs {0!r} not present'.format(xname))
            xval = self.ncattrs[xname]
            aval = ncf.getncattr(xname) 
            self.assertEqual(aval, xval,
                            'NCFile ncattrs {0!r} not correct'.format(xname))
        ncf.close()

    def test_nio_NCFile_variables(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables
        ncf.close()
        print_test_msg('NCFile.variables', actual=actual)
        self.assertEqual(len(actual), 3,
                         'NCFile variables not correct length')
        self.assertTrue('t' in actual,
                        't variable not in NCFile')
        self.assertTrue('x' in actual,
                        'x variable not in NCFile')
        self.assertTrue('v' in actual,
                        'v variable not in NCFile')
        for vn, vo in actual.iteritems():
            self.assertTrue(isinstance(vo, iobackend.NCVariable),
                            'Variable {0!r} has wrong type'.format(vn))

    def test_nc4_NCFile_variables(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables
        ncf.close()
        print_test_msg('NCFile.variables', actual=actual)
        self.assertEqual(len(actual), 3,
                         'NCFile variables not correct length')
        self.assertTrue('t' in actual,
                        't variable not in NCFile')
        self.assertTrue('x' in actual,
                        'x variable not in NCFile')
        self.assertTrue('v' in actual,
                        'v variable not in NCFile')
        for vn, vo in actual.iteritems():
            self.assertTrue(isinstance(vo, iobackend.NCVariable),
                            'Variable {0!r} has wrong type'.format(vn))

    def test_nio_NCFile_close(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.close()
        print_test_msg('NCFile.close', actual=actual)

    def test_nc4_NCFile_close(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.close()
        print_test_msg('NCFile.close', actual=actual)

    def test_nio_NCVariable_ncattrs(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].ncattrs
        expected = self.vattrs.keys()
        print_test_msg('NCVariable.ncattrs',
                       actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected),
                         'NCVariable ncattrs not correct length')
        for a in expected:
            self.assertTrue(a in actual,
                            'Attribute {0!r} not found in variable'.format(a))
            self.assertEqual(ncf.variables['v'].getncattr(a), self.vattrs[a],
                            'Attribute {0!r} not correct'.format(a))
        ncf.close()

    def test_nc4_NCVariable_ncattrs(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].ncattrs
        expected = self.vattrs.keys()
        print_test_msg('NCVariable.ncattrs',
                       actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected),
                         'NCVariable ncattrs not correct length')
        for a in expected:
            self.assertTrue(a in actual,
                            'Attribute {0!r} not found in variable'.format(a))
            self.assertEqual(ncf.variables['v'].getncattr(a), self.vattrs[a],
                            'Attribute {0!r} not correct'.format(a))
        ncf.close()

    def test_nio_NCVariable_dimensions(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].dimensions
        expected = ('t', 'x')
        ncf.close()
        print_test_msg('NCVariable.dimensions',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable dimensions not correct')

    def test_nc4_NCVariable_dimensions(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].dimensions
        expected = ('t', 'x')
        ncf.close()
        print_test_msg('NCVariable.dimensions',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable dimensions not correct')

    def test_nio_NCVariable_typecode(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].typecode
        expected = 'f'
        ncf.close()
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nc4_NCVariable_typecode(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].typecode
        expected = 'f'
        ncf.close()
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nio_NCVariable_shape(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].shape
        expected = (self.ncdims['t'], self.ncdims['x'])
        ncf.close()
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nc4_NCVariable_shape(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].shape
        expected = (self.ncdims['t'], self.ncdims['x'])
        ncf.close()
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nio_NCVariable_size(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].size
        expected = self.ncdims['t'] * self.ncdims['x']
        ncf.close()
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nc4_NCVariable_size(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'].size
        expected = self.ncdims['t'] * self.ncdims['x']
        ncf.close()
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nio_NCVariable_getitem(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'][:]
        expected = self.v[:]
        print_test_msg('NCVariable[v].__getitem__',
                       actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected)
        actual = ncf.variables['t'][:]
        expected = self.t[:]
        print_test_msg('NCVariable[t].__getitem__',
                       actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected)
        actual = ncf.variables['x'][:]
        expected = self.x[:]
        print_test_msg('NCVariable[x].__getitem__',
                       actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected)
        ncf.close()


    def test_nc4_NCVariable_getitem(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfrname)
        actual = ncf.variables['v'][:]
        expected = self.v[:]
        print_test_msg('NCVariable[v].__getitem__',
                       actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected)
        actual = ncf.variables['t'][:]
        expected = self.t[:]
        print_test_msg('NCVariable[t].__getitem__',
                       actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected)
        actual = ncf.variables['x'][:]
        expected = self.x[:]
        print_test_msg('NCVariable[x].__getitem__',
                       actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected)
        ncf.close()


#===============================================================================
# IOBackendWriteTests
#===============================================================================
class IOBackendWriteTests(unittest.TestCase):

    """
    IOBackendWriteTests Class

    This class defines all of the unit tests for the iobackend module.
    """
    
    def setUp(self):
        self.ncfwname = 'writetest.nc'
        self.ncattrs = {'a1': 'attribute 1',
                        'a2': 'attribute 2'}
        self.ncdims = {'t': 10, 'x': 5}
        self.t = np.arange(0, self.ncdims['t'], dtype='d')
        self.x = np.random.ranf(self.ncdims['x']).astype('d')
        self.v = np.random.ranf(self.ncdims['t']*self.ncdims['x']).reshape(10,5).astype('f')
        self.vattrs = {'long_name': 'variable',
                       'units': 'meters'}
    
    def tearDown(self):
        if exists(self.ncfwname):
            remove(self.ncfwname)

    def test_nio_NCFile_init_write(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        actual = type(ncf)
        ncf.close()
        expected = iobackend.NCFile
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nc4_NCFile_init_write(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        actual = type(ncf)
        ncf.close()
        expected = iobackend.NCFile
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nio_NCFile_setncattr(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for a,v in self.ncattrs.iteritems():
            ncf.setncattr(a, v)
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.attributes
        expected = self.ncattrs
        ncfr.close()
        print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
        for a,v in expected.iteritems():
            self.assertTrue(a in actual,
                            'NCFile attribute {0!r} not found'.format(a))
            self.assertEqual(actual[a], v,
                             'NCFile attribute {0!r} incorrect'.format(a))

    def test_nc4_NCFile_setncattr(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for a,v in self.ncattrs.iteritems():
            ncf.setncattr(a, v)
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.attributes
        expected = self.ncattrs
        ncfr.close()
        print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
        for a,v in expected.iteritems():
            self.assertTrue(a in actual,
                            'NCFile attribute {0!r} not found'.format(a))
            self.assertEqual(actual[a], v,
                             'NCFile attribute {0!r} incorrect'.format(a))

    def test_nio_NCFile_create_dimension(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('x', self.ncdims['x'])
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.dimensions['x']
        expected = self.ncdims['x']
        ncfr.close()
        print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile x-dimension incorrect')

    def test_nc4_NCFile_create_dimension(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('x', self.ncdims['x'])
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.dimensions['x']
        expected = self.ncdims['x']
        ncfr.close()
        print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile x-dimension incorrect')

    def test_nio_NCFile_create_dimension_unlimited(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.dimensions['t']
        expected = 0
        ncfr.close()
        print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile t-dimension incorrect')

    def test_nc4_NCFile_create_dimension_unlimited(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.dimensions['t']
        expected = 0
        ncfr.close()
        print_test_msg('NCFile.create_dimension()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile t-dimension incorrect')

    def test_nio_NCFile_create_variable(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('x', self.ncdims['x'])
        x = ncf.create_variable('x', 'd', ('x',))
        x[:] = self.x
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['x'][:]
        expected = self.x
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile x-variable incorrect')

    def test_nc4_NCFile_create_variable(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('x', self.ncdims['x'])
        x = ncf.create_variable('x', 'd', ('x',))
        x[:] = self.x
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['x'][:]
        expected = self.x
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile x-variable incorrect')

    def test_nio_NCFile_create_variable_unlimited(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        t = ncf.create_variable('t', 'd', ('t',))
        t[:] = self.t
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['t'][:]
        expected = self.t
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile t-variable incorrect')

    def test_nc4_NCFile_create_variable_unlimited(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        t = ncf.create_variable('t', 'd', ('t',))
        t[:] = self.t
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['t'][:]
        expected = self.t
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile t-variable incorrect')

    def test_nio_NCFile_create_variable_ndim(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        ncf.create_dimension('x', self.ncdims['x'])
        v = ncf.create_variable('v', 'f', ('t', 'x'))
        v[:] = self.v
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['v'][:]
        expected = self.v
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile 2d-variable incorrect')

    def test_nc4_NCFile_create_variable_ndim(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        ncf.create_dimension('x', self.ncdims['x'])
        v = ncf.create_variable('v', 'f', ('t', 'x'))
        v[:] = self.v
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['v'][:]
        expected = self.v
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile 2d-variable incorrect')

    def test_nio_NCVariable_setncattr(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        ncf.create_dimension('x', self.ncdims['x'])
        v = ncf.create_variable('v', 'f', ('t', 'x'))
        for attr,value in self.vattrs.iteritems():
            v.setncattr(attr, value)
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['v'].attributes
        expected = self.vattrs
        ncfr.close()
        print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
        for attr, value in expected.iteritems():
            self.assertTrue(attr in actual,
                            'Variable attribute {0!r} not found'.format(attr))
            self.assertEqual(actual[attr], value,
                             'Variable attribute {0!r} incorrect'.format(attr))

    def test_nc4_NCVariable_setncattr(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        ncf.create_dimension('x', self.ncdims['x'])
        v = ncf.create_variable('v', 'f', ('t', 'x'))
        for attr,value in self.vattrs.iteritems():
            v.setncattr(attr, value)
        ncf.close()
        ncfr = Nio.open_file(self.ncfwname)
        actual = ncfr.variables['v'].attributes
        expected = self.vattrs
        ncfr.close()
        print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
        for attr, value in expected.iteritems():
            self.assertTrue(attr in actual,
                            'Variable attribute {0!r} not found'.format(attr))
            self.assertEqual(actual[attr], value,
                             'Variable attribute {0!r} incorrect'.format(attr))


#===============================================================================
# IOBackendAppendTests
#===============================================================================
class IOBackendAppendTests(unittest.TestCase):

    """
    IOBackendAppendTests Class

    This class defines all of the unit tests for the iobackend module.
    """
    
    def setUp(self):
        self.ncfaname = 'appendtest.nc'
        self.ncattrs = {'a1': 'attribute 1',
                        'a2': 'attribute 2'}
        self.ncdims = {'t': 10, 'x': 5}
        self.t = np.arange(self.ncdims['t'], dtype='d')
        self.x = np.arange(self.ncdims['x'], dtype='d')
        self.v = np.arange(self.ncdims['t']*self.ncdims['x'],
                           dtype='f').reshape(self.ncdims['t'], self.ncdims['x'])
        self.vattrs = {'long_name': 'variable',
                       'units': 'meters'}
        
        self.fattrs2 = {'a3': 'attribute 3',
                        'a4': 'attribute 4'}
        self.t2 = np.arange(self.ncdims['t'], 2*self.ncdims['t'], dtype='d')
        self.v2 = np.arange(self.ncdims['t']*self.ncdims['x'],
                            dtype='f').reshape(self.ncdims['t'], self.ncdims['x'])
        self.vattrs2 = {'standard_name': 'variable'}

        ncfile = netCDF4.Dataset(self.ncfaname, 'w')
        for a,v in self.ncattrs.iteritems():
            setattr(ncfile, a, v)
        ncfile.createDimension('t')
        ncfile.createDimension('x', self.ncdims['x'])
        t = ncfile.createVariable('t', 'd', ('t',))
        t[:] = self.t
        x = ncfile.createVariable('x', 'd', ('x',))
        x[:] = self.x
        v = ncfile.createVariable('v', 'f', ('t', 'x'))
        for a,val in self.vattrs.iteritems():
            v.setncattr(a, val)
        v[:,:] = self.v
        
        ncfile.close()
    
    def tearDown(self):
        if exists(self.ncfaname):
            remove(self.ncfaname)

    def test_nio_NCFile_init_append(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        actual = type(ncf)
        ncf.close()
        expected = iobackend.NCFile
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nc4_NCFile_init_append(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        actual = type(ncf)
        ncf.close()
        expected = iobackend.NCFile
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nio_NCFile_setncattr(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        for a,v in self.fattrs2.iteritems():
            ncf.setncattr(a, v)
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.attributes
        expected = self.ncattrs
        expected.update(self.fattrs2)
        ncfr.close()
        print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
        for a,v in expected.iteritems():
            self.assertTrue(a in actual,
                            'NCFile attribute {0!r} not found'.format(a))
            self.assertEqual(actual[a], v,
                             'NCFile attribute {0!r} incorrect'.format(a))

    def test_nc4_NCFile_setncattr(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        for a,v in self.fattrs2.iteritems():
            ncf.setncattr(a, v)
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.attributes
        expected = self.ncattrs
        expected.update(self.fattrs2)
        ncfr.close()
        print_test_msg('NCFile.setncattr()', actual=actual, expected=expected)
        for a,v in expected.iteritems():
            self.assertTrue(a in actual,
                            'NCFile attribute {0!r} not found'.format(a))
            self.assertEqual(actual[a], v,
                             'NCFile attribute {0!r} incorrect'.format(a))

    def test_nio_NCFile_create_variable_ndim(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        v2 = ncf.create_variable('v2', 'f', ('t', 'x'))
        v2[:] = self.v2
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.variables['v2'][:]
        expected = self.v2
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile 2d-variable incorrect')

    def test_nc4_NCFile_create_variable_ndim(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        v2 = ncf.create_variable('v2', 'f', ('t', 'x'))
        v2[:] = self.v2
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.variables['v2'][:]
        expected = self.v2
        ncfr.close()
        print_test_msg('NCFile.create_variable()', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile 2d-variable incorrect')

    def test_nio_NCFile_variable_append(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        nt = self.ncdims['t']
        t = ncf.variables['t']
        t[nt:] = self.t2
        v = ncf.variables['v']
        v[nt:, :] = self.v2
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.variables['t'][:]
        expected = np.concatenate((self.t, self.t2))
        print_test_msg('NCVariable append', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile t-variable incorrect')
        actual = ncfr.variables['v'][:]
        expected = np.concatenate((self.v, self.v2))
        print_test_msg('NCVariable append', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile 2d-variable incorrect')
        ncfr.close()

    def test_nc4_NCFile_variable_append(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        nt = self.ncdims['t']
        t = ncf.variables['t']
        t[nt:] = self.t2
        v = ncf.variables['v']
        v[nt:, :] = self.v2
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.variables['t'][:]
        expected = np.concatenate((self.t, self.t2))
        print_test_msg('NCVariable append', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile t-variable incorrect')
        actual = ncfr.variables['v'][:]
        expected = np.concatenate((self.v, self.v2))
        print_test_msg('NCVariable append', actual=actual, expected=expected)
        npt.assert_array_equal(actual, expected,
                               'NCFile 2d-variable incorrect')
        ncfr.close()

    def test_nio_NCVariable_setncattr(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        v = ncf.variables['v']
        for attr,value in self.vattrs2.iteritems():
            v.setncattr(attr, value)
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.variables['v'].attributes
        expected = self.vattrs
        expected.update(self.vattrs2)
        ncfr.close()
        print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
        for attr, value in expected.iteritems():
            self.assertTrue(attr in actual,
                            'Variable attribute {0!r} not found'.format(attr))
            self.assertEqual(actual[attr], value,
                             'Variable attribute {0!r} incorrect'.format(attr))

    def test_nc4_NCVariable_setncattr(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        v = ncf.variables['v']
        for attr,value in self.vattrs2.iteritems():
            v.setncattr(attr, value)
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.variables['v'].attributes
        expected = self.vattrs
        expected.update(self.vattrs2)
        ncfr.close()
        print_test_msg('NCVariable.setncattr()', actual=actual, expected=expected)
        for attr, value in expected.iteritems():
            self.assertTrue(attr in actual,
                            'Variable attribute {0!r} not found'.format(attr))
            self.assertEqual(actual[attr], value,
                             'Variable attribute {0!r} incorrect'.format(attr))


#===============================================================================
# CLI
#===============================================================================
if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
