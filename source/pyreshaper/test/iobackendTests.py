"""
Unit tests for the iobackend module

Copyright 2016, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

import unittest
import numpy as np
import numpy.testing as npt
import netCDF4

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
        self.ncfname = 'test.nc'
        self.ncattrs = {'a1': 'attribute 1',
                        'a2': 'attribute 2'}
        self.ncdims = {'t': 10, 'x': 5}
        self.t = np.arange(0, self.ncdims['t'], dtype='d')
        self.x = np.random.ranf(self.ncdims['x']).astype('d')
        self.v = np.random.ranf(self.ncdims['t']*self.ncdims['x']).reshape(10,5).astype('f')
        self.vattrs = {'long_name': 'variable',
                       'units': 'meters'}

        ncfile = netCDF4.Dataset(self.ncfname, 'w')
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
        if exists(self.ncfname):
            remove(self.ncfname)

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

    def test_nio_NCFile_init_read(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
        actual = type(ncf)
        expected = iobackend.NCFile
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nc4_NCFile_init_read(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
        actual = type(ncf)
        expected = iobackend.NCFile
        print_test_msg('NCFile.__init__()', actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCFile not created with correct type')

    def test_nio_NCFile_dimensions(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.dimensions
        expected = self.ncdims
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
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.dimensions
        expected = self.ncdims
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
        ncf = iobackend.NCFile(self.ncfname)
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

    def test_nc4_NCFile_unlimited(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
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
    

    def test_nio_NCFile_attributes(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.ncattrs
        expected = self.ncattrs.keys()
        print_test_msg('NCFile.ncattrs', actual=actual, expected=expected)
        self.assertEqual(len(actual), len(expected),
                         'NCFile ncattrs not correct length')
        for dn in expected:
            self.assertTrue(dn in actual,
                            'NCFile ncattrs {0!r} not present'.format(dn))
            self.assertEqual(ncf.getncattr(dn), self.ncattrs[dn],
                            'NCFile ncattrs {0!r} not correct'.format(dn))

    def test_nc4_NCFile_attributes(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
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

    def test_nio_NCFile_variables(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.variables
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
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.variables
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

    def test_nio_NCVariable(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.close()
        print_test_msg('NCFile.close', actual=actual)

    def test_nc4_NCFile_close(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.close()
        print_test_msg('NCFile.close', actual=actual)

    def test_nio_NCVariable_attributes(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
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

    def test_nc4_NCVariable_attributes(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
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

    def test_nio_NCVariable_dimensions(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.variables['v'].dimensions
        expected = ('t', 'x')
        print_test_msg('NCVariable.dimensions',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable dimensions not correct')

    def test_nc4_NCVariable_dimensions(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.variables['v'].dimensions
        expected = ('t', 'x')
        print_test_msg('NCVariable.dimensions',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable dimensions not correct')

    def test_nio_NCVariable_typecode(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.variables['v'].typecode
        expected = 'f'
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nc4_NCVariable_typecode(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
        actual = ncf.variables['v'].typecode
        expected = 'f'
        print_test_msg('NCVariable.typecode',
                       actual=actual, expected=expected)
        self.assertEqual(actual, expected,
                         'NCVariable typecode not correct')

    def test_nio_NCVariable_getitem(self):
        iobackend.set_backend('Nio')
        ncf = iobackend.NCFile(self.ncfname)
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


    def test_nc4_NCVariable_getitem(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfname)
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


#===============================================================================
# IOBackendWriteTests
#===============================================================================
class IOBackendWriteTests(unittest.TestCase):

    """
    IOBackendWriteTests Class

    This class defines all of the unit tests for the iobackend module.
    """
    
    def test_nio_NCFile_setncattr(self):
        iobackend.set_backend('Nio')
        ncf = iobackend


#===============================================================================
# CLI
#===============================================================================
if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
