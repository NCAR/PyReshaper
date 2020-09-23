"""
Copyright 2020, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

from functools import reduce
from os import remove
from os.path import exists

import netCDF4
import Nio
import numpy as np
import numpy.testing as npt
import pytest

from pyreshaper import iobackend


@pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
def test_set_get_backend(backend):
    iobackend.set_backend(backend)
    npt.assert_equal(iobackend.get_backend(), backend)


def test_set_backend_failure():
    with pytest.raises(KeyError):
        iobackend.set_backend('x')


class ReadTests:
    @pytest.fixture(autouse=True)
    def data(self):
        self.ncfrname = 'readtest.nc'
        self.ncattrs = {'a1': 'attribute 1', 'a2': 'attribute 2'}
        self.ncdims = {'t': 10, 'x': 5, 'c': 14}
        self.vdims = {'t': ('t',), 'x': ('x',), 'v': ('t', 'x'), 's': ('c',)}
        self.vdtype = {'t': 'd', 'x': 'd', 'v': 'f', 's': 'c'}
        self.vshape = {
            't': (self.ncdims['t'],),
            'x': (self.ncdims['x'],),
            'v': (self.ncdims['t'], self.ncdims['x']),
            's': (self.ncdims['c'],),
        }
        self.ncvars = {
            't': np.arange(0, self.ncdims['t'], dtype=self.vdtype['t']),
            'x': np.random.ranf(self.ncdims['x']).astype(self.vdtype['x']),
            'v': np.random.ranf(self.ncdims['t'] * self.ncdims['x'])
            .reshape(10, 5)
            .astype(self.vdtype['v']),
            's': np.array([c for c in 'this is a stri']).astype(self.vdtype['s']),
        }
        self.vattrs = {
            't': {'long_name': u'time', 'units': u'days'},
            'x': {'long_name': u'space', 'units': u'meters'},
            'v': {'long_name': u'variable', 'units': u'kg', 'missing_value': np.float32(1e20)},
            's': {'long_name': u'string'},
        }
        self.vfill = {'t': None, 'x': None, 'v': np.float32(1e20), 's': None}

        ncfile = netCDF4.Dataset(self.ncfrname, 'w')
        for a in self.ncattrs:
            setattr(ncfile, a, self.ncattrs[a])
        for d in self.ncdims:
            if d == 't':
                ncfile.createDimension(d)
            else:
                ncfile.createDimension(d, self.ncdims[d])
        for v in self.ncvars:
            vobj = ncfile.createVariable(v, self.vdtype[v], self.vdims[v], fill_value=self.vfill[v])
            for a in self.vattrs[v]:
                vobj.setncattr(a, self.vattrs[v][a])
            vobj[:] = self.ncvars[v]
        ncfile.close()

        yield

        if exists(self.ncfrname):
            remove(self.ncfrname)

    def test_NCFile_init_mode_failure(self):
        with pytest.raises(ValueError):
            iobackend.NCFile(filename=self.ncfrname, mode='x')

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_init(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        assert isinstance(ncf, iobackend.NCFile)
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_dimensions(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        assert ncf.dimensions == self.ncdims
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_unlimited(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        assert ncf.unlimited('t') is True
        assert ncf.unlimited('x') is False
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_ncattrs(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        assert ncf.ncattrs == list(self.ncattrs.keys())
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_variables(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in ncf.variables:
            assert isinstance(ncf.variables[v], iobackend.NCVariable)
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_variable_data(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in ncf.variables:
            npt.assert_array_equal(ncf.variables[v][:], self.ncvars[v])
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_ncattrs_getncattr(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in self.vattrs:
            vattrs = {a: ncf.variables[v].getncattr(a) for a in ncf.variables[v].ncattrs}
            assert vattrs == self.vattrs[v]
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_dimensions(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in self.vdtype:
            vdims = ncf.variables[v].dimensions
            assert vdims == self.vdims[v]
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_datatype(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in self.vdtype:
            vdtype = ncf.variables[v].datatype
            assert vdtype == np.dtype(self.vdtype[v])

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_shape(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in self.vshape:
            vshape = ncf.variables[v].shape
            assert vshape == self.vshape[v]
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_fill_value(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in self.vfill:
            vfill = ncf.variables[v].fill_value
            assert vfill == self.vfill[v]
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_size(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in self.vdims:
            vsize1 = ncf.variables[v].size
            vsize2 = reduce(lambda x, y: x * y, (self.ncdims[d] for d in self.vdims[v]), 1)
            assert vsize1 == vsize2
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_getitem(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfrname)
        for v in self.vdims:
            npt.assert_array_equal(ncf.variables[v][:], self.ncvars[v])
        ncf.close()


class WriteTests:
    @pytest.fixture(autouse=True)
    def data(self):
        self.ncfwname = 'writetest.nc'
        self.ncattrs = {'a1': 'attribute 1', 'a2': 'attribute 2'}
        self.ncdims = {'t': 10, 'x': 5, 'c': 14, 'n': 2}
        self.vdims = {'t': ('t',), 'x': ('x',), 'v': ('t', 'x'), 's': ('n', 'c'), 'c': ('c',)}
        self.vdtype = {'t': 'd', 'x': 'd', 'v': 'f', 's': 'c', 'c': 'c'}
        self.vshape = {v: tuple(self.ncdims[d] for d in self.vdims[v]) for v in self.vdims}
        self.ncvars = {
            't': np.arange(0, self.ncdims['t'], dtype=self.vdtype['t']),
            'x': np.random.ranf(self.ncdims['x']).astype(self.vdtype['x']),
            'v': np.random.ranf(self.ncdims['t'] * self.ncdims['x'])
            .reshape(10, 5)
            .astype(self.vdtype['v']),
            's': np.array(['a string', 'another string'], dtype='S14')
            .view('S1')
            .reshape(self.vshape['s']),
            'c': np.array('scalar str', dtype='S14').reshape(1).view('S1'),
        }
        self.vattrs = {
            't': {'long_name': u'time', 'units': u'days'},
            'x': {'long_name': u'space', 'units': u'meters'},
            'v': {'long_name': u'variable', 'units': u'kg', 'missing_value': np.float32(1e20)},
            's': {'long_name': u'vector of strings'},
            'c': {'long_name': u'scalar string'},
        }
        self.vfill = {'t': None, 'x': None, 'v': np.float32(1e20), 's': None, 'c': None}
        self.chunks = {'t': [5], 'x': None, 'v': [2, 3], 's': None, 'c': None}

        yield

        if exists(self.ncfwname):
            remove(self.ncfwname)

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_init_write(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        assert isinstance(ncf, iobackend.NCFile)
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_setncattr(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for a in self.ncattrs:
            ncf.setncattr(a, self.ncattrs[a])
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        attrs = {a: ncfr.getncattr(a) for a in ncfr.ncattrs}
        ncfr.close()
        assert attrs == self.ncattrs

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_create_dimension(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for d in self.ncdims:
            ncf.create_dimension(d, self.ncdims[d])
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        dims = ncfr.dimensions
        ncfr.close()
        assert dims == self.ncdims

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_create_dimension_unlimited(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        ncf.create_dimension('t')
        ncf.create_dimension('x', 3)
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        assert ncfr.dimensions['t'] == 0
        assert ncfr.unlimited('t') is True
        assert ncfr.dimensions['x'] == 3
        assert ncfr.unlimited('x') is False
        ncfr.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_create_variable(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for d in self.ncdims:
            if d == 't':
                ncf.create_dimension(d)
            else:
                ncf.create_dimension(d, self.ncdims[d])
        for v in self.ncvars:
            ncf.create_variable(v, self.vdtype[v], self.vdims[v])
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        for v in self.ncvars:
            assert isinstance(ncfr.variables[v], iobackend.NCVariable)
        ncfr.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_create_variable_fill(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for d in self.ncdims:
            if d == 't':
                ncf.create_dimension(d)
            else:
                ncf.create_dimension(d, self.ncdims[d])
        for v in self.ncvars:
            ncf.create_variable(v, self.vdtype[v], self.vdims[v], fill_value=self.vfill[v])
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        for v in self.ncvars:
            vfill = None if v == 's' else self.vfill[v]
            assert ncfr.variables[v].fill_value == vfill
        ncfr.close()

    def test_NCFile_create_variable_chunksizes(self):
        iobackend.set_backend('netCDF4')
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for d in self.ncdims:
            if d == 't':
                ncf.create_dimension(d)
            else:
                ncf.create_dimension(d, self.ncdims[d])
        for v in self.ncvars:
            ncf.create_variable(v, self.vdtype[v], self.vdims[v], chunksizes=self.chunks[v])
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        for v in self.ncvars:
            chunksize = self.chunks[v] if self.chunks[v] else 'contiguous'
            assert ncfr.variables[v].chunk_sizes == chunksize
        ncfr.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_setncattr(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for d in self.ncdims:
            if d == 't':
                ncf.create_dimension(d)
            else:
                ncf.create_dimension(d, self.ncdims[d])
        for v in self.vattrs:
            vobj = ncf.create_variable(v, self.vdtype[v], self.vdims[v])
            for attr in self.vattrs[v]:
                vobj.setncattr(attr, self.vattrs[v][attr])
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        for v in self.vattrs:
            vattrs = {a: ncfr.variables[v].getncattr(a) for a in ncfr.variables[v].ncattrs}
            assert vattrs == self.vattrs[v]
        ncfr.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_setitem_getitem(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfwname, mode='w')
        for d in self.ncdims:
            if d == 't':
                ncf.create_dimension(d)
            else:
                ncf.create_dimension(d, self.ncdims[d])
        for v in self.ncvars:
            vobj = ncf.create_variable(v, self.vdtype[v], self.vdims[v])
            vobj[:] = self.ncvars[v]
        ncf.close()
        ncfr = iobackend.NCFile(self.ncfwname)
        for v in self.ncvars:
            npt.assert_array_equal(ncfr.variables[v][:], self.ncvars[v])
        ncfr.close()


class AppendTests:
    @pytest.fixture(autouse=True)
    def data(self):
        self.ncfaname = 'appendtest.nc'
        self.ncattrs = {'a1': 'attribute 1', 'a2': 'attribute 2'}
        self.ncdims = {'t': 10, 'x': 5}
        self.t = np.arange(self.ncdims['t'], dtype='d')
        self.x = np.arange(self.ncdims['x'], dtype='d')
        self.v = np.arange(self.ncdims['t'] * self.ncdims['x'], dtype='f').reshape(
            self.ncdims['t'], self.ncdims['x']
        )
        self.vattrs = {'long_name': 'variable', 'units': 'meters'}
        self.fattrs2 = {'a3': 'attribute 3', 'a4': 'attribute 4'}
        self.t2 = np.arange(self.ncdims['t'], 2 * self.ncdims['t'], dtype='d')
        self.v2 = np.arange(self.ncdims['t'] * self.ncdims['x'], dtype='f').reshape(
            self.ncdims['t'], self.ncdims['x']
        )
        self.vattrs2 = {'standard_name': 'variable'}

        ncfile = netCDF4.Dataset(self.ncfaname, 'w')
        for a, v in self.ncattrs.items():
            setattr(ncfile, a, v)
        ncfile.createDimension('t')
        ncfile.createDimension('x', self.ncdims['x'])
        t = ncfile.createVariable('t', 'd', ('t',))
        t[:] = self.t
        x = ncfile.createVariable('x', 'd', ('x',))
        x[:] = self.x
        v = ncfile.createVariable('v', 'f', ('t', 'x'))
        for a, val in self.vattrs.items():
            v.setncattr(a, val)
        v[:, :] = self.v

        ncfile.close()

        yield

        if exists(self.ncfaname):
            remove(self.ncfaname)

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_init_append(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        assert isinstance(ncf, iobackend.NCFile)
        ncf.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_setncattr(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        for a, v in self.fattrs2.items():
            ncf.setncattr(a, v)
        ncf.close()
        ncfr = netCDF4.Dataset(self.ncfaname)
        actual = {a: ncfr.getncattr(a) for a in ncfr.ncattrs()}
        ncfr.close()
        expected = self.ncattrs
        expected.update(self.fattrs2)
        assert actual == expected

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_create_variable_ndim(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        v2 = ncf.create_variable('v2', np.dtype('f'), ('t', 'x'))
        v2[:] = self.v2
        ncf.close()
        ncfr = netCDF4.Dataset(self.ncfaname)
        actual = ncfr['v2'][:]
        ncfr.close()
        npt.assert_array_equal(actual, self.v2)

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCFile_variable_append(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        nt = self.ncdims['t']
        t = ncf.variables['t']
        t[nt:] = self.t2
        v = ncf.variables['v']
        v[nt:, :] = self.v2
        ncf.close()
        ncfr = netCDF4.Dataset(self.ncfaname)
        actual = ncfr.variables['t'][:]
        expected = np.concatenate((self.t, self.t2))
        npt.assert_array_equal(actual, expected)
        actual = ncfr.variables['v'][:]
        expected = np.concatenate((self.v, self.v2))
        npt.assert_array_equal(actual, expected)
        ncfr.close()

    @pytest.mark.parametrize('backend', iobackend._AVAILABLE_)
    def test_NCVariable_setncattr(self, backend):
        iobackend.set_backend(backend)
        ncf = iobackend.NCFile(self.ncfaname, mode='a')
        v = ncf.variables['v']
        for attr, value in self.vattrs2.items():
            v.setncattr(attr, value)
        ncf.close()
        ncfr = Nio.open_file(self.ncfaname)
        actual = ncfr.variables['v'].attributes
        expected = self.vattrs
        expected.update(self.vattrs2)
        ncfr.close()
        assert actual == expected
