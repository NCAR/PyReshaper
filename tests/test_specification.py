"""
Copyright 2020, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

import os
import pickle

import pytest

from pyreshaper import specification

cwd = os.path.dirname(os.path.realpath(__file__))


def test_init():
    spec = specification.Specifier()
    assert len(spec.input_file_list) == 0
    assert spec.netcdf_format == 'netcdf4'
    assert spec.compression_level == 0
    assert spec.least_significant_digit is None
    assert spec.output_file_prefix == 'tseries.'
    assert spec.output_file_suffix == '.nc'
    assert spec.time_series is None
    assert len(spec.time_variant_metadata) == 0
    assert spec.assume_1d_time_variant_metadata is False
    assert spec.io_backend == 'netCDF4'
    assert spec.exclude_list == []
    assert spec.metadata_filename is None


def test_init_full():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf4c'
    cl = 4
    prefix = 'pre.'
    suffix = '.suf.nc'
    tseries = ['1', '2']
    metadata = ['x', 'y', 'z']
    xlist = ['g', 'h']
    meta1d = True
    metafile = 'd'
    backend = 'Nio'
    lsigfig = 3
    spec = specification.Specifier(
        infiles=in_list,
        ncfmt=fmt,
        compression=cl,
        prefix=prefix,
        suffix=suffix,
        timeseries=tseries,
        metadata=metadata,
        meta1d=meta1d,
        metafile=metafile,
        backend=backend,
        least_significant_digit=lsigfig,
        exclude_list=xlist,
    )
    for i1, i2 in zip(spec.input_file_list, in_list):
        assert i1 == i2
    assert spec.io_backend == backend
    assert spec.metadata_filename == metafile
    assert spec.netcdf_format == fmt
    assert spec.compression_level == cl
    assert spec.output_file_prefix == prefix
    assert spec.output_file_suffix == suffix
    assert spec.exclude_list == xlist
    assert spec.least_significant_digit == lsigfig
    for i1, i2 in zip(spec.time_series, tseries):
        assert i1 == i2
    for i1, i2 in zip(spec.time_variant_metadata, metadata):
        assert i1 == i2
    assert spec.assume_1d_time_variant_metadata == meta1d


def test_validate_types_defaults():
    in_list = ['a', 'b', 'c']
    spec = specification.Specifier(infiles=in_list)
    spec.validate_types()


def test_validate_types():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = 3
    prefix = 'pre.'
    suffix = '.suf.nc'
    tseries = ['1', '2']
    metadata = ['x', 'y', 'z']
    spec = specification.Specifier(
        infiles=in_list,
        ncfmt=fmt,
        compression=cl,
        prefix=prefix,
        suffix=suffix,
        timeseries=tseries,
        metadata=metadata,
        meta1d=True,
    )
    spec.validate_types()


def test_validate_types_fail_input():
    in_list = ['a', 2, 'c']
    fmt = 'netcdf'
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = ['x', 'y', 'z']
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_backend():
    in_list = ['a', 'b', 'c']
    fmt = 2342
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = ['x', 'y', 'z']
    backend = 1
    spec = specification.Specifier(
        infiles=in_list,
        ncfmt=fmt,
        compression=cl,
        prefix=prefix,
        suffix=suffix,
        metadata=metadata,
        backend=backend,
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_format():
    in_list = ['a', 'b', 'c']
    fmt = 2342
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = ['x', 'y', 'z']
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_cl():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = '6'
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = ['x', 'y', 'z']
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_prefix():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = 5
    prefix = dict()
    suffix = '.suf.nc'
    metadata = ['x', 'y', 'z']
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_suffix():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = 5
    prefix = 'pre.'
    suffix = list()
    metadata = ['x', 'y', 'z']
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_timeseries():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    tseries = ['1', 2.5]
    metadata = ['x', 'y', 'z']
    spec = specification.Specifier(
        infiles=in_list,
        ncfmt=fmt,
        compression=cl,
        prefix=prefix,
        suffix=suffix,
        timeseries=tseries,
        metadata=metadata,
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_metadata():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = ['x', 'y', 2]
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_types_fail_meta1d():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = ['x', 'y', 'z']
    meta1d = 't'
    spec = specification.Specifier(
        infiles=in_list,
        ncfmt=fmt,
        compression=cl,
        prefix=prefix,
        suffix=suffix,
        meta1d=meta1d,
        metadata=metadata,
    )
    with pytest.raises(TypeError):
        spec.validate_types()


def test_validate_values_fail_input():
    in_list = ['a', 'b', 'c']
    fmt = 'netcdf'
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = []
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    spec.validate_types()
    with pytest.raises(ValueError):
        spec.validate_values()


def test_validate_values_fail_backend():
    in_list = ['test_reshaper.py', 'test_s2srun.py']
    fmt = 'netcdf9'
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = []
    backend = 'x'
    spec = specification.Specifier(
        infiles=in_list,
        ncfmt=fmt,
        compression=cl,
        prefix=prefix,
        suffix=suffix,
        metadata=metadata,
        backend=backend,
    )
    spec.validate_types()
    with pytest.raises(ValueError):
        spec.validate_values()


def test_validate_values_fail_format():
    in_list = ['test_reshaper.py', 'test_s2srun.py']
    fmt = 'netcdf9'
    cl = 5
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = []
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    spec.validate_types()
    with pytest.raises(ValueError):
        spec.validate_values()


def test_validate_values_fail_cl():
    in_list = ['test_reshaper.py', 'test_s2srun.py']
    fmt = 'netcdf4'
    cl = 111
    prefix = 'pre.'
    suffix = '.suf.nc'
    metadata = []
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, compression=cl, prefix=prefix, suffix=suffix, metadata=metadata
    )
    spec.validate_types()
    with pytest.raises(ValueError):
        spec.validate_values()


def test_validate_values_fail_prefix():
    in_list = ['test_reshaper.py', 'test_s2srun.py']
    fmt = 'netcdf4'
    prefix = '/sfcsrytsdfv/pre.'
    suffix = '.suf.nc'
    metadata = []
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata
    )
    spec.validate_types()
    with pytest.raises(ValueError):
        spec.validate_values()


def test_validate_values_suffix():
    in_list = [cwd + '/test_specification.py']
    fmt = 'netcdf4'
    prefix = 'pre.'
    suffix = '.suf'
    metadata = []
    spec = specification.Specifier(
        infiles=in_list, ncfmt=fmt, prefix=prefix, suffix=suffix, metadata=metadata
    )
    spec.validate_types()
    spec.validate_values()
    assert spec.output_file_suffix == suffix + '.nc'


def test_write():
    in_list = ['test_specification.py']
    fmt = 'netcdf4'
    cl = 8
    prefix = 'pre.'
    suffix = '.suf.nc'
    tseries = ['1', '2']
    metadata = ['time']
    spec = specification.Specifier(
        infiles=in_list,
        ncfmt=fmt,
        compression=cl,
        prefix=prefix,
        suffix=suffix,
        timeseries=tseries,
        metadata=metadata,
    )
    fname = 'test_write.s2s'
    spec.write(fname)
    assert os.path.exists(fname), 'Specfile failed to write'
    spec2 = pickle.load(open(fname, 'rb'))
    for i1, i2 in zip(spec2.input_file_list, in_list):
        assert i1 == i2
    assert spec2.netcdf_format == fmt
    assert spec2.compression_level == cl
    assert spec2.output_file_prefix == prefix
    assert spec2.output_file_suffix == suffix
    for i1, i2 in zip(spec2.time_series, tseries):
        assert i1 == i2
    for i1, i2 in zip(spec2.time_variant_metadata, metadata):
        assert i1 == i2
    os.remove(fname)
