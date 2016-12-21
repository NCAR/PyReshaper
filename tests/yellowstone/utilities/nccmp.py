#!/usr/bin/env python
"""
NetCDF File Comparison Utility

Copyright 2016, University Corporation for Atmospheric Research
See the LICENSE.rst file for details
"""

import netCDF4 as nc
import optparse
import os

_USAGE_ = 'Usage:  %prog [options] FILE1 FILE2'
_DESC_ = """Compare two NetCDF files by comparing the header information and the
binary values of each field in the files.
"""
_PARSER_ = optparse.OptionParser(usage=_USAGE_, description=_DESC_)
_PARSER_.add_option('--header', action='store_true', default=False,
                    help='Compare only the header information and not the data values')

#==============================================================================
# Command-Line Operation
#==============================================================================
if __name__ == '__main__':
    opts, args = _PARSER_.parse_args()
    
    if len(args) != 2:
        raise ValueError('Only 2 arguments can be supplied')
    for i,f in enumerate(args):
        if not os.path.exists(f):
            raise ValueError('File {0} ({1}) does not exist'.format(i, f))
    FILE1, FILE2 = args
    
    # Record the differences
    diffs = {}
    
    # Open the files
    nc1 = nc.Dataset(FILE1)
    nc2 = nc.Dataset(FILE2)

    # Global Attribute Names
    gattdiffs = {}
    gatts1 = set(nc1.ncattrs())
    gatts2 = set(nc2.ncattrs())
    gatts1m2 = gatts1 - gatts2
    gatts2m1 = gatts2 - gatts1
    if len(gatts1m2) > 0:
        gattdiffs['names1'] = sorted(gatts1m2)
    if len(gatts1m2) > 0:
        gattdiffs['names2'] = sorted(gatts2m1)

    # Global Attribute Values
    gatts12 = set.intersection(gatts1, gatts2)
    gavldiffs = {}
    for gatt in gatts12:
        val1 = nc1.getncattr(gatt)
        val2 = nc2.getncattr(gatt)
        if val1 != val2:
            gavldiffs[gatt] = (val1, val2)
    if len(gavldiffs) > 0:
        gattdiffs['values'] = gavldiffs
    
    if len(gattdiffs) > 0:
        diffs['global'] = gattdiffs
    
    # Dimension Names
    dimdiffs = {}
    dims1 = set(nc1.dimensions.keys())
    dims2 = set(nc2.dimensions.keys())
    dims1m2 = dims1 - dims2
    dims2m1 = dims2 - dims1
    if len(dims1m2) > 0:
        dimdiffs['names1'] = sorted(dims1m2)
    if len(dims2m1) > 0:
        dimdiffs['names2'] = sorted(dims2m1)
    
    # Dimension Values
    dims12 = set.intersection(dims1, dims2)
    dvaldiffs = {}
    for dim in dims12:
        dim1 = nc1.dimensions[dim]
        dim2 = nc2.dimensions[dim]
        if len(dim1) != len(dim2) or dim1.isunlimited() != dim2.isunlimited():
            dvaldiffs[dim] = ((len(dim1), dim1.isunlimited()), (len(dim2), dim2.isunlimited()))
    if len(dvaldiffs) > 0:
        dimdiffs['values'] = dvaldiffs
    
    if len(dimdiffs) > 0:
        diffs['dimensions'] = dimdiffs
    
    # Variable Names
    vardiffs = {}
    vars1 = set(nc1.variables.keys())
    vars2 = set(nc2.variables.keys())
    vars1m2 = vars1 - vars2
    vars2m1 = vars2 - vars1
    if len(vars1m2) > 0:
        dimdiffs['names1'] = sorted(vars1m2)
    if len(vars2m1) > 0:
        dimdiffs['names2'] = sorted(vars2m1)

    print diffs
