"""
Copyright 2020, University Corporation for Atmospheric Research
See LICENSE.txt for details
"""

nlat = 19
nlon = 36
ntime = 10
nchar = 7

slices = ['input{0}.nc'.format(i) for i in range(5)]
scalars = ['scalar{0}'.format(i) for i in range(2)]
chvars = ['char{0}'.format(i) for i in range(1)]
timvars = ['tim{0}'.format(i) for i in range(2)]
xtimvars = ['tim{0}'.format(i) for i in range(2, 5)]
tvmvars = ['tvm{0}'.format(i) for i in range(2)]
tsvars = ['tsvar{0}'.format(i) for i in range(4)]
fattrs = {'attr1': 'attribute one', 'attr2': 'attribute two'}
