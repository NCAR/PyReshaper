#!/bin/bash

# Necessary modules to load
module load python
module load all-python-libs

# Assume this will be run from the test directory itself
# (i.e., tests/yellowstone/serial)
PYRDIR_REL=../../..
if [ ! -d "$PYRDIR_REL/pyreshaper" ]; then
  echo "Cannot find PyReshaper package."
  exit 1
fi

# Turn the relative path into an absolute path
PYRDIR=`readlink -f $PYRDIR_REL`
PYREXE=$PYRDIR/scripts/slice2series

# Save the old Python path and add the reshaper to it
OLD_PYTHONPATH=$PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$PYRDIR

# Set the run directory
RUNDIR=$PYRDIR/tests/yellowstone/serial
if [ ! -d "$RUNDIR" ]; then
  echo "Serial rundir does not exist."
  exit 1
fi

# Set the output directory (where output files are written)
OUTDIR=$RUNDIR/tseries
if [ -d "$OUTDIR" ]; then
  rm -rf $OUTDIR
fi
mkdir $OUTDIR

# Move to the run directory
cd $RUNDIR

# Define the input parameters to the PyReshaper script
INDIR=/glade/u/tdd/asap/data/b.e12.B1850C5CN.ne30_g16.init.ch.027/ocn/hist
INGLOB1=${INDIR}/b.e12.B1850C5CN.ne30_g16.init.ch.027.pop.h.000*.nc
INGLOB2=${INDIR}/b.e12.B1850C5CN.ne30_g16.init.ch.027.pop.h.0010*.nc
FORMAT='netcdf4c'
PREFIX="$OUTDIR/b.e12.B1850C5CN.ne30_g16.init.ch.027.pop.h."
SUFFIX='.000101-001012.nc'
MOPTS='-m time -m time_bound'

# Run the script in serial
python $PYREXE --serial --limit 3 -g "$INGLOB1" -g "$INGLOB2" -f $FORMAT -p $PREFIX -s $SUFFIX $MOPTS

# Set the Python path back to the original
export PYTHONPATH=$OLD_PYTHONPATH
