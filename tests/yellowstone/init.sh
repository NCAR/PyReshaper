#!/bin/bash

# Function to return the absolute path
function abspath {
  cd $1
  echo `pwd`
}

# Add the PyReshaper executable scripts to the path
BIN_DIR=`abspath ../../bin`
export PATH=$PATH:$BIN_DIR

# Add the PyReshaper module path to the python path
MOD_DIR=`abspath ../../`
export PYTHONPATH=$PYTHONPATH:$MOD_DIR
