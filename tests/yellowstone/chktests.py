#!/usr/bin/env python
#
#  This script is designed to run cprnc netCDF file comparisons for
#  a given set of ncReshaper tests.  It is to be run in parallel,
#  and therefore must be submitted to the queue.  This script first
#  checks if the filenames match between the "test output" and the 
#  "correct output".  If the filenames match, then it proceeds to 
#  compare each
#
#======================================================================

import os
_path = os.path.dirname(os.path.realpath(__file__)) + "/../python"
import sys
sys.path.append(_path)
import glob
import string
import subprocess
import testdata as td
from mpi4py import MPI

#_______________________________________________________________________
# Basic MPI information
#
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()

#_______________________________________________________________________
# Usage information
#
def usage():
  '''General usage information'''
  if (mpi_rank == 0):
    print "Usage: "+sys.argv[0]+" filename [testname] rootdir"
    print
    print "Options:"
    print
    print "    filename     The name of the test data file."
    print 
    print "    testname     The name of the test (in the test data file)"
    print "                 to compare.  Only compares the one test."
    print
    print "    rootdir      The root directory in which all tests exist."
    print
  
#_______________________________________________________________________
# Parse command-line arguments
#
fname = ""
tname = ""
rootd = ""
cmpOneTest = False
if (len(sys.argv) == 3):
  fname = sys.argv[1]
  rootd = sys.argv[2]
  if (mpi_rank == 0):
    print "Will compare all test results."
elif (len(sys.argv) == 3):
  fname = sys.argv[1]
  tname = sys.argv[2]
  rootd = sys.argv[3]
  cmpOneTest = True
  if (mpi_rank == 0):
    print "Will compare test results for only test name '"+tname+"'."
else:
  usage()
  sys.exit(1)

#_______________________________________________________________________
# Validate user input
#
if (not os.path.isfile(fname)):
  if (mpi_rank == 0):
    print "ERROR:  Test data file (" + fname + ") not found."
  sys.exit(1)

rootd = os.path.abspath(rootd)
if (not os.path.isdir(rootd)):
  if (mpi_rank == 0):
    print "ERROR:  Root directory (" + rootd + ") not found."
  sys.exit(1)

#_______________________________________________________________________
# Move to root directory
#
os.chdir(rootd)

#_______________________________________________________________________
# Read and validate the test data file
#
allTests = td.TestData(fname)
allTests.makeValid()

#_______________________________________________________________________
# Create the list of tests to compare
#
cmpIndices = []
if (cmpOneTest):
  # Test if testname is valid
  if (tname in allTests.test_names):
    cmpIndices = [allTests.test_names.index(tname)]
  else:
    if (mpi_rank == 0):
      print "ERROR: Test name '"+tname+"' not valid.  Not found in test data file."
    sys.exit(1)
else:
  # Compare all valid tests
  cmpIndices = allTests.complete_ind[:]

#_______________________________________________________________________
# Set the cprnc tool for netCDF comparison
#
cprnc_args = ["/glade/p/work/kpaul/installs/cprnc/bin/cprnc", "-m", "", ""]
ident_str = "files seem to be IDENTICAL"

#_______________________________________________________________________
# Loop over all tests and perform comparisons
#
num_tests = len(allTests.test_names)
num_failures = [0] * num_tests
num_files = [0] * num_tests
set_matches = [False] * num_tests

for i in cmpIndices:
  
  test_name = allTests.test_names[i]
  test_path = rootd + "/" + test_name
  
  if (not os.path.exists(test_path)):
    if (mpi_rank == 0):
      print "Test directory ("+test_path+") not found.  Skipping comparison."
  else:
    if (mpi_rank == 0):
      print "Comparing test '"+test_name+"':"
    cmp_dir = test_path + "/compare"
    if (not os.path.exists(cmp_dir)):
      if (mpi_rank == 0):
        print "   Comparison directory not found in test directory.  Skipping comparison."
    else:
      if (mpi_rank == 0):
        print "   Placing cprnc comparison failures in '" + cmp_dir + "'."

      # Define the paths to the original (old) series and new (new) series
      oldpath = allTests.series_paths[i]
      newpath = rootd + "/" + allTests.result_relpaths[i]

      # Series files and result (new series) files      
      os.chdir(oldpath)
      oldls = glob.glob("*.nc")
      os.chdir(newpath)
      newls = glob.glob("*.nc")
      
      # Record the total number of files
      num_files[i] = len(oldls)
      
      # Second comparison test: Match filenames 
      if (set(oldls) == set(newls)):
        set_matches[i] = True
        
      # If sets match, then compare files directly
      if (set_matches[i]):
        
        # Loop over file index, so striding can be done later in parallel
        for j in range(mpi_rank, len(oldls), mpi_size):
          ncfname = oldls[j]
          ncfold = oldpath + "/" + ncfname
          ncfnew = newpath + "/" + ncfname
          cprnc_args[2:] = [ncfold, ncfnew]
          cprnc_out = subprocess.check_output(cprnc_args)
          if (string.rfind(cprnc_out, ident_str) < 0):
            num_failures[i] = num_failures[i] + 1
            cmp_path = cmp_dir + "/" + ncfname + ".cprnc" 
            cmp_file = open(cmp_path,"w")
            cmp_file.write(cprnc_out)
            cmp_file.close()

# Wait for all tests to complete before gathering
mpi_comm.Barrier()
if (mpi_rank == 0):
  print 
  print "Comparison of all tests complete."

#_______________________________________________________________________
# Gather the number of failures on each processor
#
# set_matches and num_files are the same on all procs, so only
# need to gather number of failures on each procs for each test

num_fail_all = mpi_comm.gather(num_failures, root=0)
if (mpi_rank == 0):
  num_fail_transpose = map(list, zip(*num_fail_all))
  for i in range(len(num_fail_transpose)):
    num_failures[i] = sum(num_fail_transpose[i])

#_______________________________________________________________________
# Now that all tests are done, write out summary to root proc
#
if (mpi_rank == 0):
  print
  for i in cmpIndices:
    print "Test: "+allTests.test_names[i]
    if (set_matches[i]):
      print "   GOOD: Series file names match."
      if (num_failures[i] > 0):
        print "   BAD: " + str(num_failures[i]) + \
              " (of " + str(num_files[i]) + " total) file differences."
      else:
        print "   GOOD: No file differences found."
    else:
      print "   BAD: Series file names do not match."
      print "   BAD: Series files not compared due to name mismatch."

#EOF
