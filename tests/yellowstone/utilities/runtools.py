#==============================================================================
#
#  TestTools
#
#  This is a collection of functions that are useful for running the PyReshaper
#  tests on the Yellowstone compute system.
#
#==============================================================================

# Builtin Modules
import glob

# Third-Party Modules
import numpy as np
import Nio

# Package Modules
from pyreshaper import specification


#==============================================================================
# Script/Job Runner for Yellowstone
#==============================================================================
class Runner(object):
    """
    A simple class for running jobs, writing submission scripts, etc
    """
    
    def __init__(self):