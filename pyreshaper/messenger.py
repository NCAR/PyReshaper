'''
This is the Parallel Utilities class to provide basic load-balancing and
decomposition routines for the PyReshaper.  This encapsulates all of the
mpi4py calls as needed.  It also has the ability to detect if the mpi4py
module is present.  If not, it throws an exception.

_______________________
Created on Apr 30, 2014

@author: Kevin Paul <kpaul@ucar.edu>
'''

import sys


#==============================================================================
# create_messenger factory function
#==============================================================================
def create_messenger(serial=False):
    '''
    This is the factory function the creates the necessary decomposition
    object for parallel (or serial) operation.  The type must be specified
    by the user with the 'serial' argument, as there is not portable way
    of determining if the run should be assumed to be serial or parallel
    from the environment.

    @param serial  True or False, indicating whether the serial or parallel
                   decomposition utility object should be constructed and
                   returned.  DEFAULT: False (parallel operation)

    @return  A decomposition utility object.
    '''
    # Check type
    if (type(serial) is not bool):
        err_msg = "The serial argument must be a bool."
        raise TypeError(err_msg)

    # Construct and return the desired decomp utility
    if (serial):
        return Messenger()
    else:
        return MPIMessenger()


#==============================================================================
# Messenger Base Class
#==============================================================================
class Messenger(object):
    '''
    This is the base class for decomposition/parallel utilities.  This defines
    serial operation, and has no dependencies.  These methods are reimplemented
    in the derived class for parallel decomposition.
    '''

    def __init__(self):
        '''
        Constructor
        '''

        ## The type of decomp utility constructed
        self.messenger_type = 'serial'

        ## Whether this is the master process/rank
        self._is_master = True

        ## The rank of the processor
        self._mpi_rank = 0

        ## Size of the MPI communicator
        self._mpi_size = 1

        ## Indicates verbosity level
        self.verbosity = 1

    def partition(self, global_list):
        '''
        This method returns a subset of a list to be "taken" by the local
        process.  For serial execution, just returns the global list.

        @param global_list  A list of values to be distributed across all
                            processors

        @return A subset of the global list for use on the current processor
        '''
        if (not isinstance(global_list, list)):
            err_msg = "List of global variables has wrong type"
            raise TypeError(err_msg)

        return global_list

    def sync(self):
        '''
        A wrapper on the MPI Barrier method.  Forces execution to wait for
        all other processors to get to this point in the code.  Does nothing
        in serial.
        '''
        return

    def is_master(self):
        '''
        Returns True or False depending on whether this rank is the master
        rank (i.e., rank 0).  In serial, always returns True.

        @return True or False depending on whether this rank is the master
                rank (i.e., rank 0).
        '''
        return self._is_master

    def get_rank(self):
        '''
        Returns the integer associated with the local MPI processor/rank.
        In serial, it always returns 0.

        @return The integer ID associated with the local MPI processor/rank
        '''
        return self._mpi_rank

    def get_size(self):
        '''
        Returns the number of ranks in the MPI communicator.
        In serial, it always returns 1.

        @return The integer number of ranks in the MPI communicator
        '''
        return self._mpi_size

    def sum(self, data):
        '''
        This sums data across all processors.

        If the data is a dictionary (with assumed numberic values), then this
        summes the values associated with each key across all processors.  It
        returns a dictionary with the sum for each key.  Every processor
        must have the same keys in their associated dictionaries.

        If the data is not a dictionary, but is list-like (i.e., iterable),
        then this returns a single value with the sum of all values across all
        processors.

        If the data is not iterable (i.e., a single value), then it sums across
        all processors and returns one value.

        In serial, this returns just the sum on the local processor.

        @param data The data with values to be summed

        @return The sum of the data values
        '''
        if (isinstance(data, dict)):
            totals = {}
            for name in data:
                totals[name] = self.sum(data[name])
            return totals
        elif (hasattr(data, '__getitem__')):
            return reduce(lambda x, y: x + y, data)
        else:
            return data

    def max(self, data):
        '''
        This returns the maximum of the data across all processors.

        If the data is a dictionary (with assumed numberic values), then this
        find the maximum of the values associated with each key across all
        processors.  It returns a dictionary with the maximum for each key.
        Every processor must have the same keys in their associated
        dictionaries.

        If the data is not a dictionary, but is list-like (i.e., iterable),
        then this returns a single value with the maximum of all values across
        all processors.

        If the data is not iterable (i.e., a single value), then it finds the
        maximum across all processors and returns one value.

        In serial, this returns just the maximum on the local processor.

        @param data The data with values

        @return The maximum of the data values
        '''
        if (isinstance(data, dict)):
            maxima = {}
            for name in data:
                maxima[name] = self.max(data[name])
            return maxima
        elif (hasattr(data, '__getitem__')):
            return max(data)
        else:
            return data

    def print_once(self, output, vlevel=0):
        '''
        This method prints output to stdout, but only if it is the
        master processes and the verbosity level is high enough.

        @param output A string that should be printed to stdout

        @param vlevel The verbosity level associated with the message.  If
                      this level is less than the messenger's verbosity, no
                      output is generated.
        '''
        if (self.is_master() and vlevel < self.verbosity):
            print output
            sys.stdout.flush()

    def print_all(self, output, vlevel=0):
        '''
        This prints a message to stdout from every processor.

        @param output A string that should be printed to stdout

        @param vlevel The verbosity level associated with the message.  If
                      this level is less than the messenger's verbosity, no
                      output is generated.
        '''
        if (vlevel < self.verbosity):
            ostr = '[' + str(self._mpi_rank) + '/' \
                       + str(self._mpi_size) + '] ' + output
            print ostr
            sys.stdout.flush()


#==============================================================================
# MPIMessenger Class
#==============================================================================
class MPIMessenger(Messenger):
    '''
    This is the parallel-operation class for decomposition/parallel utilities.
    This is derived from the Messenger class, which defines the serial
    operation.  This defines operations using MPI.
    '''

    def __init__(self):
        '''
        Constructor
        '''

        ## Type of decomp utility constructed
        self.messenger_type = 'parallel'

        # Try to import the MPI module
        try:
            from mpi4py import MPI
        except:
            raise ImportError('Failed to import MPI.')

        ## Pointer to the MPI module
        self._mpi = MPI

        ## The rank of the processor
        self._mpi_rank = self._mpi.COMM_WORLD.Get_rank()

        ## MPI Communicator size
        self._mpi_size = self._mpi.COMM_WORLD.Get_size()

        ## Whether this is the master process/rank
        self._is_master = (self._mpi_rank == 0)

    def partition(self, global_list):
        '''
        This method returns a subset of a list to be "taken" by the local
        process.  For serial execution, just returns the global list.

        @param global_list  A list of values to be distributed across all
                            processors

        @return A subset of the global list for use on the current processor
        '''

        if (not isinstance(global_list, list)):
            err_msg = "List of global variables has wrong type"
            raise TypeError(err_msg)

        local_list = global_list[self._mpi_rank::self._mpi_size]

        return local_list

    def sync(self):
        '''
        A wrapper on the MPI Barrier method.  Forces execution to wait for
        all other processors to get to this point in the code.
        '''
        self._mpi.COMM_WORLD.Barrier()

    def sum(self, data):
        '''
        This sums data across all processors.

        If the data is a dictionary (with assumed numberic values), then this
        summes the values associated with each key across all processors.  It
        returns a dictionary with the sum for each key.  Every processor
        must have the same keys in their associated dictionaries.

        If the data is not a dictionary, but is list-like (i.e., iterable),
        then this returns a single value with the sum of all values across all
        processors.

        If the data is not iterable (i.e., a single value), then it sums across
        all processors and returns one value.

        In serial, this returns just the sum on the local processor.

        @param data The data with values to be summed

        @return The sum of the data values
        '''
        if (isinstance(data, dict)):
            totals = {}
            for name in data:
                totals[name] = self.sum(data[name])
            return totals
        elif (hasattr(data, '__len__')):
            total = Messenger.sum(self, data)
            return self.sum(total)
        else:
            return self._mpi.COMM_WORLD.allreduce(data, op=self._mpi.SUM)

    def max(self, data):
        '''
        This returns the maximum of the data across all processors.

        If the data is a dictionary (with assumed numberic values), then this
        find the maximum of the values associated with each key across all
        processors.  It returns a dictionary with the maximum for each key.
        Every processor must have the same keys in their associated
        dictionaries.

        If the data is not a dictionary, but is list-like (i.e., iterable),
        then this returns a single value with the maximum of all values across
        all processors.

        If the data is not iterable (i.e., a single value), then it finds the
        maximum across all processors and returns one value.

        In serial, this returns just the maximum on the local processor.

        @param data The data with values

        @return The maximum of the data values
        '''
        if (isinstance(data, dict)):
            maxima = {}
            for name in data:
                maxima[name] = self.max(data[name])
            return maxima
        elif (hasattr(data, '__getitem__')):
            maximum = Messenger.max(self, data)
            return self.max(maximum)
        else:
            return self._mpi.COMM_WORLD.allreduce(data, op=self._mpi.MAX)
