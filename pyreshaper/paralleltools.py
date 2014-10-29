'''
This is the Parallel Utilities class to provide basic load-balancing and
decomposition routines for the PyReshaper.  This encapsulates all of the
mpi4py calls as needed.  It also has the ability to detect if the mpi4py
module is present.  If not, it throws an exception.

_______________________
Created on Apr 30, 2014

@author: Kevin Paul <kpaul@ucar.edu>
'''


#==============================================================================
# create_partool Factory function
#==============================================================================
def create_partool(serial=False):
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
        return ParallelTool()
    else:
        return MPIParallelTool()


#==============================================================================
# ParallelTool Base Class
#==============================================================================
class ParallelTool(object):
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
        self.decomp_type = 'serial'

        ## Whether this is the master process/rank
        self._is_master = True

        ## The rank of the processor
        self._mpi_rank = 0

        ## Size of the MPI communicator
        self._mpi_size = 1

        ## Indicates whether debugging is enabled, or not
        self._debugging = False

    def enable_debugging(self):
        '''
        Turns on debugging options in parallel.
        '''
        self._debugging = True

    def partition(self, global_vars_names):
        '''
        Takes a dictionary of variables (with the keys being string names of
        the variables) and extracts a subset of the variables to be used by
        the local processor.

        For serial execution, this returns the global variables dictionary.

        @param global_vars_names  A dictionary of global time-series variables

        @return A dictionary of local time-series variables (a subset of the
                global variables dictionary)
        '''
        if (not isinstance(global_vars_names, list)):
            err_msg = "List of global variables has wrong type"
            raise TypeError(err_msg)

        return global_vars_names

    def sync(self):
        '''
        A wrapper on the MPI Barrier method.  Forces execution to wait for
        all other processors to get to this point in the code.

        Does nothing in serial.
        '''
        return

    def is_master(self):
        '''
        Returns True or False depending on whether this rank is the master
        rank (e.g., rank 0).

        In serial, always returns True.

        @return True or False depending on whether this rank is the master
                rank (e.g., rank 0).
        '''
        return self._is_master

    def get_rank(self):
        '''
        Returns the number associated with the local MPI processor/rank.

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
        Sums all of the data across all processors.

        In serial, just sums the data passed to this method.

        @param data The data with values to be summed (a dictionary)

        @return The sum of the data values
        '''
        if (not isinstance(data, dict)):
            err_msg = "Data must be a dictionary"
            raise TypeError(err_msg)
        return reduce(lambda x, y: x + y, data.values())

    def max(self, data):
        '''
        Finds the maximum value of each element in a set of data,
        across all processors.

        In serial, just returns the data passed in.

        @param data The dictionary of data in which the maximum is to be found

        @return A list with the maximum value of each element in the data
        '''
        return data

    def master_print(self, output):
        '''
        This method prints output to stdout, but only if it is the
        master processes.

        @param output A string that should be printed to stdout
        '''
        if (self.is_master()):
            print output

    def debug_print(self, output, allranks=False):
        '''
        This prints a debug message to stdout from the specified ranks,
        if debugging is enabled.

        @param output A string that should be printed to stdout

        @param ranks A list of ranks from which output should be printed
        '''
        if (self._debugging):
            ostr = str(self._mpi_rank) + '/' \
                 + str(self._mpi_size) + ': ' + output
            if (allranks or self.is_master()):
                print ostr


#==============================================================================
# MPIParallelTool Class
#==============================================================================
class MPIParallelTool(ParallelTool):
    '''
    This is the parallel-operation class for decomposition/parallel utilities.
    This is derived from the ParallelTool class, which defines the serial
    operation.
    '''

    def __init__(self):
        '''
        Constructor
        '''

        ## Type of decomp utility constructed
        self.decomp_type = 'parallel'

        # Try to import the MPI module
        try:
            from mpi4py import MPI
            self._mpi = MPI
        except:
            raise ImportError('Failed to import MPI.')

        ## The rank of the processor
        self._mpi_rank = self._mpi.COMM_WORLD.Get_rank()

        ## MPI Communicator size
        self._mpi_size = self._mpi.COMM_WORLD.Get_size()

        ## Whether this is the master process/rank
        self._is_master = (self._mpi_rank == 0)

    def partition(self, global_vars_names):
        '''
        Takes a list of variable names and extracts a subset of the
        variables to be used by the local processor.

        @param global_vars_names  A list of global time-series variable names

        @return A list of local time-series variable names (a subset of the
                global variables list)
        '''
        if (not isinstance(global_vars_names, list)):
            err_msg = "List of global variables has wrong type"
            raise TypeError(err_msg)

        rank = self._mpi_rank
        size = self._mpi_size

        local_vars_names = global_vars_names.keys()[rank::size]

        return local_vars_names

    def sync(self):
        '''
        A wrapper on the MPI Barrier method.  Forces execution to wait for
        all other processors to get to this point in the code.
        '''
        self._mpi.COMM_WORLD.Barrier()

    def sum(self, data):
        '''
        Sums all of the data across all processors.

        @param data The data with values to be summed (a dictionary)

        @return The sum of the data values
        '''
        total = ParallelTool.sum(self, data)
        return self._mpi.COMM_WORLD.allreduce(total, op=self._mpi.SUM)

    def max(self, data):
        '''
        Finds the maximum value of each element in a set of data,
        across all processors.

        @param data The list of data in which the maximum is to be found

        @return A list with the maximum value of each element in the data
        '''
        return self._mpi.COMM_WORLD.allreduce(data, op=self._mpi.MAX)
