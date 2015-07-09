#==============================================================================
#
#  TestTools
#
#  This is a collection of functions that are useful for running the PyReshaper
#  tests on the Yellowstone compute system.
#
#==============================================================================

# Builtin Modules
import os
import stat
from subprocess import Popen, PIPE, STDOUT


#==============================================================================
# Script/Job Runner for Yellowstone
#==============================================================================
class Job(object):
    """
    A simple class for running jobs, writing submission scripts, etc
    """

    def __init__(self, filename="runscript.sh", nodes=0, tiling=16,
                 minutes=120, queue="small", pcode="STDD0002", runcmds=None):
        """
        Constructor

        Parameters:
            filename (str): Name of the run script
            nodes (int): Number of nodes to request
            tiling (int): Number of processors per node to request
            minutes (int): Number of walltime minutes to request
            queue (str): Name of queue to submit the job
            runcmds (list): The list of commands to run
        """
        self.set(filename=filename, nodes=nodes, tiling=tiling,
                 minutes=minutes, queue=queue, pcode=pcode, runcmds=runcmds)

    def set(self, filename=None, nodes=None, tiling=None, minutes=None,
            queue=None, pcode=None, runcmds=None):
        """
        Set the job specific parameters

        Parameters:
            filename (str): Name of the run script
            nodes (int): Number of nodes to request
            tiling (int): Number of processors per node to request
            minutes (int): Number of walltime minutes to request
            queue (str): Name of queue to submit the job
            runcmds (list): The list of commands to run
        """
        if filename:
            self._filename = str(filename)
            self._rootname, runext = self._filename.rsplit('.', 1)
            self._logext = 'log' if runext != 'log' else 'out'
            self._rundir = os.path.dirname(os.path.realpath(self._filename))
            self._written = False
        if nodes:
            self._nodes = int(nodes)
            self._written = False
        if tiling:
            self._tiling = int(tiling)
            self._written = False
        if minutes:
            self._minutes = int(minutes)
            self._written = False
        if queue:
            self._queue = str(queue)
            self._written = False
        if pcode:
            self._pcode = str(pcode)
            self._written = False
        if runcmds:
            self._written = False
            self._runcmds = []
            if isinstance(runcmds, str):
                self._runcmds.append(runcmds)
            elif isinstance(runcmds, list):
                strcmds = [str(cmd) for cmd in runcmds]
                self._runcmds.extend(strcmds)

    def write(self):
        """
        Write the run script associated with this job

        Parameters:
            filename (str): Name of the run script file to write
        """
        # Check if written already
        if self._written:
            return

        # Check that run commands have been set
        if not self._runcmds:
            err_msg = "Cannot write run script without run commands"
            raise ValueError(err_msg)

        # Start creating the run scripts for each test
        run_script_list = ['#!/bin/bash', '']

        # If necessary, add the parallel preamble
        if (self._nodes > 0):

            # Number of processors total
            num_procs = self._nodes * self._tiling

            # Generate walltime in string form
            wtime_hours = self._minutes / 60
            if (wtime_hours > 99):
                wtime_hours = 99
                print 'Requested number of hours too large.  Limiting to', \
                    wtime_hours, '.'
            wtime_minutes = self._minutes % 60
            wtime_str = '%02d:%02d' % (wtime_hours, wtime_minutes)

            # String list representing LSF preamble
            run_script_list.extend([
                '#BSUB -n ' + str(num_procs),
                '#BSUB -R "span[ptile=' + str(self._tiling) + ']"',
                '#BSUB -q ' + self._queue,
                '#BSUB -a poe',
                '#BSUB -x',
                '#BSUB -o ' + self._rootname + '.%J.' + self._logext,
                '#BSUB -J ' + self._rootname,
                '#BSUB -P ' + self._code,
                '#BSUB -W ' + wtime_str,
                '',
                'export MP_TIMEOUT=14400',
                'export MP_PULSE=1800',
                'export MP_DEBUG_NOTIMEOUT=yes',
                ''])

        # Now create the rest of the run script
        run_script_list.extend(['# Necessary modules to load',
                                'module load python',
                                'module load all-python-libs',
                                ''])

        run_script_list.extend(self._runcmds)
        run_script_list.append('')

        # Write the script to file
        run_script_file = open(self._filename, 'w')
        run_script_file.write(os.linesep.join(run_script_list))
        run_script_file.close()

    def launch(self):
        """
        Launch the job (or submit into the queue)
        """

        # Write run script, if not done yet
        if not self._written:
            self.write()

        # Now launch the test
        cwd = os.getcwd()
        os.chdir(self._rundir)
        if (self._nodes == 0):

            # Make the script executable
            os.chmod(self._filename,
                     stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH)

            # Launch the serial job as a subprocess
            job = Popen([run_script], stdout=PIPE, stderr=STDOUT,
                        env=os.environ.copy())
            pid = str(job.pid)

            # Wait for job to finish and grab job output
            job_output = job.communicate()[0]

            # Write output to log file
            log_file = open(self._rootname + '.' + self._logext, 'w')
            log_file.write(job_output)
            log_file.close()

        else:

            # Open up the run script for input to LSF's bsub
            run_script_file = open(self._filename, 'r')

            # Launch the parallel job with LSF bsub
            job = Popen(['bsub'], stdout=PIPE, stderr=STDOUT,
                        stdin=run_script_file, env=os.environ.copy())

            # Grab the bsub output
            job_output = job.communicate()[0]

            # Close the script file and print submission info
            run_script_file.close()

        os.chdir(cwd)
