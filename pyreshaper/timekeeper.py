'''
This is a simple class to act as a time keeper for internal
performance monitoring (namely, timing given processes).

_______________________
Created on May 15, 2014

@author: Kevin Paul <kpaul@ucar.edu>
'''

import time


class TimeKeeper(object):
    '''
    Class to keep timing recordings, start/stop/reset timers.
    '''

    def __init__(self):
        '''
        Constructor
        '''

        ## Dictionary of start times associated with a string name
        self._start_times = {}

        ## Dictionary of accumulated times associated with a string name
        self._accumulated_times = {}

        ## List containing the order of the timers
        #  (when added to the dictionaries)
        self._added_order = []

    def reset(self, name):
        '''
        Method to reset a timer associated with a given name.  If the name
        has never been used before, the timer is created and the accumulated
        time is set to 0.  If the timer has been used before, the accumulated
        time is set to 0.

        @param name  The name or ID of the timer to reset
        '''

        # Reset the named timer (creates it if it doesn't exist yet)
        if (name not in self._added_order):
            self._added_order.append(name)
        self._accumulated_times[name] = 0.0
        self._start_times[name] = time.time()

    def start(self, name):
        '''
        Method to start a timer associated with a given name.  If the name
        has never been used before, the timer is created and the accumulated
        time is set to 0.

        @param name  The name or ID of the timer to start
        '''

        # Start the named timer (creates it if it doesn't exist yet)
        if (name not in self._accumulated_times):
            self.reset(name)
        else:
            self._start_times[name] = time.time()

    def stop(self, name):
        '''
        Method to stop a timer associated with a given name, and adds
        the accumulated time to the timer when stopped.  If the given timer
        name has never been used before (either by calling reset() or start()),
        the timer is created and the accumulated time is set to 0.

        @param name  The name or ID of the timer to stop
        '''

        # Stop the named timer, add to accumulated time
        if (name not in self._accumulated_times):
            self.reset(name)
        else:
            self._accumulated_times[name] += \
                time.time() - self._start_times[name]

    def get_time(self, name):
        '''
        Returns the accumulated time of the given timer.  If the given timer
        name has never been created, it is created and the accumulated time
        is set to zero before returning.

        @param name  The name or ID of the timer to stop

        @return  The accumulated time of the named timer (or 0.0 if the
                 named timer has never been created before).
        '''

        # Get the accumulated time
        if (name not in self._accumulated_times):
            self.reset(name)
        return self._accumulated_times[name]

    def print_times(self):
        '''
        Returns the accumulated times of all of the timers in a string
        with easy-to-read formating.  Printed in the order in which the
        timers were created.

        @return  An easy-to-read string with the accumulated times of
                 all of the timers
        '''

        # Header line with Title
        hline = '-' * 40 + '\n'
        ostr = hline + ' TIMING DATA:\n' + hline

        # Determine the longest timer name
        # and thus computer the column to line up values
        valcol = 0
        for name in self._added_order:
            if (len(str(name)) > valcol):
                valcol = len(str(name))
        valcol += 2

        # Print out the timer names and the accumulated times
        for name in self._added_order:
            spacer = ' ' * (valcol - len(str(name)))
            ostr += str(name) + ':' + spacer
            ostr += str(self._accumulated_times[name]) + '\n'

        return ostr
