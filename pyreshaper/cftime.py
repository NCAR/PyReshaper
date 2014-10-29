'''
This is module of CF v1.6 time-specific compliant classes and helper
functions for dealing with CF-compliant time data in NetCDF files.

Presumably, most (if not all) of this functionality has already been
implemented in the cf-python package (http://cfpython.bitbucket.org/).
Hence, I leave it to future developers of the PyReshaper package to
decide if adding the extra dependencies required by the cf-python
package is worth the hastle.

My personal opinion is that this will be worth the effort if/when
(1) the required dependencies are already reliably installed on all necessary
platforms and machines and (2) Python usage throughout the community becomes
wide-spread enough to warrant in installation of the cf-python package.

Currently, the dependencies of the cf-python package are:

1. NetCDF4 v0.9.7+ (reliable everywhere)
2. UDUNITS-2 (reliable everywhere?)
3. Python v2.6+ (but < v3.0) (reliable everywhere?)
4. Python Numpy v1.6+ (reliable everywhere?)

_______________________
Created on Apr 30, 2014

@author Kevin Paul <kpaul@ucar.edu>
'''


#==============================================================================
# Factory function for creating CFCalendar objects by the CF convention name
#==============================================================================
def create_calendar(cal_attr=None):
    '''
    Factory function for converting string calendar specifications (i.e.,
    the value of the 'time:calendar' attribute) into a CFCalendar class.

    @param cal_attr: String containing the 'time:calendar' attribute value
                     (if any, could be None)
    @return A CFCalendar class object of the associated type
    '''


#==============================================================================
# CF Time Units Helper Class
#==============================================================================
class CFTimeUnits(object):
    '''
    A simple helper class for time unit conversion according to the CF 1.6
    convention.  This is limited to the 6 common units for time in UDUNITS:
    year, month, day, hour, minute, and second.
    '''

    def __init__(self):
        '''
        Constructor

        @return      An instance of a CFTimeUnits class
        '''

        ## Simple data structure allowing index-to-name conversion (using the
        # [index] syntax of the list class) and name-to-index conversion
        # (using the .index(name) syntax of the list class)
        self.__unit_names = ['year', 'month', 'day',
                             'hour', 'minute', 'second']

        ## List of all synonyms for units and their corresponding unit ID
        self.__unit_synonyms = {'yr': 0, 'year': 0, 'years': 0,
                                'month': 1, 'months': 1,
                                'day': 2, 'days': 2,
                                'h': 3, 'hr': 3, 'hrs': 3,
                                'hour': 3, 'hours': 3,
                                'min': 4, 'mins': 4, 'minute': 4, 'minutes': 4,
                                's': 5, 'sec': 5, 'secs': 5,
                                'second': 5, 'seconds': 5}

        ## Convertion table (first index is 'from' id, second index is 'to'
        #  id).  Floating point numbers indicate multiplication by that number,
        #  and text...
        self.__convertion_table = [[1., 12., 365.],
                                   [],
                                   [],
                                   [],
                                   [],
                                   []]

    def get_unit_name(self, units):
        '''
        Get the base name of the unit from a string or unit ID number

        @param units  Either a string containing a name of a unit (could be a
                      synonym), or an integer specifying the unit ID number
        '''


#==============================================================================
# CF Convention Date-Time Helper Class & Derived Types
#==============================================================================
class CFDateTime(object):
    '''
    This is a simple helper class to deal with date and time objects according
    to the CF 1.6 conventions.  This class is useful for converting between
    date-time offsets (e.g., date and time measured in units of days from a
    specific date and time) and date-time arrays (i.e., calendar year, month,
    and day and time in that day).
    '''

    def __init__(self, time_attr):
        '''
        Constructor

        @param time_attr  Dictionary-like data structure with the ability to
                          access attributes by string name.  Must contain the
                          'units' attribute.
        @return An instance of a CFDateTime class
        '''

        # Check for correct type, or raise error
        if (not hasattr(time_attr, '__getitem__')):
            raise TypeError('Time attributes must be dictionary-like')

        # Save the units_attr string for future use
        self.__units_attr_str = units_attr

        # Split the units attribute string into parts
        units_attr_parts = units_attr.split()

        # Extract the parts ("unit since YEAR-MONTH-DAY HOUR:MIN:SEC ZONE")
        # NOTE: Time zones are ignored in this class
        self.__unit = units_attr_parts[0]
        self.__date_start = map(int, units_attr_parts[2].split('-'))
        self.__time_start = [0.0, 0.0, 0.0]
        if (len(units_attr_parts) > 3):
            self.__time_start = map(float, units_attr_parts[3].split(':'))
        self.__timezone = 'UTC'
        if (len(units_attr_parts) > 4):
            self.__timezone = units_attr_parts[4]

        # Current datetime is initialized to the start date-time
        self.__date = list(self.__date_start)
        self.__time = list(self.__time_start)

    def units(self):
        '''
        Method for returning the units used for the date-time value

        @return A string specifying the units of the date-time values
        '''
        return self.__unit

    def year(self):
        '''
        Method for returning the datetime year

        @return The year of the date-time object
        '''
        return self.date[0]

    def month(self):
        '''
        Method for returning the datetime month

        @return: The month of the date-time object
        '''
        return self.date[1]

    def day(self):
        '''
        Method for returning the datetime day

        @return: The day of the date-time object
        '''
        return self.date[2]

    def hour(self):
        '''
        Method for returning the datetime hour

        @return: The hour of the date-time object
        '''
        return self.time[0]

    def minute(self):
        '''
        Method for returning the datetime minute

        @return: The minute of the date-time object
        '''
        return self.time[1]

    def second(self):
        '''
        Method for returning the datetime second

        @return: The second of the date-time object
        '''
        return self.time[2]

    def offset(self, dtoffset):
        '''
        Method for setting the date-time offset in units specified by the
        'time:units' attribute.

        @param dtoffset
        '''


#==============================================================================
# CF Convention Calendar Helper Class & Derived Types
#==============================================================================
class CFCalendar(object):
    '''
    CFCalendar Class

    This is a helper class containing the CF 1.6 specifications associated
    with the time:calendar attribute.  This is useful for determining the
    number of days in a month, leap years, etc., as defined by a given
    calendar.

    According to the CF 1.6 Convention, user-defined calendars are also an
    option.  This class also allows the definition of a user-defined
    calendar, per the CF convention specifications.

    The specifications of the base class correspond to the CF 1.6 calendar
    option, 'time:calendar = julian', which is the Julian calendar.  The
    Julian calendar has 12 months, each with the standard number of days
    in non-leap years.  A leap year is defined as any year evenly divisible
    by 4.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.__year_length = 365
        self.__month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        self.__leap_month = 2

    def is_leap_year(self, datetime):
        '''
        This function determines if a given year is a leap year, according to
        the calendar.  For the Julian calendar, any year evenly divisible by 4
        is a leap year.

        @param datetime  A CFDateTime object specifying the date and time
                         when a leap year is to be tested
        @return: True if the given year is a leap year, False otherwise.
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        if (datetime.year() % 4 == 0):
            return True
        else:
            return False

    def days_in_month(self, datetime):
        '''
        A function for computing the number of days in a given month of a
        given year.  This takes into account the leap days, if the given year
        is a leap year (according to the given calendar).

        @param datetime  A CFDateTime object specifying the date and time
                         when the number of days in the month is requested
        @return: Integer number of days in the datetime object's month
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        days = self.__month_lengths[datetime.month()]
        if (datetime.month() == self.__leap_month and
            self.is_leap_year(datetime)):
            days += 1
        return days

    def days_in_year(self, datetime):
        '''
        A method for computing the number of days in a given year, depending
        on whether the year is a leap year.

        @param datetime  A CFDateTime object specifying the date and time
                         when the number of days in the year is requested
        @return: Integer number of days in the datetime object's year
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        days = self.__year_length
        if (self.is_leap_year(datetime)):
            days += 1
        return days


class CFProlepticGregorianCalendar(CFCalendar):
    '''
    CFProlepticGregorianCalendar Class

    This class implements the 'time:calendar = proleptic_gregorian' CF 1.6
    calendar option.  This corresponds to a purely Gregorian calendar, assumed
    to exist prior to the 15 October 1582 inception.  For the Gregorian
    calendar, leap years are defined as:

        1. Any year evenly divisible by 4,
        2. BUT NOT if it is evenly divisible by 100,
        3. UNLESS it is also evenly divisible by 400.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(CFProlepticGregorianCalendar, self).__init__()

    def is_leap_year(self, datetime):
        '''
        This function determines if a given year is a leap year, according to
        the calendar.  For the Julian calendar, any year evenly divisible by 4
        is a leap year.

        @param datetime  A CFDateTime object specifying the date and time
                         when a leap year is to be tested
        @return: True if the given year is a leap year, False otherwise.
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        if (datetime.year() % 4 == 0):
            if (datetime.year() % 100 == 0):
                return (datetime.year() % 400 == 0)
            else:
                return True
        else:
            return False


class CFNoLeapCalendar(CFCalendar):
    '''
    CFNoLeapCalendar Class

    This class implements the 'time:calendar = noleap' (or, synonymously,
    'time:calendar = 365_day') CF 1.6 calendar option.  This corresponds to
    365-day year without leap years.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(CFNoLeapCalendar, self).__init__()

    def is_leap_year(self, datetime):
        '''
        This function determines if a given year is a leap year.  Since the
        'noleap' (or '365_day') calendars have no leap years, this always
        returns False.

        @param datetime  A CFDateTime object specifying the date and time
                         when a leap year is to be tested
        @return False
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        return False


class CFAllLeapCalendar(CFCalendar):
    '''
    CFAllLeapCalendar Class

    This class implements the 'time:calendar = allleap' (or, synonymously,
    'time:calendar = 366_day') CF 1.6 calendar option.  This corresponds to
    366-day year (i.e., every year is a leap year).
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(CFCalendar, self).__init__()

    def is_leap_year(self, datetime):
        '''
        This function determines if a given year is a leap year.  Since the
        'allleap' (or '366_day') calendars has only leap years, this always
        returns True.

        @param datetime  A CFDateTime object specifying the date and time
                         when a leap year is to be tested
        @return True
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        return True


class CF360DayCalendar(CFCalendar):
    '''
    CF360DayCalendar Class

    This class implements the 'time:calendar = 360_day' CF 1.6 calendar
    option.  This corresponds to a 360-day year (i.e., 30 days per month).
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(CF360DayCalendar, self).__init__()
        self.__year_length = 360
        self.__month_lengths = [30] * 12

    def is_leap_year(self, datetime):
        '''
        This function determines if a given year is a leap year.  Since the
        '360_day' calendar has no leap years, this always returns False.

        @param datetime A CFDateTime object specifying the date and time
                        when a leap year is to be tested
        @return False
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        return False


class CFGregorianCalendar(CFCalendar):
    '''
    CFGregorianCalendar Class

    This class implements the 'time:calendar = standard' or
    'time:calendar = gregorian' CF 1.6 calendar option.  This corresponds
    to the Julian calendar for dates before 15 Oct 1582, and the Gregorian
    calendar for all other dates.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(CFGregorianCalendar, self).__init__()
        self.__gregorian_inception_year = 1582
        self.__gregorian_inception_month = 10
        self.__gregorian_inception_day = 15

    def is_leap_year(self, datetime):
        '''
        This function determines if a given year is a leap year.  This obeys
        the rules of the Julian calendar for dates before 15 Oct 1582, and
        the rules of the Gregorian calendar for all other dates.

        @param datetime A CFDateTime object specifying the date and time
                        when a leap year is to be tested
        @return True if the year is a leap year, False otherwise.
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        if (datetime.year() < self.__gregorian_inception_year):
            return CFCalendar.is_leap_year(self, datetime)
        elif (datetime.year() == self.__gregorian_inception_year):
            if (datetime.month() < self.__gregorian_inception_month):
                return CFCalendar.is_leap_year(self, datetime)
            elif (datetime.month() == self.__gregorian_inception_month):
                if (datetime.day() < self.__gregorian_inception_day):
                    return CFCalendar.is_leap_year(self, datetime)
                else:
                    return CFProlepticGregorianCalendar.is_leap_year(self,
                                                                     datetime)
            else:
                return CFProlepticGregorianCalendar.is_leap_year(self,
                                                                 datetime)
        else:
            return CFProlepticGregorianCalendar.is_leap_year(self, datetime)


class CFNoneCalendar(CFCalendar):
    '''
    CFNoneCalendar Class

    This class implements the 'time:calendar = none' CF 1.6 calendar
    option.  This corresponds to a calendar with a constant repetition of the
    first day and time.  (This is like a periodic-bounded time, with the
    length of the time period equal to one 'unit' of time.)
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.__year_length = None
        self.__month_lengths = None
        self.__leap_month = None

    def is_leap_year(self, dtime):
        '''
        This function determines if a given year is a leap year.  Since the
        'none' calendar has no leap years, this always returns False.

        @param dtime: A datetime object specifying the date (and, optionally,
                      time) when a leap year is to be tested
        @return: False
        '''
        if (type(dtime) is not CFDateTime):
            raise TypeError('dtime must be a datetime object')
        return False

    def days_in_month(self, datetime):
        '''
        A function for computing the number of days in a given month of a
        given year.  For the 'none' calendar, this is returned as 0.

        @param datetime  A CFDateTime object specifying the date and time
                         when the number of months are requested
        @return 0
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a datetime object')
        return 0

    def days_in_year(self, datetime):
        '''
        A method for computing the number of days in a given year, depending
        on whether the year is a leap year.  For the 'none' calendar, this
        returns 0.

        @param datetime  A CFDateTime object specifying the date and time
                         when the number of days in the year is requested
        @return 0
        '''
        if (type(datetime) is not CFDateTime):
            raise TypeError('datetime must be a CFDateTime object')
        return 0
