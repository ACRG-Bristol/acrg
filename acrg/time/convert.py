# -*- coding: utf-8 -*-
"""

FUNCTIONS

sec2time: Convert "seconds since YYYY-MM-DD HH:MM" format (e.g. from CF files)
    to datetime

    Example: 
        time = \
        acrg_time.convert.sec2time(seconds, "2000-01-01 00:00")

time2sec: Calculate seconds since some reference time. Assumes you want 
    a reference time at the beginning of that year, unless you specify a 
    particular reference time

    Example:
        seconds_since, reference = acrg_time.convert.time2sec(datetime)


Created on Fri Nov 21 10:48:30 2014

@author: chxmr
"""
from __future__ import division
import sys
if sys.version_info[0] == 2:
    from builtins import str

import datetime as dt
import time as tm
import calendar
import dateutil
from matplotlib.dates import (julian2num, num2date)

def check_iter(var):
    
    if not hasattr(var, '__iter__'):
        var = [var]
        notIter = True
    else:
        notIter = False

    return var, notIter


def return_iter(var, notIter):
    
    if notIter:
        return var[0]
    else:
        return var


def reftime(time_reference):
    
    time_reference, notIter = check_iter(time_reference)
    time_reference = return_iter(time_reference, notIter)
    #If reference time is a string, assume it's in CF convention 
    # and convert to datetime
    #if type(time_reference[0]) is str or type(time_reference[0]) is str:
    if isinstance(time_reference,str):
        time_ref=dateutil.parser.parse(time_reference)
    else:
        time_ref=time_reference
    
    return time_ref


def sec2time(seconds, time_reference):

    seconds, notIter = check_iter(seconds)

    time_ref = reftime(time_reference)

    return return_iter([time_ref + 
        dt.timedelta(seconds=int(s)) for s in seconds], notIter)


def min2time(minutes, time_reference):
    
    minutes, notIter = check_iter(minutes)
    
    time_ref = reftime(time_reference)

    return return_iter([time_ref +
        dt.timedelta(minutes=m) for m in minutes], notIter)

def hours2time(hours, time_reference):
    
    hours, notIter = check_iter(hours)
    
    time_ref = reftime(time_reference)

    return return_iter([time_ref +
        dt.timedelta(hours=m) for m in hours], notIter)

def day2time(days, time_reference):
    
    days, notIter = check_iter(days)

    time_ref = reftime(time_reference)
    
    return return_iter([time_ref + dt.timedelta(days=d) for d in days], 
                        notIter)


def time2sec(time, time_reference=None):
    
    time, notIter = check_iter(time)
    
    if time_reference is None:
        time_reference=dt.datetime(min(time).year, 1, 1, 0, 0)

    time_seconds=[\
        (t.replace(tzinfo=None)-time_reference).total_seconds() \
        for t in time]

    return return_iter(time_seconds, notIter), time_reference


def time2decimal(dates):
    
    def sinceEpoch(date): # returns seconds since epoch
        return tm.mktime(date.timetuple())
    s = sinceEpoch
    
    dates, notIter = check_iter(dates)
    
    frac=[]
    for date in dates:
        year = date.year
        startOfThisYear = dt.datetime(year=year, month=1, day=1)
        startOfNextYear = dt.datetime(year=year+1, month=1, day=1)
    
        yearElapsed = s(date) - s(startOfThisYear)
        yearDuration = s(startOfNextYear) - s(startOfThisYear)
        fraction = old_div(yearElapsed,yearDuration)

        frac.append(date.year + fraction)
        
    return return_iter(frac, notIter)
    
    
def decimal2time(frac):
    
    frac, notIter = check_iter(frac)

    dates = []
    for f in frac:
        year = int(f)
        yeardatetime = dt.datetime(year, 1, 1)
        daysPerYear = 365 + calendar.leapdays(year, year+1)
        dates.append(yeardatetime + dt.timedelta(days = daysPerYear*(f - year)))
    
    return return_iter(dates, notIter)


def julian2time(dates):
    '''
    Convert Julian dates (e.g. from IDL) to datetime
    '''
    
    dates, notIter = check_iter(dates)

    dates_julian = []
    for date in dates:
        dates_julian.append(num2date(julian2num(date)))
    
    return return_iter(dates_julian, notIter)
        
def convert_to_hours(time):
    '''
    Convert to hours
    
    Returns in the input provided, float or list of floats
    '''
    hours_per_unit = {"H": 1, "D": 24, "W": 168, "M": 732, "Y":8760}
    if type(time) is list:
        time_hrs_list = []
        for ii in range(len(time)):
            time_hrs = float(time[ii][:-1]) * hours_per_unit[time[ii][-1]]
            time_hrs_list.append(time_hrs)   
        return time_hrs_list    
    else:
        time_hrs = float(time[:-1]) * hours_per_unit[time[-1]]
        return time_hrs