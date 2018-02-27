"""
Author: matth
Date Created: 28/02/2017 12:06 PM
"""

from __future__ import division

class wells_db(object):
    """
    class to hold input kwargs for rd.sql
    """
    monthly = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'DTW_MYaverage'}
    yearly = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'DTW_YEARLY_AVERAGES_VIEW'}
    daily = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'DTW_READINGS'}
    base ={'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'VIEW_WELL'}
    well_details = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'WELL_DETAILS'}
    WMCR_Zones = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'WMCR_Zones'}
    aquifer_tests = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'AQUIFER_TESTS'}
    screen_details = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'SCREEN_DETAILS'}
    allo_zones = {'server': 'SQL2012PROD05', 'database': 'Wells', 'table': 'AllocationZones'}