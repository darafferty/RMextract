#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 11:47:31 2018

@author: mevius
"""

import os
from RMextract import getIONEX as ionex
import numpy as np
from datetime import date
import EMM.EMM as EMM
import math


def getRM(*args,**kwargs):
    
    '''Dummy function for compatibility'''
    logging.warning("getRM is depreciated, please use extract(\"RM\") instead")
    RMextract.extract("RM",*args,**kwargs)
    

def _bookkeeping():
    '''Define which function to call with what parameters'''
    pass



def _get_RM(times, pointing, position, tecinfo, apply_earth_rotation=0):
    '''Get interpolated vTEC, piercepoint ocordinates and parallel B-Field for an array of times
    
    optionally correcting for earth rotation.
    Args:
        times (np.array) : times in MJD seconds.
        pointing (np.array) : ra, dec in radians.
        position (np.rray) : antenna xyz position (ITRF) in meters.
        tecinfo Tuple[np.array, np.array, np.array, np.array, np.array] :
            output of getIonex.read_tec
        apply_earth_rotation (float) :  specify (with a number between 0 and 1)
            how much of the earth rotaion is taken in to account in the
            TEC interpolation step.
    '''
    pass
    
def _getRM_azel(time,azels,position,tecinfo, apply_earth_rotation=0):
    '''Get interpolated vTEC, piercepoint ocordinates and parallel B-Field for an array of azimuth elevation pairs, single station, single time
    
    optionally correcting for earth rotation.
    Args:
        time (float) : times in MJD seconds.
        azels (np.array) : az,el pairs in radians.
        position (np.array) : antenna xyz position (ITRF) in meters.
        tecinfo Tuple[np.array, np.array, np.array, np.array, np.array] :
            output of getIonex.read_tec
        apply_earth_rotation (float) :  specify (with a number between 0 and 1)
            how much of the earth rotaion is taken in to account in the
            TEC interpolation step.
    '''
    pass
