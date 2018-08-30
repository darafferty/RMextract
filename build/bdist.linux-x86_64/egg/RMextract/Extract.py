#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 12:14:31 2018

@author: mevius
"""

import RMextract.getIONEX as gt

def extract(etype="RM", MS=None, server = "ftp://ftp.aiub.unibe.ch/CODE/", prefix="CODG",*args,**kwargs):
    '''Extract RM or TEC
    
    Use MS to extract pointing, station positions and times, or specify with Keyword arguments.
    Keyword arguments will overrule metadata found in MS 
    
    Args: 
        etype (string) : TEC or RM
        MS (string) : (optional) path to the MS file
    
    Returns:
        dict : dictionary with all results
    '''
    times, lonlatpp = _bookkeeping(*args,**kwargs)
    ionexdata = gt.getTECdata(times,lontlatpp,server,prefix)

    

def _bookkeeping(MS=None, 
                 _timerange = 0,
                 _start_time = 0,
                 _end_time = 0,
                 _time_step = 0,
                 stat_pos = [PosTools.posCS002],
                  **kwargs):
    '''Find the timerange and pointing given input parameters'''
    
    if MS is not None:
      timerange,timestep,pointing,stat_names,stat_pos=PosTools.getMSinfo(MS)
    if _timerange != 0:
        start_time = _timerange[0]
        end_time = _timerange[1]
        reference_time = _start_time
        timerange[0] = _start_time - _timestep
        timerange[1] = _end_time + _timestep
    else:
        timerange,reference_time=PosTools.get_time_range(_start_time,_end_time,_timestep,True,0)
    if not stat_names:
        stat_names =['st%d'%(i+1) for i in range(len(stat_pos))]
    