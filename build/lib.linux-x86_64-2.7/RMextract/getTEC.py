#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 11:54:23 2018

@author: mevius
"""
import logging
import RMextract.extract

def getTEC(*args, **kwargs):
    '''dummy function for compatibility'''
    logging.warning("getTEC is depreciated, please use extract(\"TEC\") \
    instead")
    RMextract.extract("TEC", *args, **kwargs)