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
import logging

def getRM(*args,**kwargs):
    
    '''Dummy function for compatibility'''
    logging.warning("getRM is depreciated, please use extract(RM) instead")
    RMextract.extract("RM",*args,**kwargs)
    

