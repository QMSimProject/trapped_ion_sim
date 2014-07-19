#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    12.06.2013 10:31:13 EDT
# File:    helper.py

import os
from .color import *

#------------------- helper ------------------- 
def readable(name):
    return os.access(name, os.R_OK)

#------------------- zip with index  ------------------- 
def zipi(l):
    return zip(l, range(len(l)))

#------------------- converter ------------------- 
def to_number(string):
    try:
        res = int(string)
        return res
    except:
        pass
    try:
        res = float(string)
        return res
    except:
        return string


#------------------- type checks ------------------- 
def is_list(obj):
    return isinstance(obj, list)
    
def is_int(obj):
    return isinstance(obj, int)

def is_float(obj):
    return isinstance(obj, float)

def is_number(obj):
    return is_int(obj) or is_float(obj)

#------------------- depth of a list ------------------- 
def depth(l):
    if is_list(l):
        subdepth = [depth(item) for item in l]
        if subdepth == []:
            return 1
        else:
            return 1 + max(subdepth)
    else:
        return 0

#------------------- ranges ------------------- 
def drange(start, end, step):
    return [step * i for i in range(int(start/step), int(end/step))]

#------------------- debug / error ------------------- 
def ERROR(text):
    raise Exception("{0}error: {1}{2}{3}".format(REDB_, RED_, text, NONE_))

def WARNING(text):
    print("{0}warning: {1}{2}{3}".format(YELLOWB_, YELLOW_, text, NONE_))
