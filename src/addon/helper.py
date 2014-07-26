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
    """
    Checks if a file called name is readable. Uses os.
    """
    return os.access(name, os.R_OK)

#------------------- zip with index  ------------------- 
def zipi(l):
    """
    Shorthand for zip(list, range(len(list))), if one needs the index and the content of a list.
    """
    return zip(l, range(len(l)))

#------------------- converter ------------------- 
def to_number(string):
    """
    Tries to convert the input string into an int, if that doesn't work into a float and if that also fails, returns the string again.
    """
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
    """
    Checks if the obj is a list.
    """
    return isinstance(obj, list)
    
def is_int(obj):
    """
    Checks if the obj is an int.
    """
    return isinstance(obj, int)

def is_float(obj):
    """
    Checks if the obj is a float.
    """
    return isinstance(obj, float)

def is_number(obj):
    """
    Checks if the obj is an int or a float.
    """
    return is_int(obj) or is_float(obj)

#------------------- depth of a list ------------------- 
def depth(l):
    """
    Retruns the maximal depth of a nested list system. Is recursive and searches the whole "tree", might be slow.
    """
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
    """
    Returns a range of floating point numbers.
    """
    return [step * i for i in range(int(start/step), int(end/step))]

#------------------- debug / error ------------------- 
def ERROR(text):
    """
    Raises an exception and outputs text in red color.
    """
    raise Exception("{0}error: {1}{2}{3}".format(REDB_, RED_, text, NONE_))

def WARNING(text):
    """
    Just prints the text in yellow as a warning.
    """
    print("{0}warning: {1}{2}{3}".format(YELLOWB_, YELLOW_, text, NONE_))

def ASSERT(cond, text = ""):
    """
    Will output the text in red color before calling pythons assert if the condition is not True.
    """
    if not cond:
        print("{0}assert failed: {1}{2}{3}".format(REDB_, RED_, text, NONE_))
        assert(cond)
