#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    19.07.2014 12:32:03 CEST
# File:    laser_class.py

from .src_import import *

class laser:
    """
    this class stores information for a laser.
    
    :ivar pattern: the firing pattern of the laser i.e. intervalls when the laser in on.
    :ivar freq: the frequency of the laser
    :ivar phase: optional. the phase of the laser. is 0 if no value is given.
    """
    def __init__(self, pattern, freq, phase = 0):
        """
        constructor of the laser
        
        :param pattern: the firing pattern of the laser. needs to have the form [[0, 1], [2, 3], ...].
        :param freq: the frequency of the laser
        :param phase: the phase of the laser
        :returns: None
        """
        self.pattern = pattern
        self.freq = freq
        self.phase = phase
