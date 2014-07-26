#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    04.03.2014 11:44:14 CET
# File:    optical.py

from .helper import *

#------------------- progress ------------------- 
def progress_bar(p):
    """
    p is a number between 0 and 1 and progress_bar will return a string that contains a progress bar representing that progress.
    """
    size = 50
    
    bars = int(min(1, p) * size);
    bar = "<"

    for i in range(bars):
        bar += "|"
    for i in range(bars, size):
        bar += " "
    bar += ">"
    end = "{0:3}%".format(int(p*100))
    
    if(p < .33):
        return RED_ + bar + NONE_ + end
    elif(p < .66):
        return YELLOW_ + bar + NONE_ + end
    elif(p < 1):
        return GREEN_ + bar + NONE_ + end
    elif(p == 1):
        return CYAN_ + bar + NONE_ + end
    else:
        warning("progress > 100%")
        return CYAN_ + bar + NONE_ + end
        
