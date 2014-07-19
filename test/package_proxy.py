#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    01.05.2014 19:44:38 CEST
# File:    package_proxy.py

# adds the src path in order to import it
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../src")

# importing the trapped_ion_sim library
from src import *
import addon

# just the usual imports
import qutip as q
import numpy as np
import pylab as pl
