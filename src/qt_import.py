#!/usr/bin/python3.2
# -*- coding: cp1252 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    16.01.2012 17:35:43 CET
# File:    qt_import.py

from .addon import *

# this tries to import PySide ot PyQt4 
# (they are almost always interchangeable)
# both are python binding to the c++ QT library
# I slightly prefer pyside but PyQt4 is on the same level today
qt_binding = "none"

try:
    GREEN("PySide loaded")
    from PySide.QtCore import *
    from PySide.QtGui import *
    qt_binding = "PySide"
except ImportError:
    try:
        from PyQt4.QtCore import *
        from PyQt4.QtGui import *
        qt_binding = "PyQt4"
        GREEN("PyQt4 loaded")
    except ImportError:
        ERROR("no PySide or PyQt4 module found")
