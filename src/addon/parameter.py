#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    06.06.2013 20:49:45 EDT
# File:    parameter.py

import sys
from .color import *
from .helper import *
import os
import re
import subprocess

#------------------- parameter class ------------------- 
class parameter_class(dict):
    def __init__(self):
        super(dict, self).__init__()
        self.reserved_names = ["arg", "flag"]
        
        #if a key has an "_" at the end it is treated as "hidden" in the sense that it isn't printed
        self["print_"] = False;
    
    def __str__(self):
        out = ""
        for key in sorted(self.keys()):
            if not self["print_"]:
                if key[-1] == "_":
                    continue
            if key != sorted(self.keys())[0]: #prevent newline at the start
                 out += "\n"
            out += GREENB_ + str(key) + ":\t" + GREEN_ + str(self[key]) + NONE_
            
        return out
    
    def read(self, argv):
        pas = False
        
        #------------------- regex for = notation ------------------- 
        # "p=1", "p= 1", "p =1" and "p = 1" is valid notation, here I transform all to the form "p=1"
        # "p= " and "p =" is invalid notation 
        
        text = " ".join(argv)
        text = re.sub(" ?= ?", "=", text)
        
        argv = text.split(" ")
        
        #------------------- some nice errors ------------------- 
        for w in argv:
            if w[0] == "-" and w.find("=") != -1:
                error("flags cannot be assigned with a equal sigh")
            if w.find("=") != -1 and w.find("=") != w.rfind("="):
                error("too many = signs, check syntax")
            if w[-1] == "=":
                error("no assignment after equal sign")
            if w[0] == "-" and len(w) == 1:
                error("syntax not correct, check -")
                
        
        #------------------- normal set doesn't work bc of reserved_names check ------------------- 
        dict.__setitem__(self, "arg", [])
        dict.__setitem__(self, "flag", [])
        
        for i in range(1, len(argv)):
            w = argv[i]
            #checking if = sign
            if w.find("=") != -1:
                key, val = w.split("=")
                if self.param_set(key):
                        warning("parameter {0} already set to {1} -> overwrite to {2}".format(key, self[key], val))
                self[key] = to_number(val)
                continue
                
            if w[0] == '-' and len(w) > 1:
                if i + 1 < len(argv) and argv[i+1][0] != '-' and argv[i+1].find("=") == -1: #parameter
                    
                    #------------------- just checking for false input ------------------- 
                    if self.param_set(w[1:]):
                        warning("parameter {0} already set to {1} -> overwrite to {2}".format(w[1:], self[w[1:]], argv[i+1]))
                    #------------------- setting the parameter ------------------- 
                    self[w[1:]] = to_number(argv[i+1])
                    pas = True
                else: #flag
                    #------------------- just checking for false input ------------------- 
                    if self.flag_set(w[1:]):
                        warning("flag {0} was already set".format(w[1:]))
                    else:
                        #------------------- setting the flag ------------------- 
                        self["flag"].append(w[1:])
            else:
                if pas:
                    pas = False
                else: #arg
                    #------------------- just checking for false input ------------------- 
                    if self.arg_set(w):
                        warning("arg {0} was already set".format(w))
                    else:
                        #------------------- adding the arg ------------------- 
                        self["arg"].append(w)
                    
    def arg_set(self, arg):
        return arg in self["arg"]
        
    def flag_set(self, flag):
        return flag in self["flag"]
    
    def param_set(self, param):
        return param in self.keys()
    
    def contains(self, key):
        return self.arg_set(key) or self.flag_set(key) or self.param_set(key)
    
    def __setitem__(self, key, val):
        if key in self.reserved_names: #guard reserved names
            error("do not use the following names: {0}".format(self.reserved_names))
        else:
            dict.__setitem__(self, key, val)
    
parameter = parameter_class()

#------------------- parameter action ------------------- 
def bash_if(flag, action):
    if parameter.contains(flag):
        if type(action) == type(" "): #normal bash cmd
            bash(action)
        else: #fct call
            CYAN("called function: ", "")
            action()
    return 0

def bash(cmd):
    CYAN(cmd)
    os.system(cmd)

def popen(cmd):
    CYAN(cmd)
    part = cmd.split(" ")
    return subprocess.check_output([part[0], " ".join(part[1:])])
