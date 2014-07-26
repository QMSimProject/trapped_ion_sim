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
    """
    An useful class to parse argv. Is derived from dict.
    """
    def __init__(self):
        """
        Initializes the class. Do not used self.reserved_names as keys for the dict. It will raise an error.
        """
        super(dict, self).__init__()
        self.reserved_names = ["arg", "flag"]
        
        #if a key has an "_" at the end it is treated as "hidden" in the sense that it isn't printed
        self["print_"] = False
        self["warn_"] = True
    
    def __str__(self):
        """
        String conversion for printing. If they key "print_" is set to True all keys with an trailing underscore will be printed as well. Good for hiding technical/private keys.
        """
        out = ""
        for key in sorted(self.keys()):
            if not self["print_"]:
                if key[-1] == "_":
                    continue
            if key != sorted(self.keys())[0]: #prevent newline at the start
                 out += "\n"
            out += GREENB_ + str(key) + ":\t" + GREEN_ + str(self[key]) + NONE_
            
        return out
    
    def warn(self, text):
        """
        If the key "warn_" is True, the user will be warned if a key is overwritten or a flag set a second time.
        """
        if self["warn_"]:
            WARNING(text)
    
    def read(self, argv):
        """
        The passed argv is parsed and stored in the right places. Legal notation is:\n 
        - -flag
        - -param value
        - param = value
        - arg
        """
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
                ERROR("flags cannot be assigned with a equal sigh")
            if w.find("=") != -1 and w.find("=") != w.rfind("="):
                ERROR("too many = signs, check syntax")
            if w[-1] == "=":
                ERROR("no assignment after equal sign")
            if w[0] == "-" and len(w) == 1:
                ERROR("syntax not correct, check -")
                
        
        #------------------- normal set doesn't work bc of reserved_names check ------------------- 
        dict.__setitem__(self, "arg", [])
        dict.__setitem__(self, "flag", [])
        
        for i in range(1, len(argv)):
            w = argv[i]
            #checking if = sign
            if w.find("=") != -1:
                key, val = w.split("=")
                if self.has_param(key):
                    self.warn("parameter {0} already set to {1} -> overwrite to {2}".format(key, self[key], val))
                self[key] = to_number(val)
                continue
                
            if w[0] == '-' and len(w) > 1:
                if i + 1 < len(argv) and argv[i+1][0] != '-' and argv[i+1].find("=") == -1: #parameter
                    
                    #------------------- just checking for false input ------------------- 
                    if self.has_param(w[1:]):
                        self.warn("parameter {0} already set to {1} -> overwrite to {2}".format(w[1:], self[w[1:]], argv[i+1]))
                    #------------------- setting the parameter ------------------- 
                    self[w[1:]] = to_number(argv[i+1])
                    pas = True
                else: #flag
                    #------------------- just checking for false input ------------------- 
                    if self.has_flag(w[1:]):
                        self.warn("flag {0} was already set".format(w[1:]))
                    else:
                        #------------------- setting the flag ------------------- 
                        self["flag"].append(w[1:])
            else:
                if pas:
                    pas = False
                else: #arg
                    #------------------- just checking for false input ------------------- 
                    if self.has_arg(w):
                        self.warn("arg {0} was already set".format(w))
                    else:
                        #------------------- adding the arg ------------------- 
                        self["arg"].append(w)
                    
    def has_arg(self, arg):
        """
        Checks if arg is in the parameter_class. An arg is an entry without a value, like a filename.
        """
        return arg in self["arg"]
        
    def has_flag(self, flag):
        """
        Checks if flag is set. A flag does not have a value. It is eighter on (has_flag -> True) or off (has_flag -> False)
        """
        return flag in self["flag"]
    
    def has_param(self, param):
        """
        Checks if the key param is set. A param is a key with a value.
        """
        return param in self.keys()
    
    def has_key(self, key):
        """
        Checks if key is a parameter, flag or arg. It does not the same as the dict.has_key function, since flags and args aren't technically stored as keys.
        """
        return self.has_arg(key) or self.has_flag(key) or self.param_set(key)
    
    def __setitem__(self, key, val):
        """
        Forwards to the dict.__setitem__ but makes sure that the reserved_names aren't used
        """
        if key in self.reserved_names: #guard reserved names
            ERROR("do not use the following names: {0}".format(self.reserved_names))
        else:
            dict.__setitem__(self, key, val)
    
parameter = parameter_class()

#------------------- parameter action ------------------- 
def bash_if(flag, action):
    """
    If the flag is in the parameter instance of parameter_class, the action will be executed by as a bash command if it is a string and be called otherwise (assumption action == python function with no args)
    """
    if parameter.contains(flag):
        if type(action) == type(" "): #normal bash cmd
            bash(action)
        else: #fct call
            CYAN("called function: ", "")
            action()
    return 0

def bash(cmd):
    """
    Just calls os.system and outputs the command.
    """
    CYAN(cmd)
    os.system(cmd)

def popen(cmd):
    """
    If one needs the output of the bash-command, this function can provide it. Works exactly like bash(cmd) but returns the output as a string.
    """
    CYAN(cmd)
    part = cmd.split(" ")
    return subprocess.check_output([part[0], " ".join(part[1:])])
