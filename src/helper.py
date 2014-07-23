#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    23.07.2014 23:35:05 CEST
# File:    helper.py

from .src_import import *
from .atom_class import *
from .vibron_class import *

#=================== helper ===================
def hold_args(names, expr):
    for n in names:
        pos = 0
        pos = expr.find(n)
        
        def find_matching_bracket(expr, pos_open):
            assert(expr[pos_open] == "(")
            pos = pos_open + 1
            
            opened = 1
            while True:
                next_open = expr.find("(", pos)
                next_close = expr.find(")", pos)
                
                if next_open == -1 or next_close < next_open:
                    opened -= 1
                    if opened == 0:
                        pos_closed = next_close
                        break
                    else:
                        pos = next_close + 1
                else:
                    opened += 1
                    pos = next_open + 1
            
            return expr[:pos_open + 1], expr[pos_open + 1: pos_closed], expr[pos_closed:]
        
        while pos != -1:
            ipos = pos + len(n)
            f = find_matching_bracket(expr, ipos)
            expr = f[0] + "'" + f[1] + "'" + f[2]
            pos = expr.find(n, pos + 2)
    
    return expr

class tag_class():
    all_sys = []
    def __init__ (self, name, idx, offset = 0):
        self.name = name
        self.idx = idx
        self.off = offset
        self.N = self.__class__.all_sys[self.idx + self.off].N
    def __int__(self):
        return self.idx + self.off
    def __call__(self, op):
        from qutip import fock_dm
        from qutip import coherent_dm
        from qutip import thermal_dm
        from qutip import squeeze
        from qutip import destroy
        from qutip import create
        from qutip import sigmax
        from qutip import sigmay
        from qutip import sigmaz
        N = self.N
        return [self.idx + self.off, eval(op)]

    
def tensor(*ops):
    if isinstance(ops[0], q.Qobj): #tensor(tensor(x)) == tensor(x)
        return ops[0]
    sys_list = []
    for d, i_d in zipi(atom_holder_.current_dims + vibron_holder_.current_dims):
        sys_list.append(q.qeye(d))
    for op in ops:
        sys_list[op[0]] = op[1]
    
    return q.tensor(sys_list)
