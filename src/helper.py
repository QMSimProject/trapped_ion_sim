#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    23.07.2014 23:35:05 CEST
# File:    helper.py

from .src_import import *
from .atom_class import *
from .vibron_class import *
from .useable_fct import *

#=================== helper ===================
def hold_args(names, expr):
    """
    helper function for :func:`create_obs` and :func:`create_collpase_ops`
    
    the function converts the pattern "name(....)" to "name('....')" if name is in the list names.
    
    :param names: a list of names who's argument should be stringified
    :param expr: and exprssion in string form
    :returns: the modified expression
    """
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
    """
    this class is used by :func:`create_obs` and :func:`create_collpase_ops`
    
    :ivar idx: the index of an atom or vibron inside the atom-holder or vibron holder
    :ivar offset: is 0 for atoms and len(atoms) for vibrons
    :ivar N: is the dimension of the system with index idx in the global state
    :ivar all_sys: is a static list that contains all atoms and vibrons. needs to be set before using tag_class.
    """
    all_sys = []
    def __init__ (self, name, idx, offset = 0):
        """
        the constructor just sets `idx`, `offset` and `N`
        
        
        :param name: the name of the atom/vibron (not used yet, but perhaps useful in the future)
        :param idx: the index of the syste,
        :param offset: is set to zero if not specified
        :returns: None
        """
        self.idx = idx
        self.off = offset
        self.N = self.__class__.all_sys[self.idx + self.off].N
    def __int__(self):
        """
        an int cast
        
        :returns: `idx` + `offset`
        """
        return self.idx + self.off
    def __call__(self, op):
        """
        takes an operator (Qobj) but in string form. The operator can contain N as a variable that this function will provide.
        
        :returns: a list with [`idx` + `offset`, evaluated operator (Qobj)]
        """
        if op == '':
            op = "0"
        
        N = self.N
        return [self.idx + self.off, eval(op)]

    
def tensor(*ops):
    """
    a tensor function. if the input is a Qobj nothing happens, thus tensor(tensor(x)) == tensor(x).
    
    :param ops: a list of return values of `tag_class.__call__` or a Qobj
    :returns: an operator that acts on the full state containing the inserted `ops`
    """
    if isinstance(ops[0], q.Qobj): #tensor(tensor(x)) == tensor(x)
        return ops[0]
    sys_list = []
    for d, i_d in zipi(atom_holder_.current_dims + vibron_holder_.current_dims):
        sys_list.append(q.qeye(d))
    for op in ops:
        sys_list[op[0]] = op[1]
    
    return q.tensor(sys_list)
