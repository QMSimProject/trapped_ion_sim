#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    19.07.2014 12:30:05 CEST
# File:    vibron_class.py

from .src_import import *

class vibron_holder:
    """
    global class that keeps track of all created vibrons.
    a global instance is created in this file.
    
    :ivar vibron_list: contains all vibrons. the vibrons add themselves automatic to `vibrons_holder_` when constructed.
    :ivar current_dims: contains the the dimensions of the vibrons.
    :ivar N: contains the total dimensions of the vibron space (i.e. multiply all entries of :attr:`current_dims`).
    """
    def __init__(self):
        """
        the constructor just calls :meth:`clear`
        
        :returns: None
        """
        self.clear()
    
    def destroy_correction(self, vibron):
        """
        the vibrons don't know how many other atoms there are and cannot provide a correct annihilator that acts on the global state. this  function fixes this by applying the correct tensor-product with the identities of all other vibrons.
        
        :param vibron: an vibron of type :class:`vibron`
        :returns: the corrected :attr:`vibron.destroy`
        """
        #apply all identities to the new annihilator
        for d in reversed(self.current_dims):
            vibron.destroy = q.tensor(q.qeye(d), vibron.destroy)
            
        #correct all current vibrons by adding (tensor) the identity of the new vibron
        for vb in self.vibron_list:
            vb.destroy = q.tensor(vb.destroy, q.qeye(vibron.N))
        
        #update internal
        self.current_dims.append(vibron.N)
        self.N *= vibron.N
        self.vibron_list.append(vibron)
        return vibron.destroy
    
    def state(self):
        """
        this function returns the global state of all vibrons. it takes all single states and tensors them together.
        
        :returns: a `qutip.Qobj` containing the density matrix of the collective vibron state
        """
        if len(self.vibron_list) == 0:
            return 1
        s = self.vibron_list[0].state
        for v in self.vibron_list[1:]:
            s = q.tensor(s, v.state)
        return s
    
    def clear(self):
        """
        clears the internal list of atoms and dimensions.
        
        :returns: None
        """
        self.vibron_list = []
        self.current_dims = []
        self.N = 1

vibron_holder_ = vibron_holder()
"""
the only instance of atom_holder that should exist.
"""

class vibron:
    """
    this class stores information for a vibron.
    
    :ivar freq: the frequency of the vibrational mode
    :ivar state: the density matrix of the initial internal state of the vibron
    :ivar N: the dimension of the vibron-space
    :ivar destroy: the annihilator operator of the vibron
    """
    global vibron_holder_
    
    def __init__(self, freq, N, state = "none"):
        """
        the vibron-constructor
        
        :param freq: the frequency of the vibrational mode
        :param N: the dimension of the vibron (if you work with thermal states always check that this is large enough)
        :param state: optional parameter of the initial density matrix. needs to be convertable to `qutip.Qobj`. if no state is given, it will be set to the pure groundstate.
        
        :returns: None
        """
        self.freq = freq
        self.N = N
        if state == "none":
            state = q.fock_dm(self.N, 0)
        
        self.state = q.Qobj(state)
        self.destroy = q.destroy(self.N)
        
        #correct and register the vibron according to the other vibrons
        self.destroy = vibron_holder_.destroy_correction(self)
