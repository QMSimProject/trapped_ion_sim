#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    19.07.2014 12:27:41 CEST
# File:    atom_class.py

from .src_import import *

class atom_holder:
    """
    global class that keeps track of all created atoms.
    a global instance is created in this file.
    
    :ivar atom_list: contains all atoms. the atoms add themselves automatic to `atom_holder_` when constructed.
    :ivar current_dims: contains the the dimensions of the atoms.
    """
    def __init__(self):
        """
        the constructor just calls :meth:`clear`
        
        :returns: None
        """
        self.clear()
    
    def trans_correction(self, atom):
        """
        the atoms don't know how many other atoms there are and cannot provide a correct transition matrix that acts on the global state. this  function fixes this by applying the correct tensor-product with the identities of all other atoms.
        
        :param atom: an atom of type :class:`atom`
        :returns: the corrected :attr:`atom.trans` list
        """
        #apply all identities to the new trans list
        for d in reversed(self.current_dims):
            for t in atom.trans:
                t[1] = q.tensor(q.qeye(d), t[1])
        
        #correct all current atoms by adding (tensor) the identity of the new atom
        for at in self.atom_list:
            for i in range(len(at.trans)):
                at.trans[i][1] = q.tensor(at.trans[i][1], q.qeye(atom.N))
        
        #update internal
        self.current_dims.append(atom.N)
        self.atom_list.append(atom)
        
        return atom.trans
    
    def state(self):
        """
        this function returns the global state of all particles. it takes all single states and tensors them together.
        
        :returns: a `qutip.Qobj` containing the density matrix of the collective atom state
        """
        s = self.atom_list[0].state
        for a in self.atom_list[1:]:
            s = q.tensor(s, a.state)
        return s
    
    def clear(self):
        """
        clears the internal list of atoms and dimensions.
        
        :returns: None
        """
        self.atom_list = []
        self.current_dims = []

atom_holder_ = atom_holder()
"""
the only instance of atom_holder that should exist.
"""


class atom:
    """
    this class stores information for an atom.
    
    :ivar energy_levels: a list of increasing energy levels
    :ivar state: the density matrix of the initial internal state of the atom
    :ivar N: the dimension of the atom-space
    :ivar trans: a list of pairs for all transitions possible. each pair contains [transition frequency, transition matrix], where the transition matrix is 0 except for one entry that is 1. There are :attr:`N` * (:attr:`N` - 1) / 2 transitions.
    """
    global atom_holder_
    
    def __init__(self, levels, state = "none"):
        """
        the atom-constructor
        
        :param levels: a list of increasing energy levels
        :param state: optional parameter of the initial density matrix. needs to be convertable to `qutip.Qobj`. if no state is given, it will be set to the pure groundstate.
        
        :returns: None
        """
        self.energy_levels = levels
        self.N  = len(self.energy_levels)
        
        if state == "none":
            state = q.fock_dm(self.N, 0)
        self.state = q.Qobj(state)
        
        self.trans = []
        
        #create all transitions
        for i in range(self.N):
            for j in range(i + 1, self.N):
                trans_energy = np.abs(self.energy_levels[i] - self.energy_levels[j])
                m = np.zeros((self.N, self.N))
                m[i, j] = 1
                self.trans.append([trans_energy, q.Qobj(m)])
        #and correct the matrix with the propper spaces of other atoms
        #register the atom in the global atom_holder_
        self.trans = atom_holder_.trans_correction(self)
        
    def __str__(self):
        """
        used by the print() function for nice output.
        
        :returns: a string representation of the atom
        """
        r = ""
        r += "atom\n"
        r += str(self.N) + " enegry levels: " + str(self.energy_levels) + "\n"
        r += "current state: "
        r += str(self.state)
        return r
    
    def transition(self, nr):
        """
        get function
        
        :param nr: the position of the required pair
        :returns: the pair nr of the list :attr:`trans`
        """
        return self.trans[nr]
