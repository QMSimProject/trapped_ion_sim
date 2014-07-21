#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    21.06.2014 15:37:17 CEST
# File:    test_integration.py

from package_proxy import *

if __name__ == "__main__":
    cf = {}
    
    cf["atom"] = []
    cf["atom"].append(["Ca", [0, 100], np.diag([1, 0])])
    cf["atom"].append(["atom2", [0, 100], np.diag([1, 0])])
    
    cf["vibron"] = []
    cf["vibron"].append(["vib1", 10, 3, np.diag([.3, .7, 0])])
    
    cf["laser"] = []
    cf["laser"].append(["laser1", 110+1, 0, [[0*np.pi, 100*np.pi]]])
    cf["laser"].append(["laser2", 90-1,  0, [[0*np.pi, 100*np.pi]]])
    
    cf["rabi"] = [[[0.2 for k in range(len(j[1]) *(len(j[1])-1) / 2)] for j in cf["atom"]] for i2 in cf["laser"]]
    cf["eta"] = [[0.5 for i in cf["vibron"]] for j in cf["atom"]]
    
    cf["lower"] = 0
    cf["upper"] = 100 * np.pi
    cf["measure"] = 1000
    cf["cutoff"] = 5
    cf["plot"] = [["atom1(dest)", "atom1(fock_dm(2, 1))"], ["atom2(fock_dm(2, 0))", "atom2(fock_dm(2, 1))"], ["vib1(fock_dm(3, 0))", "vib1(fock_dm(3, 1))", "vib1(fock_dm(3, 2))"]]
    
    #~ , ["fidelity(tensor(atom2(fock_dm(2, 0)), atom1(fock_dm(2, 0))), state)"], ["laser1"]
    
    #~ integrate_cf(cf) #old method
    t, e = integrate_cf_eff(cf)
    plot(t, e, cf)
