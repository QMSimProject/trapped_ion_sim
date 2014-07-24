#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    21.06.2014 15:37:17 CEST
# File:    test_integration.py

from package_proxy import *

if __name__ == "__main__":
    cf = {}
    
    # setting the planck constant
    cf["hbar"] = 1
    
    # setting the atoms [name, energy-levels, starting density matrix]
    cf["atom"] = []
    cf["atom"].append(["atom1", [0, 100], np.diag([1, 0])])
    
    # setting up the vibrons [name, frequence, dimension, starting DM]
    cf["vibron"] = []
    cf["vibron"].append(["vib1", 10, 3, np.diag([1, 0, 0])])
    
    # setting up the lasers [name, freq, phase, firing pattern]
    cf["laser"] = []
    cf["laser"].append(["laser1", 110, 0, [[0.5*np.pi, 1.5*np.pi]]])
    
    # setting the rabi freq: list[#lasers][#atoms][#transitions per atom]
    cf["rabi"] = [[[1 for k in range(len(j[1]) *(len(j[1])-1) / 2)] for j in cf["atom"]] for i2 in cf["laser"]]
    
    # setting the lamb-dicke parameter: list[#atoms][#vibrons]
    cf["eta"] = [[1 for i in cf["vibron"]] for j in cf["atom"]]
    
    # setting integration details
    cf["lower"] = 0
    cf["upper"] = 2 * np.pi
    cf["measure"] = 1000
    cf["cutoff"] = 5
    
    # setting plot options
    cf["plot"] = [[["atom1(fock_dm(N, 0))", "GS1"], ["atom1(fock_dm(N, 1))", "ES1"]], [["vib1(fock_dm(N, 0))", "V0"], ["vib1(fock_dm(N, 1))", "V1"], ["vib1(fock_dm(N, 2))", "V2"]], [["laser1(.5)", "laser1"]]]
    
    # setting collapse ops
    cf["collapse"] = []
    
    t, e = integrate_cf(cf)
    plot(t, e, cf)
