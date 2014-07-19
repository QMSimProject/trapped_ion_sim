#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    19.07.2014 11:57:08 CEST
# File:    01_basic_integration.py

from qutip import *
from pylab import *

def run_example():
    
    # we have a spin 1/2 system that is in the up state
    # psi can either be a a vector or a density matrix
    psi = fock_dm(2, 0) #dimensions / state

    # a time independent Hamiltonian can just be given a an operator
    H = [sigmay()]
    
    # the observables that are mesured during the integration
    obs = [sigmax(), sigmay(), sigmaz()]
    
    # collaps operators if needed
    col = [0.3 * sigmax()]
    
    # the integration time (0 to 10) and the amount of measurement points (for the observables)
    tlist = linspace(0, 3*pi, 201)
    
    # the integration (uses ode45)
    res = mesolve(H, psi, tlist, col, obs)
    
    # get the measured values of the observables
    exp = res.expect
    
    # plot observables
    for ex in exp:
        plot(res.times, ex)
        
    xlabel("Time")
    ylabel("exp val")
    legend(["x", "y", "z"])
    show()
    
    # plot bloch sphere
    bloch = Bloch()
    bloch.add_points(exp[0:3]) # [0:3] are x, y, z
    bloch.make_sphere()
    show()
    
if __name__ == "__main__":
    print("01_basic_integration.py")
    run_example()
