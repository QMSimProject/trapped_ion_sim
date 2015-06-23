#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    21.06.2014 23:49:11 CEST
# File:    all_to_all.py

from .src_import import *
from .atom_class import *
from .vibron_class import *
from .laser_class import *
from .helper import hold_args, tag_class, tensor

#=================== global variables ===================
default_cutoff = 100000
"""
the default cutoff frequency used by :func:`create_hamiltonian`
"""

def create_phi():
    """
    creates the initial global state phi by adding (tensor) the atom and vibron states
    
    :returns: a `qutip.Qobj` that contains the density matrix of the global state
    """
    if len(vibron_holder_.current_dims) == 0:
        return atom_holder_.state()
    return q.tensor(atom_holder_.state(), vibron_holder_.state()).unit()

class res_collector:
    """
    the mesovle as well as other solve methods of qutip take an list of Qobj matrices as observables or a callback function. this class contains the callback functions as well as the datastructure that collects the data. it is done in the same fashion as in mesolve.py of the qutip framework. the callback is only used if there are measurements that aren't represented by the form expect(operator, state).
    """
    def __init__(self):
        """
        the constructor initializes the class with one timestep and no functions
        
        :returns: None
        """
        self.init([], 1)
    def init(self, fct, n_tsteps):
        """
        init is used by :func:`create_obs` if there are non expectation value (see class description) measurement functions.
        
        :param fct: is a list of python functions with the signature: f(time, state) and return a float
        :param n_tsteps: has to be the same as cf["measure"] (represents how many measurements are taken during the integration)
        
        :returns: None
        """
        self.t_idx = 0
        self.times = np.zeros(n_tsteps)
        self.expect = []
        self.fct = fct
        for n in range(len(self.fct)):
            self.expect.append(np.zeros(n_tsteps))
    def measure(self, t, state):
        """
        this function will call all functions in `fct` and store the return value.
        
        :param t: the time at which the measurement is taken
        :param state: the density matrix representing the state at time t
        
        :returns: None
        """
        i = self.t_idx
        self.times[i] = t
        for f, i_f in zipi(self.fct):
            self.expect[i_f][i] = f(t, state)
        
        self.t_idx += 1

res = res_collector()
"""
the only instance that should be used in the framework.
"""


def callback_fct(t, state):
    """
    this is the callback function given to the solver routines. it calls `res.measure`
    
    :param t: the time at which the measurement is taken
    :param state: the density matrix representing the state at time t
    
    :returns: None
    """
    global res
    res.measure(t, q.Qobj(state))
    
def create_obs(cf):
    """
    creates the observables specified in cf["plot"]. one can use a simplified language to formulate the observables in an more convenient way. if a certain form is not supported it can be expanded to do so. For now valid discription code is:
    
    - "fidelity(atom1(...), state); label" important is the "; " before writing the label (how the plot line is called)
    - "atom1(...); label" this will produce an operator. the expectation value will be measured 
    - "tensor(atom1(...), atom2(...)); label"
    - "laser1(.3); label" the argument in laser1 just determines the y-location in the plot. it has no computational meaning
    
    :param cf: canonical form of the integration problem
    :return: obs, laser_obs, callback
    - `obs`: is eighter a list of operators or a python function (callback, see QuTip mesolve)
    - `laser_obs`: a list of functions that are calculated after the integration, since the lasers are independent of the state
    - `callback`: a boolean value that states if the callback or the obs_list is returned. this is done for speedup, since the callback is slower and only used if needed.
    """
    global res
    flatten = lambda l: [y for x in l for y in x]
    
    tag_class.all_sys = atom_holder_.atom_list + vibron_holder_.vibron_list

    pi = flatten(cf["plot"])
    #------------------- collect all names ------------------- 
    a_nm = [a[0] for a in cf["atom"]]
    v_nm = [a[0] for a in cf["vibron"]]
    l_nm = [a[0] for a in cf["laser"]]
    
    pi = [hold_args(a_nm + v_nm, x[0]) for x in pi]
    
    all_nm = a_nm + v_nm + l_nm
    
    for a, i_a in zipi(a_nm):
        exec(a + " =  tag_class('" + a + "', " + str(i_a) + ")") in locals(), globals()
    
    for v, i_v in zipi(v_nm):
        exec(v + " =  tag_class('" + v + "', " + str(i_v) + ", " + str(len(a_nm)) + ")") in locals(), globals()
    
    class laser_helper():
        def __init__(self, pattern):
            self.pattern = pattern
        def __call__(self, m):
            self.magn = m
            def closure(t):
                for l in self.pattern:
                    if t < l[1] and t >= l[0]:
                        return self.magn
                return 0
            return closure
    
    for l, i_l in zipi(l_nm):
        exec(l + " = laser_helper(" + str(cf["laser"][i_l][3]) + ")") in locals()
    
    def fidelity(ref, s):
        op = ref[1]
        sys = [ref[0]]
        
        def closure(t, state):
            return q.fidelity(op, state.ptrace(sys))
        return closure
    
    from qutip import expect
    
    callback = False
    fct = []
    op_map = []
    laser_fct = []
    
    for p, i_p in zipi(pi):
        is_laser = False
        for l in l_nm:
            if p.find(l) != -1:
                is_laser = True
                break
        
        if is_laser:
            exec("f = " + p) in locals(), globals()
            laser_fct.append([f, i_p])
        elif p.find("fidelity") != -1:
            callback = True
            state = 1
            exec("f = " + p) in locals(), globals()
            fct.append(f)
        else:
            op = tensor(eval(p))
            op_map.append(op)
            def closure(i):
                def f(t, state):
                    return expect(op_map[i], state)
                return f
            
            fct.append(closure(len(op_map) - 1))
    
    if not callback:
        return op_map, laser_fct, callback
    else:
        res.init(fct, cf["measure"])
        return callback_fct, laser_fct, callback

def create_collpase_ops(cf):
    """
    this function creates the collpase operators based on the string inputs on cf["collapse"]. examples:
    
    - "atom1(fock_dm(N, 2))" inside the atom1 function any valid expression that results in a Qobj of the right dimension can be given
    - "0.5*atom1(...)" doesn't work, since atom1(...) is not an operator from the programms perspective
    - "0.5*tensor(atom1(...))" works since tensor transforms the atom1 obj into a Qobj, where 0.5* is valid code
    - "tensor(atom1(...), atom5(...))" the tensor automatically produces the right global state
    - timedependent collapse ops not possible yet, but the function can be expanded
    
    :param cf: canonical form of the integration problem
    :return: a list of operators
    """
    col = cf["collapse"]
    
    a_nm = [a[0] for a in cf["atom"]]
    v_nm = [a[0] for a in cf["vibron"]]
    
    col = [hold_args(a_nm + v_nm, x) for x in col]
    
    for a, i_a in zipi(a_nm):
        exec(a + " =  tag_class('" + a + "', " + str(i_a) + ")") in locals(), globals()
    
    for v, i_v in zipi(v_nm):
        exec(v + " =  tag_class('" + v + "', " + str(i_v) + ", " + str(len(a_nm)) + ")") in locals(), globals()
    
    col_ops = []
    
    for c in col:
        op = tensor(eval(c))
        col_ops.append(op)
    
    return col_ops
    
def create_hamiltonian(atoms, vibrons, lasers, _rabi, _eta, **kwargs):
    """
    :param atoms: a list of :class:`atom`
    :param vibrons: a list of :class:`vibron`
    :param lasers: a list of :class:`laser`
    :param _rabi: a list containing the rabi-frequency with the following dimensions: list[N_lasers][N_atoms][N_transitions_per_atom]. since the amount of transitions can be different for each atom, the last dimension will depend on the atom and is not the same for all atoms.
    :param _eta: a list containing the lamb-dicke parameter with the following dimensions: list[N_atoms][N_vibrons]
    :param **kwargs: key word arguments for additional optional parameters like cutoff. if no cutoff is given, it will be set to `default_cutoff`. hbar is an mandatory kwarg and is the value for the reduced planck constant.
    :returns: a list with [operator, pattern, time dependency] entries for each component in the hamiltonian as well as a dictionary that holds all variables that can be found in the time dependency cython strings.
    """
    global vibron_holder_
    
    if "cutoff" not in kwargs:
        kwargs["cutoff"] = default_cutoff
    cutoff = kwargs["cutoff"]
    
    H = []
    col_args = {}
    col_args["hb"] = kwargs["hbar"]
    
    rwa_takes = [0, 0]
    
    for l, i_l in zipi(lasers):
        for a, i_a in zipi(atoms):
            for t, i_t in zipi(a.trans):
                trans_label = "w0_A" + str(i_a) + "T" + str(i_t)
                rabi_label = "Ohm_L" + str(i_l) + "A" + str(i_a) + "T" + str(i_t)
                laser_label = "wL_L" + str(i_l)
                phase_label = "fL_L" + str(i_l)
                #------------------- set up the operators ------------------- 
                col_args[trans_label] = t[0]
                col_args[rabi_label]  = _rabi[i_l][i_a][i_t]
                col_args[laser_label] = l.freq
                col_args[phase_label] = l.phase
                
                d_internal = t[1]
                ci = q.tensor(d_internal.dag(), q.qeye(vibron_holder_.N))
                di = q.tensor(d_internal      , q.qeye(vibron_holder_.N))
                
                Hc= [ [ci, l.pattern, " 0.5 * hb * {0} * ( exp( 1j * (({1} - {2}) * t + {3})) )".format(rabi_label,     trans_label, laser_label, phase_label), [i_a]]
                    , [di, l.pattern, " 0.5 * hb * {0} * ( exp(-1j * (({1} - {2}) * t + {3})) )".format(rabi_label, trans_label, laser_label, phase_label), [i_a]]
                    ]

                #------------------- rotating wave approx ------------------- 
                rwa_c = [ abs(col_args[trans_label] - col_args[laser_label])
                        , abs(col_args[trans_label] - col_args[laser_label])
                        ]
                #------------------- add carrier terms ------------------- 
                if col_args[rabi_label] == 0: #then the term will be zero anyway
                    continue
                    
                rwa_takes[1] += len(rwa_c)
                for h, diff in zip(Hc, rwa_c):
                    if diff < cutoff:
                        H.append(h)
                        rwa_takes[0] += 1
                    
                
                for v, i_v in zipi(vibrons):
                    eta_label = "eta_A" + str(i_a) + "V" + str(i_v)
                    vib_label = "wz_V" + str(i_v)
                    #------------------- set up the constants ------------------- 
                    col_args[eta_label] = _eta[i_a][i_v]
                    col_args[vib_label] = v.freq
                    d_osz = v.destroy
                    
                    #------------------- set up the operators ------------------- 
                    cc = q.tensor(d_internal.dag(), d_osz.dag())
                    cd = q.tensor(d_internal.dag(), d_osz)
                    dc = q.tensor(d_internal      , d_osz.dag())
                    dd = q.tensor(d_internal      , d_osz)
                    
                    Hs= [ [cc, l.pattern, " 0.5 * hb * {0} * ( exp( 1j * (({1} - {2}) * t + {3})) ) * 1j * {4} * exp( 1j * {5} * t)".format(rabi_label, trans_label, laser_label, phase_label, eta_label, vib_label), [i_a, len(atoms) + i_v]]
                        , [cd, l.pattern, " 0.5 * hb * {0} * ( exp( 1j * (({1} - {2}) * t + {3})) ) * 1j * {4} * exp(-1j * {5} * t)".format(rabi_label, trans_label, laser_label, phase_label, eta_label, vib_label), [i_a, len(atoms) + i_v]]
                        , [dc, l.pattern, "-0.5 * hb * {0} * ( exp(-1j * (({1} - {2}) * t + {3})) ) * 1j * {4} * exp( 1j * {5} * t)".format(rabi_label, trans_label, laser_label, phase_label, eta_label, vib_label), [i_a, len(atoms) + i_v]]
                        , [dd, l.pattern, "-0.5 * hb * {0} * ( exp(-1j * (({1} - {2}) * t + {3})) ) * 1j * {4} * exp(-1j * {5} * t)".format(rabi_label, trans_label, laser_label, phase_label, eta_label, vib_label), [i_a, len(atoms) + i_v]]
                        ]
                    
                    #------------------- rotating wave approx ------------------- 
                    rwa_s = [ abs(col_args[trans_label] - col_args[laser_label] + col_args[vib_label])
                            , abs(col_args[trans_label] - col_args[laser_label] - col_args[vib_label])
                            , abs(col_args[trans_label] - col_args[laser_label] - col_args[vib_label])
                            , abs(col_args[trans_label] - col_args[laser_label] + col_args[vib_label])
                            ]
                    #------------------- add sideband terms ------------------- 
                    if col_args[eta_label] == 0:
                        continue
                        
                    rwa_takes[1] += len(rwa_s)
                    for h, diff in zip(Hs, rwa_s):
                        if diff < cutoff:
                            H.append(h)
                            rwa_takes[0] += 1
    GREEN("{1} of {2} terms survived the RWA (cutoff: {0})".format(cutoff, *rwa_takes))
    return H, col_args

import time
import datetime

def integrate_cf(cf, **kwargs):
    """
    this function implements the event integration for the `cf`. the event integration recalculates the hamiltonian everytime a laser is switched on or off. this avoids to recalculate parts of the hamiltonian, that are zero anyways during a certain time frame.
    
    :param cf: the canonical form of the integration problem
    :return: (times, exp) a list with the times and a list of lists for every observed property exp[#obs][#measure]
    """
    global res
    GREEN("start integration")
    start_time = time.time()
    
    #------------------- reset holders ------------------- 
    atom_holder_.clear()
    vibron_holder_.clear()
    #------------------- checking kwargs ------------------- 
    at = []
    for a in cf["atom"]:
        at.append(atom(a[1], a[2]))
    
    vb = []
    for v in cf["vibron"]:
        vb.append(vibron(v[1], v[2], v[3]))
    
    ls = []
    events = set([0])
    for l in cf["laser"]:
        ls.append(laser(l[3], l[1], l[2]))
        for intervall in l[3]:
            events.add(intervall[0])
            events.add(intervall[1])
    events = sorted(list(events))
        
    rabi = cf["rabi"] #matrix
    eta = cf["eta"] #matrix
    
    H, args = create_hamiltonian(at, vb, ls, rabi, eta, cutoff = cf["cutoff"], hbar = cf["hbar"])
    
    #------------------- set up event structure ------------------- 
    act = set()
    col_event = []
    for e in events:
        for h, i_h in zipi(H):
            for it in h[1]:
                if it[0] == e:
                    act.add(i_h)
                if it[1] == e:
                    act.remove(i_h)
        col_event.append([e, list(act)])
    
    #col_event is a list with the following content:
    #col_event[N][0]: event time
    #col_event[N][1]: systems active (bc of laser pulses)
    
    #------------------- assign tlist parts to col_res ------------------- 
    tlist = np.linspace(cf["lower"], cf["upper"], cf["measure"])
    space = cf["upper"] - cf["lower"]
    
    for i_e in range(len(col_event)):
        if col_event[i_e][0] <= cf["lower"]:
            start = i_e
        if col_event[i_e][0] < cf["upper"]:
            end = i_e
    
    for i_e in range(start, end + 1):
        g = lambda s, e: [s, e, max(2, int((e - s) / space * cf["measure"]))]
        
        if i_e == start:
            if i_e == end:
                col_event[i_e][0] = g(cf["lower"], cf["upper"])
            else:
                col_event[i_e][0] = g(cf["lower"], col_event[i_e + 1][0])
        elif i_e == end:
            col_event[i_e][0] = g(col_event[i_e][0], cf["upper"])
        else:
            col_event[i_e][0] = g(col_event[i_e][0], col_event[i_e + 1][0])
        
    col_event = col_event[start: end + 1]
    
    #------------------- create states/obs ------------------- 
    
    if cf["adv_phi"] == "":
        phi = create_phi()
    else:
        phi = q.Qobj(cf["adv_phi"])
    
    obs, laser_ops, callback = create_obs(cf)
    
    col_ops = create_collpase_ops(cf)
    
    opts = q.Odeoptions()
    opts.store_final_state = True
    
    #------------------- integration routine ------------------- 
    col_res = [] #collected results
    
    #=================== event integration ===================
    for e in col_event:
        tlist_event = np.linspace(e[0][0], e[0][1], e[0][2])
        h_idx = e[1]
        H_event = []
        for i_h in h_idx:
            H_event.append([H[i_h][0], H[i_h][2]])
        if callback:
            data = q.mesolve(H_event, phi, tlist_event, col_ops, obs, options = opts, args = args)
            phi = data.final_state
        else:
            data = q.mesolve(H_event, phi, tlist_event, col_ops, obs, options = opts, args = args)
            phi = data.final_state
            col_res.append(data)
    
    end_time = time.time()
    GREEN("integration done in {}".format(datetime.timedelta(seconds = int(end_time - start_time))))
    
    #------------------- collect result helper ------------------- 
    collect_times = lambda *args: [a for res in args for a in res.times]
    collect_expect = lambda *args: [[a for res in args for a in res.expect[j]] for j in range(len(args[0].expect))]
    
    if not callback:
        res.times, res.expect = collect_times(*col_res), collect_expect(*col_res)
    
    for l in laser_ops:
        l_res = np.zeros(len(res.times))
        for t, i_t in zipi(res.times):
            l_res[i_t] = l[0](t)
        res.expect.insert(l[1], l_res)
    
    while res.times[-1] == 0:
        res.times = res.times[:-1]
    
    for i in range(len(res.expect)): #small fix
        res.expect[i] = res.expect[i][:len(res.times)]
    
    return res.times, res.expect
