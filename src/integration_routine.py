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

from addon import *

#=================== global variables ===================
default_cutoff = 100000
"""
the default cutoff frequency used by :func:`create_hamiltonian`
"""

#=================== classes ===================
def syn(expr, pat, args):
    """
    a helper function of the :func:`create_hamiltonian`. fuses the physical equation part with the laser-pattern part to a new python function
    
    :param expr: the physical equation in string representation
    :param args: the dictionary that holds the values to the physical expression
    :param pat: the laser pattern, see :attr:`laser.pattern`
    :returns: a python function that takes the time and returns the physical expression if t is in the active pattern and 0 else
    """
    from pylab import exp
    from pylab import pi
    
    for k, v in args.iteritems():
        exec(k + " = v") in locals()
    
    def fct(t, pat):
        for p in pat:
            if t > p[0] and t < p[1]:
                return 1
        return 0
        
    exec("g = lambda t, args: (" + expr + ") * fct(t, pat)") in locals()
    return g

def create_hamiltonian(atoms, vibrons, lasers, _rabi, _eta, **kwargs):
    """
    :param atoms: a list of :class:`atom`
    :param vibrons: a list of :class:`vibron`
    :param lasers: a list of :class:`laser`
    :param _rabi: a list containing the rabi-frequency with the following dimensions: list[N_lasers][N_atoms][N_transitions_per_atom]. since the amount of transitions can be different for each atom, the last dimension will depend on the atom and is not the same for all atoms.
    :param _eta: a list containing the lamb-dicke parameter with the following dimensions: list[N_atoms][N_vibrons]
    :param **kwargs: key word arguments for additional optional parameters like cutoff. If no cutoff is given, it will be set to `default_cutoff`.
    :returns: the hamiltonian for the integration
    """
    global vibron_holder_
    
    if "cutoff" not in kwargs:
        kwargs["cutoff"] = default_cutoff
    cutoff = kwargs["cutoff"]
    
    H = []
    
    rwa_takes = [0, 0]
    
    for l, i_l in zipi(lasers):
        for a, i_a in zipi(atoms):
            for t, i_t in zipi(a.trans):
                args = {}
                args["hb"] = 1
                args["Ohm"] = _rabi[i_l][i_a][i_t]
                args["w0"] = t[0]
                args["wL"] = l.freq
                args["fL"] = l.phase
                
                #------------------- set up the operators ------------------- 
                d_internal = t[1]
                ci = q.tensor(d_internal.dag(), q.qeye(vibron_holder_.N))
                di = q.tensor(d_internal      , q.qeye(vibron_holder_.N))
                
                Hc= [ [ci, syn(" 0.5 * hb * Ohm * ( exp( 1j * ((w0 - wL) * t + fL)) )", l.pattern, args)]
                    , [di, syn(" 0.5 * hb * Ohm * ( exp(-1j * ((w0 - wL) * t + fL)) )", l.pattern, args)]
                    ]
                #------------------- rotating wave approx ------------------- 
                rwa_c = [ abs(args["w0"] - args["wL"])
                        , abs(args["w0"] - args["wL"])
                        ]

                #------------------- add carrier terms ------------------- 
                if args["Ohm"] == 0: #then the term will be zero anyway
                    continue
                    
                rwa_takes[1] += len(rwa_c)
                for h, diff in zip(Hc, rwa_c):
                    if diff < cutoff:
                        H.append(h)
                        rwa_takes[0] += 1
                    
                for v, i_v in zipi(vibrons):
                    #------------------- set up the constants ------------------- 
                    args["eta"] = _eta[i_a][i_v]
                    args["wz"] = v.freq
                    d_osz = v.destroy
                    
                    #------------------- set up the operators ------------------- 
                    cc = q.tensor(d_internal.dag(), d_osz.dag())
                    cd = q.tensor(d_internal.dag(), d_osz)
                    dc = q.tensor(d_internal      , d_osz.dag())
                    dd = q.tensor(d_internal      , d_osz)
                    
                    Hs= [ [cc, syn(" 0.5 * hb * Ohm * ( exp( 1j * ((w0 - wL) * t + fL)) ) * 1j * eta * exp( 1j * wz * t)", l.pattern, args)]
                        , [cd, syn(" 0.5 * hb * Ohm * ( exp( 1j * ((w0 - wL) * t + fL)) ) * 1j * eta * exp(-1j * wz * t)", l.pattern, args)]
                        , [dc, syn("-0.5 * hb * Ohm * ( exp(-1j * ((w0 - wL) * t + fL)) ) * 1j * eta * exp( 1j * wz * t)", l.pattern, args)]
                        , [dd, syn("-0.5 * hb * Ohm * ( exp(-1j * ((w0 - wL) * t + fL)) ) * 1j * eta * exp(-1j * wz * t)", l.pattern, args)]
                        ]
                    
                    #------------------- rotating wave approx ------------------- 
                    rwa_s = [ abs(args["w0"] - args["wL"] + args["wz"])
                            , abs(args["w0"] - args["wL"] - args["wz"])
                            , abs(args["w0"] - args["wL"] - args["wz"])
                            , abs(args["w0"] - args["wL"] + args["wz"])
                            ]
                    #------------------- add sideband terms ------------------- 
                    if args["eta"] == 0:
                        continue
                        
                    rwa_takes[1] += len(rwa_s)
                    for h, diff in zip(Hs, rwa_s):
                        if diff < cutoff:
                            H.append(h)
                            rwa_takes[0] += 1
    print("{1} of {2} terms survived the RWA (cutoff: {0})".format(cutoff, *rwa_takes))
    return H

def create_phi():
    """
    creates the initial global state phi by adding (tensor) the atom and vibron states
    
    :returns: a `qutip.Qobj` that contains the density matrix of the global state
    """
    return q.tensor(atom_holder_.state(), vibron_holder_.state()).unit()

def create_obs():
    """
    creates the default observables. they just measure the population in each state of the atoms and vibrons.
    
    :returns: all observables for the populations of each atom and vibron state
    """
    obs = []
    
    dims = atom_holder_.current_dims + vibron_holder_.current_dims
    
    for d, i_d in zipi(dims):
        for l in range(d):
            obs_in = q.fock_dm(d, l)
                
            for k in dims[:i_d]:
                obs_in = q.tensor(q.qeye(k), obs_in)
                
                
            for k in dims[i_d+1:]:
                obs_in = q.tensor(obs_in, q.qeye(k))
                
            obs.append(obs_in)
    
    return obs
    
def create_labels(cf):
    """
    helper function for :func:`plot`
    
    just specifies the color and linestyle of the plots and extracts the propper names from the cf
    
    :param cf: canonical form of the integration problem
    :returns: the label as well as color list for plotting
    """
    lbl = []
    col = []
    
    color = ["b-", "g-", "c-", "m-", "r-", "y-", "b:", "g:", "c:", "m:", "r:", "y:"]
    
    n_atoms = len(atom_holder_.current_dims)
    
    for a, i_a in zipi(atom_holder_.current_dims):
        for lvl in range(a):
            if cf == "none":
                lbl.append("atom " + str(i_a) + " s" +  str(lvl))
            else:
                lbl.append(cf["atom"][i_a][0] + " s" +  str(lvl))
            col.append(color[lvl])
    
    for v, i_v in zipi(vibron_holder_.current_dims):
        for lvl in range(v):
            if cf == "none":
                lbl.append("vib " + str(i_v) + " s" +  str(lvl))
            else:
                lbl.append(cf["vibron"][i_v][0] + " s" +  str(lvl))
            col.append(color[lvl])
    
#~ '-' 	solid line style
#~ '--'	dashed line style
#~ '-.'	dash-dot line style
#~ ':' 	dotted line style
#~ '.' 	point marker
#~ ',' 	pixel marker
#~ 'o' 	circle marker
#~ 'v' 	triangle_down marker
#~ '^' 	triangle_up marker
#~ '<' 	triangle_left marker
#~ '>' 	triangle_right marker
#~ '1' 	tri_down marker
#~ '2' 	tri_up marker
#~ '3' 	tri_left marker
#~ '4' 	tri_right marker
#~ 's' 	square marker
#~ 'p' 	pentagon marker
#~ '*' 	star marker
#~ 'h' 	hexagon1 marker
#~ 'H' 	hexagon2 marker
#~ '+' 	plus marker
#~ 'x' 	x marker
#~ 'D' 	diamond marker
#~ 'd' 	thin_diamond marker
#~ '|' 	vline marker
#~ '_' 	hline marker
#~ 
#~ The following color abbreviations are supported:
#~ character 	color
#~ ‘b’ 	blue
#~ ‘g’ 	green
#~ ‘r’ 	red
#~ ‘c’ 	cyan
#~ ‘m’ 	magenta
#~ ‘y’ 	yellow
#~ ‘k’ 	black
#~ ‘w’ 	white
    
    return lbl, col

def plot(times, exp, cf = "none"):
    #=================== plot ===================
    
    lbl, col = create_labels(cf)
    
    
    dims = atom_holder_.current_dims + vibron_holder_.current_dims
    
    f, ax = pl.subplots(len(dims))
    
    cur = 0
    idx = 0
    pl.xlabel("Time")
    pl.ylabel("exp val")
    pl.subplot(len(dims), 1, idx)
    for ex in range(len(exp)):
        if ex >= (cur + dims[idx]):
            cur += dims[idx]
            idx += 1
            pl.subplot(len(dims), 1, idx)
        
        pl.plot(times, exp[ex], col[ex])
        pl.legend(lbl[cur:(cur + dims[idx])], loc = 5)
        pl.ylim([-0.005, 1.005])
        
    pl.show()

import time
import datetime
def integrate_cf(cf):
    GREEN("start integration")
    start_time = time.time()
    at = []
    for a in cf["atom"]:
        at.append(atom(a[1], a[2]))
    
    vb = []
    for v in cf["vibron"]:
        vb.append(vibron(v[1], v[2], v[3]))
    
    ls = []
    for l in cf["laser"]:
        ls.append(laser(l[3], l[1], l[2]))
        
    rabi = cf["rabi"]
    eta = cf["eta"]
    
    H = create_hamiltonian(at, vb, ls, rabi, eta, cutoff = cf["cutoff"])
    phi = create_phi()
    obs = create_obs()
    
    tlist = np.linspace(cf["lower"], cf["upper"], cf["measure"])
    opts = q.Odeoptions(max_step = 0.01)
    
    res = q.mesolve(H, phi, tlist, [], obs, options = opts)
    
    end_time = time.time()
    GREEN("integration done in {}".format(datetime.timedelta(seconds = int(end_time - start_time))))
    
    plot(res.times, res.expect, cf)
    
    atom_holder_.clear()
    vibron_holder_.clear()

def create_hamiltonian_eff(atoms, vibrons, lasers, _rabi, _eta, **kwargs):
    """
    :param atoms: a list of :class:`atom`
    :param vibrons: a list of :class:`vibron`
    :param lasers: a list of :class:`laser`
    :param _rabi: a list containing the rabi-frequency with the following dimensions: list[N_lasers][N_atoms][N_transitions_per_atom]. since the amount of transitions can be different for each atom, the last dimension will depend on the atom and is not the same for all atoms.
    :param _eta: a list containing the lamb-dicke parameter with the following dimensions: list[N_atoms][N_vibrons]
    :param **kwargs: key word arguments for additional optional parameters like cutoff. If no cutoff is given, it will be set to `default_cutoff`.
    :returns: the hamiltonian for the integration
    """
    global vibron_holder_
    
    if "cutoff" not in kwargs:
        kwargs["cutoff"] = default_cutoff
    cutoff = kwargs["cutoff"]
    
    H = []
    col_args = {}
    col_args["hb"] = 1
    
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

def integrate_cf_eff(cf, **kwargs):
    GREEN("start integration")
    start_time = time.time()
    
    #------------------- checking kwargs ------------------- 
    #the two possibilities are "event" and "subspace"
    if not "opt" in kwargs:
        kwargs["opt"] = "event"
    
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
        
    rabi = cf["rabi"]
    eta = cf["eta"]
    
    H, args = create_hamiltonian_eff(at, vb, ls, rabi, eta, cutoff = cf["cutoff"])
    
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
    
    #------------------- reduce H-space as much as possible for each event ------------------- 
    all_dims = H[0][0].shape[0]
    
    for i_e in range(len(col_event)):
        e = col_event[i_e][0]
        i_h_list = col_event[i_e][1]
        #in the beginning it is assumed that each hamiltonian only affects one subsystem
        sys = [set([i]) for i in range(len(H[0][0].dims[0]))]
        affecetd_sys = []
        for i_h in i_h_list:
            affecetd_sys.append(set(H[i_h][3]))
        
        #now the list of sets sys gets corrected, bc hamiltonians can affect more than one subsystem (RSB/BSB)
        for d in affecetd_sys:
            found = False
            for s in sys:
                if d <= s:
                    found = True
                    break
            if not found:
                temp = set()
                for el in d:
                    for s in sys:
                        if el in s:
                            temp = temp | s
                            sys.remove(s)
                            break
                sys.append(temp)
        
        #and the sets in sys that don't have a hamiltonian get eliminated
        used = set([])
        for d in affecetd_sys:
            used = used | d
        
        for s in reversed(sys):
            if not s <= used:
                sys.remove(s)
        #convert back to list
        for i_s in range(len(sys)):
            sys[i_s] = list(sys[i_s])
        
        split_act = [[] for i in range(len(sys))]
        
        for i_h in i_h_list:
            for i_s in range(len(sys)):
                if set(H[i_h][3]) <= set(sys[i_s]):
                    split_act[i_s].append(i_h)
                    break
        #now the i_h_list is split into independent subspaces (e.g. integrating N/2 dim twice is way faster than N once)
        col_event[i_e][1] = zip(sys, split_act)
    
    #col_event is a list with the following content:
    #col_event[N][0]: event time
    #col_event[N][1]: a list of which subsystems are affected by which hamiltonians
    #col_event[N][1][M][0]: the affected subsystems
    #col_event[N][1][M][1]: the acting hamiltonians (the index, the hamiltonian is H[index])
    
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
    
    phi = create_phi()
    obs = create_obs()
    
    if "obs" in cf.keys():
        for ob in cf["obs"]:
            fct_name = ob.split("(")[0]
            fct_name = fct_name.split("def ")[1]
            
            exec(ob) in globals()
            
            #~ print(eval(fct_name)(q.tensor(q.qeye(2), q.qeye(2), q.fock_dm(5, 0)).unit()))
            
    opts = q.Odeoptions()
    opts.rec = 1
    #------------------- collect result helper ------------------- 
    #a list with the subsystem dimensions
    full_dims = H[0][0].dims[0]
    #the starting position in obs for each system
    full_start = list(np.add.accumulate(full_dims))[:-1]
    full_start.insert(0, 0)
    #just a set with all subsystems
    full_space = set(range(len(full_dims)))
    #multiplies dimensions of subsystems that are not in sub (needed bc of ptrace) / if full_space == sub the return 1
    dim_corr = lambda sub: np.multiply.accumulate([[full_dims[i] for i in (full_space - set(sub))], [1]][full_space == set(sub)])[-1]
    
    collect_times = lambda *args: [a for res in args for a in res.times]
    collect_expect = lambda *args: [[a for res in args for a in res.expect[j]] for j in range(len(args[0].expect))]
    
    #since phi is updated on a subspace, one has to update the global state before processing the next event
    def fix_phi(sub_space, new_state, old_state):
        old = [old_state.ptrace(s) for s in full_space]
        new = [new_state.ptrace(s) for s in range(len(sub_space))]
        for i in range(len(sub_space)):
            old[sub_space[i]] = new[i]
        return q.tensor(old)
    
    #------------------- integration routine ------------------- 
    col_res = []
    
    for e in col_event:
        tlist_event = np.linspace(e[0][0], e[0][1], e[0][2])
        
        if kwargs["opt"] == "subspace":
            #~ #=================== improved subspace integration ===================
            obs_val = [[q.expect(o, phi)]*e[0][2] for o in obs]
            if len(e[1]) == 0:
                res = q.mesolve([], phi, tlist_event, [], [], options = opts, args = args)
                phi = res.states[-1]
            else:
                for sub in e[1]:
                    sub_space = sub[0]
                    h_idx = sub[1]
                    H_event = []
                    for i_h in h_idx:
                        H_event.append([H[i_h][0].ptrace(sub_space) / dim_corr(sub_space), H[i_h][2]])
                    
                    obs_event = []
                    for sub_sys in sub_space:
                        for j in range(full_dims[sub_sys]):
                            obs_event.append(obs[full_start[sub_sys] + j].ptrace(sub_space) / dim_corr(sub_space))
                    phi_event = phi.ptrace(sub_space)
                    res = q.mesolve(H_event, phi_event, tlist_event, [], obs_event, options = opts, args = args)
                    
                    n = 0
                    for sub_sys in sub_space:
                        for j in range(full_dims[sub_sys]):
                            obs_val[full_start[sub_sys] + j] = res.expect[n]
                            n += 1
                    
                    phi = fix_phi(sub_space, res.states[-1], phi)
            
            res.expect = obs_val
        elif kwargs["opt"] == "event":
            #=================== event integration ===================
            h_idx = [y for x in e[1] for y in x[1]]
            H_event = []
            for i_h in h_idx:
                H_event.append([H[i_h][0], H[i_h][2]])
            res = q.mesolve(H_event, phi, tlist_event, [], obs, options = opts, args = args)
            phi = res.states[-1]
        else:
            print("wrong kwarg for opt given")
        col_res.append(res)
    
    end_time = time.time()
    GREEN("integration done in {}".format(datetime.timedelta(seconds = int(end_time - start_time))))
    
    plot(collect_times(*col_res), collect_expect(*col_res), cf)
    
    atom_holder_.clear()
    vibron_holder_.clear()
