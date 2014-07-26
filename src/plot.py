#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    23.07.2014 23:39:41 CEST
# File:    plot.py

from .src_import import *

def create_labels(cf):
    """
    helper function for :func:`plot`
    
    just specifies the color and linestyle of the plots and extracts the propper names from the cf
    
    :param cf: canonical form of the integration problem
    :returns: the label as well as color list for plotting
    """
    flatten = lambda l: [y for x in l for y in x]
    
    lbl = []
    col = []
    
    color = ["b-", "g-", "c-", "m-", "r-", "y-", "b:", "g:", "c:", "m:", "r:", "y:"]
    for subplot in cf["plot"]:
        lbl += [x[1] for x in subplot]
        col += color[:len(subplot)]
    
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

def plot(times, exp, cf):
    """
    plots the data
    
    just specifies the color and linestyle of the plots and extracts the propper names from the cf
    
    :param times: time points measured: list[#measurements]
    :param exp: for each observed/measured property this list contains a list with the results, list[#obs][#measurements]
    :param cf: canonical form of the integration problem. especially the key "plot" contains the information what obs are plotted together and which need an own subplot
    :returns: None
    """
    #=================== plot ===================
    
    lbl, col = create_labels(cf)
    
    
    dims = [len(x) for x in cf["plot"]]
    
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
