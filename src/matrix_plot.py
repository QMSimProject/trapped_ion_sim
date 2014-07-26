#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    03.06.2014 18:54:37 CEST
# File:    matrix_plot.py

#extracted/modified from the Bell state example of the QuTip framework

# Creation and manipulation of a Bell state with
# 3D histogram plot output.

from .src_import import *
import matplotlib as mpl
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

def l_bra(s):
    return '$\left|' + s + '\\rangle\\right.$'

def l_ket(s):
    return '$\left\\langle' + s + '|\\right.$'

def qubit_hist(Q, xlabels = "auto", ylabels = "auto", title = "", save_name = "none"):
    #------------------- generate labels automatically ------------------- 
    d = Q.dims[0]
    autolabels = []
    label = ""
    def recursion(lst, current_label):
        if len(lst) == 1:
            for i in range(lst[0]):
                autolabels.append(l_bra(current_label + str(i)))
        else:
            for i in range(lst[0]):
                recursion(lst[1:], current_label + str(i))
            
    recursion(d, label)
    
    #------------------- assing if not provided ------------------- 
    if xlabels == "auto":
        xlabels = autolabels
    if ylabels == "auto":
        ylabels = autolabels
    
    
    # Plots density matrix 3D histogram from Qobj
    # xlabels and ylabels are list of strings for axes tick labels
    num_elem = np.prod(Q.shape)       # num. of elements to plot
    xpos, ypos = np.meshgrid(range(Q.shape[0]), range(Q.shape[1]))
    xpos = xpos.T.flatten() - 0.5  # center bars on integer value of x-axis
    ypos = ypos.T.flatten() - 0.5  # center bars on integer value of y-axis
    zpos = np.zeros(num_elem)         # all bars start at z=0
    dx = 0.8 * np.ones(num_elem)      # width of bars in x-direction
    dy = dx.copy()  # width of bars in y-direction (same as x-dir here)
    dz = np.real(Q.full().flatten())  # height of bars from density matrix

    # generate list of colors for probability data
    # add +-0.1 in case all elements are the same (colorbar will fail).
    nrm = mpl.colors.Normalize(0, 1)
    #~ nrm = mpl.colors.Normalize(min(dz) - 0.1, max(dz) + 0.1)
    colors = cm.jet(nrm(dz))

    # plot figure
    fig = pl.figure(1)
    ax = Axes3D(fig, azim=-50, elev=55)
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors)

    # set x-axis tick labels and label font size. set x-ticks to integers
    ax.axes.w_xaxis.set_major_locator(pl.IndexLocator(1, -0.5))
    ax.set_xticklabels(xlabels)
    ax.tick_params(axis='x', labelsize=18)

    # set y-axis tick labels and label font size. set y-ticks to integers
    ax.axes.w_yaxis.set_major_locator(pl.IndexLocator(1, -0.5))
    ax.set_yticklabels(ylabels)
    ax.tick_params(axis='y', labelsize=18)

    # remove z-axis tick labels by moving them outside the plot range
    ax.axes.w_zaxis.set_major_locator(pl.IndexLocator(2, 2))
    # set the plot range in the z-direction to fit data
    ax.set_zlim3d([0, 1])
    #~ ax.set_zlim3d([min(dz) - 0.1, max(dz) + 0.1])
    pl.title(title)
    # add colorbar with color range normalized to data
    cax, kw = mpl.colorbar.make_axes(ax, shrink=.75, pad=.02)
    cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cm.jet, norm=nrm)
    cb1.set_label("Probability", fontsize=14)
    if save_name == "none":
        pl.show()
        pl.close()
    else:
        pl.savefig(save_name)
