#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 23:35:11 2022
This is a supplementary script that visualizes a 3D cubic lattice, using data from main.py  
@author: mrahm32
"""

import itertools
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import color as color

def generate_fig_parameters(cluster_num):
    colors = ["#b6174b","#cbf3d2","#ff0000", "#4e937a","#ffcc00", "#7b287d","#37c871", "#b7c0ee","#0066ff", "#957DAD"]
    colors = ['#A5be00', '#d81e5b','#331832', '#b892ff', '#29bf12','#08bdbd', '#ff9914', '#17bebb', '#003049', '#2a0800', '#314cb6', "#b6174b","#cbf3d2","#ff0000", "#4e937a","#ffcc00", "#7b287d","#37c871", "#b7c0ee","#0066ff", "#957DAD"]
    color_gradient = color.polylinear_gradient(colors,n=10)
    color_gradient = color_gradient['hex']
    markers = ['o', 's', 'v', '*', 'x']
    return colors, markers

def make_vis(networks, L, path):
    '''parameters: networks = dictionary of networks with their labels and 3D coordinates of the members. 
    L = dimension of cube of sized L^3
    path = all the points visited in order to determine that there is a spanning network. this is for just 
    visual confirmation that there is indeed a path.'''
    
    path = np.array(path)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    network_num = len(list(networks.keys()))
    colors, markers = generate_fig_parameters(network_num)
    for label in networks.keys():
        if label != 0:
            indx = list(networks.keys()).index(label)-1
            triplets = networks[label]
            triplets = np.array(triplets)
            col = colors[indx]
            ax.scatter(triplets[:,0], triplets[:,1], triplets[:,2], marker = 's', s=150,c = col)
            #ax.scatter(triplets[:,0], triplets[:,1], triplets[:,2], marker = 's', s=150,c ='black')
            ax.scatter(path.T[0],path.T[1],path.T[2], marker = 'o', c='black', s= 100)
            ax.set_xlim3d(0,L-1)
            ax.set_ylim3d(0,L-1)
            ax.set_zlim3d(0,L-1)
    plt.show()
    fig.savefig(f'percolation_network.pdf')
