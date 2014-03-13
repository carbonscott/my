#!/usr/bin/env python
"""Topological flat model.

For simple usage, change parameters in bolcks surrounded by '#'. This script
can be used to:
    plot energy spectrum
    color edge states
    calculate flatness

TODO:
    plot DOS
    calculate berry phase
    calculate chern number

Created by Wei-Wei Luo, lww5879@gmail.com
"""

# Copyright under GNU General Public License 2010, 2012
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)

#from pythtb import * # import TB model class
import pythtb as tb  # import TB model class
import numpy as np
import matplotlib.pyplot as plt

#parameter############################################
# information used for title, figure name.
info = ['kagome lattice model', 'KG']
# define lattice vectors
lat=[[1.0,0.0],[0.5,np.sqrt(3.0)/2.0]]
# define coordinates of orbitals
orb=[[0.0, 0.0],[0.0, 0.5],[0.5, 0.0]]
# make two dimensional tight-binding model
my_model=tb.tb_model(2,2,lat,orb)

# set model parameters
# situation 1
t1 = -1.0*np.exp( 0.7j*np.pi)
t2 = 0.48*np.exp( 0.4j*np.pi)
# situation 2
t1 = -1.0*np.exp(0.0j*np.pi)
t2 = 0.44*np.exp(0.24j*np.pi)
# situation 3
t1 = -1.0*np.exp(0.22j*np.pi)
t2 = 0.19*np.exp(0.0j*np.pi)

# calculate bulk or edge spectrum
glue_edgs = True
#glue_edgs = False

# set on-site energies
my_model.set_onsite([0.0]*len(orb))
######################################################

# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
#hoppping#############################################
# add first neighbour hoppings
my_model.set_hop(t1, 0, 1, [ 0, 0])
my_model.set_hop(t1, 0, 1, [ 0,-1])
my_model.set_hop(t1, 1, 2, [ 0, 0])
my_model.set_hop(t1, 1, 2, [-1, 1])
my_model.set_hop(t1, 2, 0, [ 0, 0])
my_model.set_hop(t1, 2, 0, [ 1, 0])
# add second neighbour hoppings
my_model.set_hop(t2, 0, 1, [-1, 0])
my_model.set_hop(t2, 0, 1, [ 1,-1])
my_model.set_hop(t2, 1, 2, [ 0, 1])
my_model.set_hop(t2, 1, 2, [-1, 0])
my_model.set_hop(t2, 2, 0, [ 0, 1])
my_model.set_hop(t2, 2, 0, [ 1,-1])
######################################################

def plot_energybands(bulkE, edgeE):
    """Plot both bulk and edge energy bands, respectively.

    Parameters:
    ===========
    bulkE : ndarray
    edgeE : ndarray
    returns : figure
    """
#   fig, (ax1, ax2) = plt.subplots(1,2, figsize=(8,6))
#   for idx, (ax, e) in enumerate(zip((ax1, ax2), (bulkE, edgeE))):
#       for eigenvalues in e:
#           ax.plot(eigenvalues, 'b-')
#           if  0<idx<N-1 and (idx*n%N==0 or (idx+1)*n%N==0):
#               ax.plot(eigenvalues, "r-")
#   plt.tight_layout()

#   fig = plt.figure(figsize=(8,6))
#   ax1 = fig.add_axes([0.0, 0.0, 0.5, 1.0])
#   ax2 = fig.add_axes([0.5, 0.0, 0.5, 1.0])
#   fig, (ax1, ax2) = plt.subplots(1,2, figsize=(9,6))

    fig = plt.figure(figsize=(8,5))
    ax1 = plt.subplot(121)
    for level in bulkE:
        ax1.plot(kpts, level, 'b-')
#   plt.title(title)
    plt.title('(a)')
    plt.xlim(kpts.min(), kpts.max())
    plt.xlabel("$k$")
    plt.ylabel("$E$")

    ax2 = plt.subplot(122)
    for idx, level in enumerate(edgeE):
        if idx in [N/3, N/3-1,N*2/3, N*2/3-1]:
            ax2.plot(kpts, level, 'r-')
        else:
            ax2.plot(kpts, level, 'b-')
#   plt.title(title)
    plt.title('(b)')
    plt.xlim(kpts.min(), kpts.max())
    plt.xlabel("$k$")
#   plt.ylabel("$E$")

    fig.tight_layout()

    return fig



# make the supercell of the model
#vector = [[1,0],[0,1]] # original
#sc_model=my_model.make_supercell(vector,to_home=True)
sc_model=my_model

# now make a slab of the supercell
import copy
sc_model2 = copy.deepcopy(sc_model)
slab_model=sc_model.cut_piece(30,0,glue_edgs=glue_edgs)
slab_model2=sc_model2.cut_piece(30,0,glue_edgs=(glue_edgs+1)%2)

# visualize slab unit cell
#(fig,ax)=slab_model.visualize(0,1)
#ax.set_title("Graphene, arbitrary surface")
#ax.set_xlabel("x coordinate")
#ax.set_ylabel("y coordinate")
##fig.savefig("supercell_vis{}.pdf".format(vector))
#fig.savefig("supercell_vis.pdf")
##fig.close()

# compute the band structure in the entire band
#path=[0.0,1.0]
path=[-0.5,.5]
kpts=tb.k_path(path,100)
evals=slab_model.solve_all(kpts)
evals2=slab_model2.solve_all(kpts)

# plotting of band structure
print 'Plotting bandstructure...'

# First make a figure object
fig, ax = plt.subplots()
# plot all bands
n, N = len(orb), len(evals)
for idx, level in enumerate(evals):
    plt.plot(kpts, level, "b-")
# color edge states.
    if not glue_edgs and 0<idx<N-1 and (idx*n%N==0 or (idx+1)*n%N==0):
        plt.plot(kpts, level, "r-")

# label and title
# flatness = Gap/BandWidth
flatness = (evals[N/n+1].min()-evals[N/n-2].max())/(
    evals[N/n-2].max()-evals[0].min())
title = [
        "Edge states of {}".format(info[0]),
        'Energy band of {}, flatness = {:.1f}'.format(info[0],flatness),
        ][glue_edgs]
plt.title(title)
plt.xlim(kpts.min(), kpts.max())
plt.xlabel("$k$")
plt.ylabel("$E$")

figs = plot_energybands(evals, evals2)
#plt.show()
figs.savefig('KG_spectrum.pdf')

# make an PDF figure of a plot
figname = '{}_{}.pdf'.format(info[1], 'bulk' if glue_edgs else 'edge')
#figname = '{}_{}.png'.format(info[1], 'bulk' if glue_edgs else 'edge')
plt.savefig(figname, dpi=200,  bbox_inches='tight')
plt.close()

print 'Done.\n'
