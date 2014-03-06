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
info = ['triangular lattice model', 'TR']
# define lattice vectors
lat=[[np.sqrt(3.0)/2.0,0.5], [0.,1.]]
# define coordinates of orbitals
orb=[[0,0], [2./3.,-1./3.], [1./3.,1./3.]]
# make two dimensional tight-binding model
my_model=tb.tb_model(2,2,lat,orb)

# set model parameters
t1 = 1.0
t2 = 1.0/4
phi = np.pi/3

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
my_model.set_hop( t1                  , 0, 1, [ 0, 0])
my_model.set_hop(-t1*np.exp(1j* 2*phi), 0, 2, [ 0, 0])
my_model.set_hop( t1                  , 0, 1, [-1, 1])
my_model.set_hop(-t1                  , 0, 2, [-1, 0])
my_model.set_hop( t1                  , 0, 1, [-1, 0])
my_model.set_hop(-t1*np.exp(1j*-2*phi), 0, 2, [ 0,-1])
#
my_model.set_hop( t1                  , 1, 2, [ 1,-1])
my_model.set_hop( t1*np.exp(1j*-2*phi), 1, 2, [ 0, 0])
my_model.set_hop( t1*np.exp(1j* 2*phi), 1, 2, [ 0,-1])
# add second neighbour hoppings
my_model.set_hop( t2*np.exp(1j* 1*phi), 0, 0, [ 1, 0])
my_model.set_hop( t2*np.exp(1j* 1*phi), 0, 0, [ 0,-1])
my_model.set_hop( t2*np.exp(1j* 1*phi), 0, 0, [-1, 1])
#
my_model.set_hop( t2*np.exp(1j* 1*phi), 1, 1, [-1, 0])
my_model.set_hop( t2*np.exp(1j* 1*phi), 1, 1, [ 0, 1])
my_model.set_hop( t2*np.exp(1j* 1*phi), 1, 1, [ 1,-1])
#
my_model.set_hop(-t2                  , 2, 2, [ 1, 0])
my_model.set_hop(-t2                  , 2, 2, [ 1,-1])
my_model.set_hop(-t2                  , 2, 2, [ 0,-1])
######################################################


# make the supercell of the model
#vector = [[1,0],[0,1]] # original
#sc_model=my_model.make_supercell(vector,to_home=True)
sc_model=my_model

# now make a slab of the supercell
slab_model=sc_model.cut_piece(30,0,glue_edgs=glue_edgs)

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

# make an PDF figure of a plot
figname = '{}_{}.pdf'.format(info[1], 'bulk' if glue_edgs else 'edge')
figname = '{}_{}.png'.format(info[1], 'bulk' if glue_edgs else 'edge')
plt.savefig(figname, dpi=200,  bbox_inches='tight')
plt.close()

print 'Done.\n'
