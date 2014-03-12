#!/usr/bin/env python
"""
Prediction of energy degeneracy of edge states in each momentum
sector.

Based on:

    * `nboson` particles occupying `N` orbitals(lowest band) in flat-band limit.

    * generalized pauli principle that no more than `NBOSON` bosons
      occupy `NORBITAL` consecutive orbitals.

"""

NORBITAL = 2
NBOSON = 1
#N = 18
N = 30
nboson = 3

def config(*args):
    lst = [0]*N
    for idx in args:
        lst[idx] += 1
    return lst

import itertools
configs = (config(*sites) for sites in
        itertools.product(range(N), repeat=nboson))


import numpy as np
def genPauli(cfg):
    if max(cfg) > NBOSON: return False
    a = np.array(cfg)
    if np.all(sum(np.roll(a, i) for i in range(NORBITAL)) <= NBOSON):
        return True

#cfgs = [cfg for cfg in configs if genPauli(cfg)]
cfgs = list({tuple(cfg) for cfg in configs if genPauli(cfg)})


cfgs = np.array(cfgs)
l = np.arange(5, 5+N)
momenta = np.sum(cfgs*l, axis=1)
sorted_cfgs = cfgs[np.argsort(momenta)]
for L in range(momenta.min(), momenta.max()+1):
    print '# of states in L={:>2} sector is {:>2}'.format(L,
            np.count_nonzero(momenta == L))


