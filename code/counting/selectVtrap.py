"""Select suitable parameter 'vtrap' for investigation of edge states in FCI.

source data format:
        vtrap sector energy

"""

import sys
import numpy as np
src = sys.argv[1:2] or 'dat'
dat = np.genfromtxt(src)

# select vtrap where the lowest energy sequence should be 1,0,5,4,3.2.
traps = []
#for vtrap in np.arange(101)*0.0005:
#for vtrap in np.arange(3,31)*0.001:
# For floating point comparison, use `abs(a-b)<1e-10` instead of `a==b` (often
# return False)
#   data = dat[abs(dat[:,0]-vtrap) < 1e-10]
#   s0 = data[data[:,1] == 0][0, 2]
#   s1 = data[data[:,1] == 1][0, 2]
#   s2 = data[data[:,1] == 2][0, 2]
#   s3 = data[data[:,1] == 3][0, 2]
#   s4 = data[data[:,1] == 4][0, 2]
#   s5 = data[data[:,1] == 5][0, 2]
#   if s1 < s0 < s5 < s4 < s3 < s2:
#       s0 = data[data[:,1] == 0][1, 2]
#       s1 = data[data[:,1] == 1][1, 2]
#       s2 = data[data[:,1] == 2][1, 2]
#       s3 = data[data[:,1] == 3][1, 2]
#       s4 = data[data[:,1] == 4][1, 2]
#       s5 = data[data[:,1] == 5][1, 2]
#       if s1 < s0 < s5 < s4 < s3 < s2:
#           s0 = data[data[:,1] == 0][2, 2]
#           s1 = data[data[:,1] == 1][2, 2]
#           s2 = data[data[:,1] == 2][2, 2]
#           s3 = data[data[:,1] == 3][2, 2]
#           s4 = data[data[:,1] == 4][2, 2]
#           s5 = data[data[:,1] == 5][2, 2]
#           if s1 < s0 < s5 < s4 < s3 < s2:
#              traps.append(vtrap)
#              traps.append(round(vtrap, 4))
#
#print traps
#print '# of selected traps is {}'.format(len(traps))

#if not raw_input('plot figures? (Enter to quit)'):
#    sys.exit()

import matplotlib.pyplot as plt
vtrap, sector, energy = dat[:,0], dat[:,1], dat[:,2]
for trap in np.arange(5,6)*0.001:
#for trap in np.arange(4,10)*0.001:
#for trap in np.arange(10,20)*0.001:
#for trap in traps:
    x, y = sector[abs(vtrap-trap)<1e-8], energy[abs(vtrap-trap)<1e-8]
    if x.size > 200*6:
        raise Exception("\n\nToo many energy states to handle!")


#   x, y = sector[abs(vtrap-trap)<1e-8], energy[abs(vtrap-trap)<1e-8]
    fig, ax = plt.subplots(
                           figsize=(6,7)
                           )

# for each sector, annote levels with index with large gap below it.
#   for i in range(6):
#       labels = (y[x==i][1:]-y[x==i][:-1]).argsort()[-4:] + 1
#       for idx, (x0, y0) in enumerate(zip(x[x==i][labels], y[x==i][labels])):
#           ax.text(x0, y0, '{:d}'.format(labels[idx]+1))

#   ax.plot(x, y, '_', ms = 12, mew = .5, alpha=.8)
#   ax.plot(x, y, '_', ms = 12, mew = 1)
    ax.plot((6-x)%6, y, '_', ms = 12, mew = 1)
#   ax.plot((6-x)%6, y, 'ro', ms=1)
#   ax.set_title(fin)
    ax.set_xlim(-0.2,5.2)
    ax.set_ylim(y.min()-0.005, y.max())
#   ax.set_ylim(-6.94,-6.87)
    ax.set_xlabel(r'$k$')
    ax.set_ylabel(r'$E$')
#ax.grid(True)

    fig.tight_layout()
#fig.show()
#plt.show() #This should be used in non-interactive env for display.

#   fig.savefig('text_{:.3f}.pdf'.format(trap) )
    fig.savefig('lapack_edge.pdf')
#   fig.savefig('text_{:.3f}.png'.format(trap) )
#   fig.savefig('vtrap{:.3f}.png'.format(trap) )
    plt.close()



