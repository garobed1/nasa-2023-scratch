import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal
from scipy.sparse import diags
from scipy.interpolate import interp1d
import copy
import sys, os
import json
'''

Find a way to parametrize a scalar profile that spans an envelope of uncertainty
in some sense. Inspired by KL expansion work in Quest, which produces a random
field for inlet boundary conditions.

Nominal function: Some best fit for the given data, profile interpolating mean
at each point?

'''

'''
Functions
'''

def unbiased_sample_cov_estimator(datat, norm = False):
    N, Ndat = datat.shape

    fac = N-1
    if norm:
        fac = N

    mat = np.zeros([N,N])
    means = np.mean(datat, axis=1)

    for j in range(Ndat):
        work = datat[:,j] - means
        mat += np.outer(work, work)/fac

    return means, mat

# produce 
def truncated_karhunen_loeve_expansion(z, eigval, eigvec):
    trunc, Ngen = z.shape
    N = copy.deepcopy(eigval.shape[0])

    # reorder
    order = eigval.argsort()[::-1]
    eigvalo = eigval[order]
    eigveco = eigvec[:, order]
    # truncate
    eigvalt = eigvalo[:trunc]
    eigvect = eigveco[:, :trunc]

    W = np.zeros([N, Ngen])
    for i in range(trunc):
        # import pdb; pdb.set_trace()
        W += np.outer(z[i,:], np.sqrt(eigvalt[i])*eigvect[:, i]).T

    return W

# import matplotlib
# matplotlib.use('tkagg') 
# generate distributed data from these distributions

name = 'testprob'
N = 33
trunc = 15
max_alt = 10000
altitudes = np.linspace(0, max_alt, N)
Ndat = 100
Ndatplot = 100
Ngen = 10
means = np.log(np.linspace(1., 50., N))
stdvs = np.linspace(3., 0.5, N)
corrfrac = 1.
# corr = max_alt/1. # correlation between adjacent points

# options if we're pulling data
prop = 'TEMP'
if len(sys.argv) > 1:
    dataname = sys.argv[1]
    name = dataname.split('/')[-1]

    with open(dataname) as fj:
        fulldata = json.load(fj)

    datap_raw = np.array(fulldata[prop]).T
    altitudes_raw = np.array(fulldata['altitude']).T
    Ndat = datap_raw.shape[1]

    # find max alt
    max_alt = np.max(altitudes_raw)

    # need to interpolate
    altitudes = np.linspace(0, max_alt, N)
    # import pdb; pdb.set_trace()
    datap = np.zeros([N, Ndat])
    for i in range(Ndat):
        datap[:,i] = np.interp(altitudes, altitudes_raw[:,i], datap_raw[:,i])

    # get means, stdvs
    means = np.mean(datap, axis=1)
    stdvs = np.std(datap, axis=1)

    datat = copy.deepcopy(datap)



else: # generate own data
    cov = np.zeros([N,N])
    for i in range(N):
        for j in range(N):
            cov[i,j] = np.exp(-abs(altitudes[i]-altitudes[j])/corr)

    datap = multivariate_normal.rvs(mean=np.zeros(N), cov=cov, size=Ndat).T
    datat = copy.deepcopy(datap)

    # scale by varying data stdvs
    for i in range(N):
        datat[i,:] *= stdvs[i]
        datat[i,:] += means[i]
        
plt.plot(datat[:,:Ndatplot], altitudes)
plt.xlabel(prop)
plt.ylabel('Altitude (1000 ft)')
plt.savefig(f'{name}_{prop}_pathsinit.png', bbox_inches="tight", dpi=500)
plt.clf()


# now we find the K-L expansion coefficients
# SPECIFICALLY FOR K(x1, x2) = e^-|x1 - x2|/corr, WE HAVE THE FOLLOWING EXACT
# QUEST MANUAL Sec 3.5.6

# ALTERNATIVELY, NUMERICALLY
meansest, covest = unbiased_sample_cov_estimator(datat, norm=False)
eigval, eigvec = np.linalg.eig(covest)

# generate independent random gaussian vectors
datag = np.zeros([trunc, Ngen])
for i in range(trunc):
    datag[i,:] = norm.rvs(size=Ngen)

# import pdb; pdb.set_trace()

pathsgen = truncated_karhunen_loeve_expansion(datag, eigval, eigvec)
pathsgent = copy.deepcopy(pathsgen)

# shift by means
for i in range(N):
    pathsgent[i,:] /= stdvs[i]
    pathsgent[i,:] += means[i]

plt.plot(pathsgent, altitudes)
plt.xlabel(prop)
plt.ylabel('Altitude (1000 ft)')
plt.savefig(f'{name}_{prop}_t{trunc}_pathsgen.png', bbox_inches="tight", dpi=500)
plt.clf()

import pdb; pdb.set_trace()



