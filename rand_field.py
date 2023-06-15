import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal
from scipy.sparse import diags
import copy
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
N = 100
trunc = 5
max_alt = 10000
altitudes = np.linspace(0, max_alt, N)
Ndat = 100
Ngen = 10
means = np.log(np.linspace(1., 50., N))
stdvs = np.linspace(3., 0.5, N)
# corr = max_alt/1. # correlation between adjacent points
corrfrac = 1.5

corr = corrfrac*max_alt
data = np.zeros([N, Ndat])
for i in range(N):
    data[i,:] = norm.rvs(loc=means[i], scale=stdvs[i], size=Ndat)
# import pdb; pdb.set_trace()
plt.scatter(data, np.tile(altitudes, Ndat).reshape(Ndat, N).T)
plt.savefig('altscatter.png', bbox_inches="tight")
plt.clf()

# now generate correlated paths from the same distributions
# two point correlation function K(x1, x2) = e^-|x1 - x2|/corr
# corrents = np.zeros(N-1)
# for i in range(N-1):
#     corrents[i] = np.exp(-abs(altitudes[i+1]-altitudes[i])/corr)
# cov = diags([stdvs, corrents, corrents], [0, -1, 1]).toarray()
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
# import pdb; pdb.set_trace()
plt.plot(datat, altitudes)
plt.savefig('altpathsinit.png', bbox_inches="tight", dpi=500)
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
    pathsgent[i,:] += means[i]

plt.plot(pathsgent, altitudes)
plt.savefig('altpathsgen.png', bbox_inches="tight", dpi=500)
plt.clf()

# import pdb; pdb.set_trace()



