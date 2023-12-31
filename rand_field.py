import numpy as np
from scipy.stats import norm, multivariate_normal, lognorm
from scipy.sparse import diags
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde

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

def unbiased_sample_cov_estimator(datat, norm = True):
    N, Ndat = datat.shape

    # was being silly here
    fac = Ndat-1
    if norm:
        fac = Ndat

    mat = np.zeros([N,N])
    means = np.mean(datat, axis=1)

    for j in range(Ndat):
        work = datat[:,j] - means
        outer = np.outer(work, work)
        mat += outer/fac
        # import pdb; pdb.set_trace()

    # if 0:
    #     import matplotlib.pyplot as plt
    #     from matplotlib import colors
    #     from matplotlib import cm
    #     plt.rcParams['font.size'] = 16

    #     # make indices altitude?
    #     cmap = cm.coolwarm
    #     plt.matshow(mat, cmap=cmap, norm=colors.CenteredNorm())
    #     plt.xlabel('Altitude Index')
    #     plt.ylabel('Altitude Index')
    #     plt.colorbar()
    #     plt.savefig('covmat.png', bbox_inches='tight')
    #     plt.clf()
    #     plt.scatter(np.arange(0,N), np.linalg.eig(mat)[0])
    #     plt.xlabel('Sorted Eigenvalue Index')
    #     plt.ylabel('Eigenvalue')
    #     plt.savefig('eigdecay.png', bbox_inches='tight')
    #     plt.clf()
    #     import pdb; pdb.set_trace()
    return means, mat

# produce 
def truncated_karhunen_loeve_expansion(z, means, eigval, eigvec, prop=None):
    # prop allows for specific variable transformation/post processing

    trunc, Ngen = z.shape
    N = copy.deepcopy(eigvec.shape[0])

    # reorder
    order = eigval.argsort()[::-1]
    eigvalo = eigval[order]
    eigveco = eigvec[:, order]
    # truncate
    eigvalt = eigvalo[:trunc]
    eigvect = eigveco[:, :trunc]

    W = np.zeros([N, Ngen])
    for i in range(trunc):
        W += np.outer(z[i,:], np.sqrt(eigvalt[i])*eigvect[:, i]).T

    # shift by means
    for i in range(N):
        W[i,:] += means[i]
    
    # # humidity post processing
    # if prop == 'HUMIDITY':
    #     W *= W
    #     W = np.clip(W, a_min=None, a_max=100.)

    return W



def preprocess_data(datapath, prop, N=0):
    
    name = datapath.split('/')[-1]

    with open(datapath) as fj:
        fulldata = json.load(fj)

    datap_raw = np.array(fulldata[prop]).T
    altitudes_raw = np.array(fulldata['altitude']).T
    Ndat = datap_raw.shape[1]


    # find max alt
    max_alt = np.max(altitudes_raw)

    # need to interpolate
    if N:
        altitudes = np.linspace(0, max_alt, N)
    else:
        #EDIT: instead of linspace, take the average of each measured altitude
        altitudes = np.mean(altitudes_raw, axis=1)
        N = altitudes.shape[0]

    # import pdb; pdb.set_trace()
    datap = np.zeros([N, Ndat])
    for i in range(Ndat):
        datap[:,i] = np.interp(altitudes, altitudes_raw[:,i], datap_raw[:,i])
        datap[:,i] = interp1d(altitudes_raw[:,i], datap_raw[:,i], fill_value='extrapolate')(altitudes)

    # import pdb; pdb.set_trace()
    # transform to sqrt(HUMIDITY) for KL purposes
    # if prop == "HUMIDITY":
    #     datap = np.clip(datap, 0., None)
    #     datap = np.sqrt(datap)
    # get means, stdvs
    means = np.mean(datap, axis=1)
    stdvs = np.std(datap, axis=1)

    return altitudes, datap, means, stdvs, name

        

def get_kl_coefficients(datat, norm=False):
    # now we find the K-L expansion coefficients
    # SPECIFICALLY FOR K(x1, x2) = e^-|x1 - x2|/corr, WE HAVE THE FOLLOWING EXACT
    # QUEST MANUAL Sec 3.5.6

    # ALTERNATIVELY, NUMERICALLY
    meansest, covest = unbiased_sample_cov_estimator(datat, norm=norm)
    eigval, eigvec = np.linalg.eig(covest)
    return eigval, eigvec, covest


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.rcParams['font.size'] = 16

    # needed parameters
    name = 'testprob'
    N = 33
    trunc = 4
    # options if we're pulling data
    prop = 'HUMIDITY' #'TEMP'
    #AND DATA AS COMMAND LINE ARGUMENT

    Ndatplot = 200
    Ngen = 5

    # manufactured data
    exclude = 4
    Ndat = 100
    max_alt = 10000
    corrfrac = 1.
    altitudes = np.linspace(0, max_alt, N)
    means = np.log(np.linspace(1., 50., N))
    stdvs = np.linspace(3., 0.5, N)


    # corr = max_alt/1. # correlation between adjacent points




    # options for generating new paths
    # will be determined by QUEST in practice
    gstd = 1.0 # std deviation of the random KL coeffs
    # or do list of length trunc
    # gstd = [0.2, 0.2, 0.4, 0.5] # seems ok for temp, trunc 4?
    # gstd = [0.8, 0.9, 0.1, 0.1]


    # preprocess data
    if len(sys.argv) > 1:
        datapath = sys.argv[1]
        altitudes, datat, means, stdvs, name = preprocess_data(datapath, prop)
        if exclude > 0:
            altitudes = altitudes[:-exclude]
            datat = datat[:-exclude, :]
            means = means[:-exclude]
            stdvs = stdvs[:-exclude]
        Ndat = datat.shape[1]
        N = datat.shape[0]
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

    # get kl coefficients
    eigval, eigvec = get_kl_coefficients(datat, norm=False)

    # generate independent random gaussian vectors for this example
    datag = np.zeros([trunc, Ngen])
    if not isinstance(gstd, list):
        gstd = trunc*[gstd]
    for i in range(trunc):
        datag[i,:] = norm.rvs(size=Ngen, scale=gstd[i])
        # datag[i,:] = lognorm.rvs(size=Ngen, s=gstd[i])

    # generate new paths
    pathsgent = truncated_karhunen_loeve_expansion(datag, means, eigval, eigvec, prop)

    # transform some humidity back from square root space
    # if prop == 'HUMIDITY':
    #     datat *= datat
    #     means = np.mean(datat, axis=1)
    #     stdvs = np.std(datat, axis=1)

    meanpaths = np.mean(pathsgent, axis=1)
    stdvpaths = np.std(pathsgent, axis=1)

    pprop2label = {
        "TEMP": "Temperature (F)",
        "HUMIDITY": "Relative Humidity (%)",
        "DEWPOINT": "Dew Point Temperature (F)",
    }

    mpm1s = np.array([[means+stdvs], [means-stdvs]]).T.reshape([N, 2])
    mpm1spaths = np.array([[meanpaths+stdvpaths], [meanpaths-stdvpaths]]).T.reshape([N, 2])
    # plotting starts here
    datplot = np.random.randint(0, Ndat, size = Ndatplot)
    plt.plot(datat[:,datplot], altitudes, '-', linewidth=1.0, color = '0.3', alpha = 0.2, solid_capstyle='projecting')
    plt.plot([], [], '-', linewidth=1.0,  label = 'Original Data')
    plt.plot(means, altitudes, 'k-', linewidth=1.6)
    plt.plot([], [], 'k-', linewidth=1.6, label = r'$\mu$')
    plt.plot(mpm1s, altitudes, 'k--', linewidth=1.6)
    plt.plot([], [], 'k--', linewidth=1.6, label = r'$\mu \pm 1\sigma$')
    plt.legend()
    plt.xlabel(pprop2label[prop])
    plt.ylabel('Altitude (1000 ft)')
    plt.savefig(f'{name}_{prop}_pathsinit.png', bbox_inches="tight", dpi=500)
    plt.clf()

    plt.plot(pathsgent, altitudes,  '-', linewidth=1.0)
    plt.plot([], [], '-', linewidth=1.0,  label = 'Synthetic Data')
    plt.plot(means, altitudes, 'k-', linewidth=1.6)
    plt.plot([], [], 'k-', linewidth=1.6, label = r'$\mu$ (Original)')
    plt.plot(mpm1s, altitudes, 'k--', linewidth=1.6)
    plt.plot([], [], 'k--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Original)')
    plt.plot(meanpaths, altitudes, 'b-', linewidth=1.6)
    plt.plot([], [], 'b-', linewidth=1.6, label = r'$\mu$ (Synthetic)')
    plt.plot(mpm1spaths, altitudes, 'b--', linewidth=1.6)
    plt.plot([], [], 'b--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Synthetic)')
    plt.xlabel(pprop2label[prop])
    plt.ylabel('Altitude (1000 ft)')
    plt.legend()
    plt.savefig(f'{name}_{prop}_t{trunc}_pathsgen.png', bbox_inches="tight", dpi=500)
    plt.clf()
    # now the mean, std, pdf errors
    # sig = 3
    # nsamp = 2000

    # kdetv = np.zeros([N, nsamp])
    # xs = np.zeros([N,nsamp])
    # for i in range(N):
    #     x = np.linspace(means[i]-sig*stdvs[i], means[i]+sig*stdvs[i], nsamp)
    #     kde = gaussian_kde(datat[i,:], bw_method='silverman')
    #     kdetv[i,:] = kde.pdf(x)
    #     xs[i,:] = x

    # # now get realization stats (MC ONLY)
    # meansm = np.mean(pathsgent, axis=1)
    # stdvsm = np.std(pathsgent, axis=1)

    # kdemv = np.zeros([N,nsamp])
    # for i in range(N):
    #     kde = gaussian_kde(pathsgent[i,:], bw_method='silverman')
    #     kdemv[i,:] = kde.pdf(xs[i,:])
    
    # means_err = abs(means-meansm)
    # stdvs_err = abs(stdvs-stdvsm)
    # pdfs_err = np.sum(abs(kdetv-kdemv), axis=1)
    # print('Means Error')
    # print(means_err)
    # print('Stdvs Error')
    # print(stdvs_err)
    # print('PDF Error')
    # print(pdfs_err)

    # import matplotlib as mpl
    # import matplotlib.pyplot as plt
    # colourmap = mpl.colormaps['rainbow']
    # # plt.figure().set_figheight(9.6)
    # plt.figure().set_figheight(4.8)
    # # plt.plot(meansm, altitudes, 'b-', linewidth=1.1)
    # # plt.plot(meansm + stdvsm, altitudes, 'b--', linewidth=1.1)
    # # plt.plot(meansm - stdvsm, altitudes, 'b--', linewidth=1.1)
    # altind = 4
    # for j in [altind]:#range(N):
    #     maxer = np.max(kdetv[j,:])
    #     plt.plot(xs[j,:], kdetv[j,:], 'k-', label='true pdf') #/maxer + altitudes[j]
    #     for i in range(nsamp - 1):
    #         # y2 = [altitudes[j], altitudes[j]]
    #         # y1 = [kdemv[j,i]/maxer+y2[0], kdemv[j,i+1]/maxer+y2[1]]
    #         y1 = [kdemv[j,i], kdemv[j,i+1]]
    #         plt.fill_between([xs[j,i], xs[j,i+1]],
    #                          y1,
    #                         #  y2,
    #                          color=colourmap(kdemv[j,i]/np.max(kdemv[3,:]))
    #                          ,alpha=0.6)
    # plt.title(f'{prop} PDF at {altitudes[altind]:.1f} thousand ft')
    # plt.ylabel(f'PDF({prop})')
    # plt.xlabel(prop)
    # plt.savefig(f'{name}_{prop}_t{trunc}_cases_PDF_comp.png', bbox_inches="tight", dpi=500)





# import pdb; pdb.set_trace()

# 1. preprocess_data
# 2. get_kl_coefficients
# 3. truncated_karhunen_loeve_expansion