import numpy as np
import os, sys
import argparse
from rand_field import preprocess_data, get_kl_coefficients, truncated_karhunen_loeve_expansion
from scipy.stats import gaussian_kde

"""
Compare statistics of KL fields to that of original data after processing
through QUEST



"""



        # mean l2 error
        # means is true
        # stdvs is true
        # get true pdf
        sig = 3
        nsamp = 2000

        kdetv = np.zeros([N, nsamp])
        xs = np.zeros([N,nsamp])
        for i in range(N):
            x = np.linspace(means[i]-sig*stdvs[i], means[i]+sig*stdvs[i], nsamp)
            kde = gaussian_kde(datat[i,:], bw_method='silverman')
            kdetv[i,:] = kde.pdf(x)
            xs[i,:] = x

        # now get realization stats (MC ONLY)
        meansm = np.mean(pathgent, axis=1)
        stdvsm = np.std(pathgent, axis=1)

        kdemv = np.zeros([N,nsamp])
        for i in range(N):
            kde = gaussian_kde(pathgent[i,:], bw_method='silverman')
            kdemv[i,:] = kde.pdf(xs[i,:])
        
        means_err = abs(means-meansm)
        stdvs_err = abs(stdvs-stdvsm)
        pdfs_err = np.sum(abs(kdetv-kdemv), axis=1)
        print('Means Error')
        print(means_err)
        print('Stdvs Error')
        print(stdvs_err)
        print('PDF Error')
        print(pdfs_err)

        # write to file
        with open(root + f'/{prop}_stat_errs.txt' , 'w') as wf:
                for i in range(N):
                    wf.write(f'{altitudes[i]} {means_err[i]} {stdvs_err[i]} {pdfs_err[i]}\n')
