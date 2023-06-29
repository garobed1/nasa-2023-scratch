"""
Generate PDF plots, mean, standard deviation of direct atmospheric data outputs

PDF generated from kernel density estimation (NOTE: Using Default Bandwidth est)
Same as basic Quest, see 5.15

Also unlike Quest, no weights, assuming 'Monte Carlo' distribution of data
"""
import os, sys
import numpy as np
import json
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import gaussian_kde

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true') 
parser.add_argument('-d', '--datafile', action='store', default='[\'TEMP\']_data_loud.json')
colourmap = mpl.colormaps['rainbow']

args = parser.parse_args()
verbose = args.verbose
datafile = args.datafile
root = os.getcwd()    

with open(datafile) as fj:
    fulldata = json.load(fj)

props = datafile.split('/')[-1]
props = props.split('_')[0]

data_stats = {}
data_pdfs = {}
for key, value in fulldata.items():
    if key == 'n':
        data_stats['n'] = value
        continue

    stats = np.zeros(2)
    stats[0] = np.mean(value)
    stats[1] = np.std(value)

    data_stats[key] = stats

    # kernel density estimation (NOTE: Using Default Bandwidth est)
    if all(value) == 0.0:
        continue
    kde = gaussian_kde(value, bw_method='silverman') # No weights, presuming monte carlo estimates
 
    # estimate over \pm 3 sigma
    sig = 3
    nsamp = 1000
    x = np.linspace(stats[0]-sig*stats[1], stats[0]+sig*stats[1], nsamp)
    y = kde.pdf(x)
    data_pdfs[key] = (x, y)

# just print out mu and sigma for each output
if verbose:
    print('Quantity    | Mean         | Std. Dev.      ')
    for key, value in fulldata.items():
        if key == 'n':
            continue

        stats = data_stats[key]
        print(key + f'  {stats[0]}  {stats[1]}')

        if all(value) == 0.0:
            continue

        x = data_pdfs[key][0]
        y = data_pdfs[key][1]
        y = y/np.max(y)
        plt.grid()
        plt.plot(x, y)
        # fill with color gradient
        for i in range(nsamp - 1):
            plt.fill_between([x[i], x[i+1]],
                             [y[i], y[i+1]],
                             color=colourmap(y[i])
                             ,alpha=0.6)
        plt.xlabel(key)
        plt.ylabel(f'PDF({key}) (normalized)')
        plt.savefig(f'{key}_{props}_pdf_est.png', bbox_inches='tight')
        plt.clf()