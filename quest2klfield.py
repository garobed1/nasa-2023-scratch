import numpy as np
import os, sys
import argparse
from rand_field import preprocess_data, get_kl_coefficients, truncated_karhunen_loeve_expansion
"""
Quest will output realizations of the KL coefficients to a database 
file. Put those realizations through the rand_field.py code to produce
characteristic paths. This script converts database values to a form 
readable as arguments to rand_field.py. The actual path itself depends
on the originating data, which is completely separate from Quest.
"""

"""
Follow this convention: 
# Var_Name Var_Value
TEMP_path_0 T0
TEMP_path_1 T1

HUMIDITY_path_0 H0
HUMIDITY_path_1 H1
"""

# start with one property for now

# before anything, first estimate the covariance of the data
# and find its eigenpairs

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true') 
parser.add_argument('-p', '--plot', action='store_true') 
parser.add_argument('-d', '--datadir', action='store')
# MODEL HUMIDITY FROM TEMPERATURE/DEW POINT
# parser.add_argument('-H', '--modelhumidity', action='store_true')
parser.add_argument('-n', '--numpoints', default=0) # if 0, use mean of altitudes 
args = parser.parse_args()
verbose = args.verbose
Ngrid = args.numpoints
pflag = args.plot
datadir = args.datadir
# mhflag = args.modelhumidity

proplist = ['TEMP', 'HUMIDITY', 'PRESSURE', 'WINDX', 'WINDY']

root = os.getcwd()
quest_file = 'QUEST.dat'

data_file = 'fulldata.json'
if datadir is not None:
    data_file = datadir
case_dir = f'{root}/cases'

if pflag:
    import matplotlib.pyplot as plt

if verbose:
    print(f"Root dir: {root}")
    # print(f" o Working in {db}")

# first, generate the KL expansion
if verbose:
    print(f"Preprocessing data ...")

# begin loop over properties
for prop in proplist:
    trunc = None # different properties may call for different numbers of vars

    
    altitudes, datat, means, stdvs, name = preprocess_data(data_file, prop, Ngrid)
    Ndat = datat.shape[1]
    N = datat.shape[0]

    if verbose:
        print(f"Processing {prop}, {N} altitudes")

    # get kl coefficients
    eigval, eigvec = get_kl_coefficients(datat, norm=False)

    if verbose:
        print(f"Generating {prop} profiles for all cases")

    casestrlist = []
    casenumlist = []
    casecounter = 0
    for case in os.listdir(case_dir):
        casecounter += 1
        if os.path.isdir(case_dir + '/' + case):
            if verbose:
                print(f"Generating {prop} profile for {case}")
            casenum = int(case.split('.')[-1])
            cur_file = case_dir + '/' + case + '/' + quest_file

            datag = None
            with open(cur_file) as cf:
                # useful flags
                lc = 0
                newvar = False
                propname = None
                # begin line loop to count number of variables
                if trunc is None:
                    for line in cf:
                        if line.startswith(f'KL_{prop}'): # only touch relevant lines
                            lc += 1
                        else:
                            pass
                    trunc = lc

                if trunc == 0:
                    continue

                datag = np.zeros([trunc, 1])
                for line in cf:
                    if line.startswith(f'KL_{prop}'): # only touch relevant lines
                        data = line.split(' ')
                        ind = int(data[0].split('_')[-1])
                        val = float(data[1])
                        datag[ind] = val
                    else:
                        pass
            
            # now produce the path
            pathgen = truncated_karhunen_loeve_expansion(datag, means, eigval, eigvec, prop)
            
            # write to file
            with open(case_dir + '/' + case + '/' + f'{prop}_profile.txt' , 'w') as wf:
                for i in range(N):
                    wf.write(f'{altitudes[i]} {pathgen[i][0]}\n')
            
            if pflag and trunc:
                if all(datag == 0.):
                    plt.plot(pathgen, altitudes, 'k-', linewidth=1.6, zorder=999999)
                else:
                    plt.plot(pathgen, altitudes,  '-', linewidth=1.0)

    if pflag and trunc:
        plt.plot([], [], '-', linewidth=1.0,  label = 'Synthetic Data (QUEST)')
        plt.plot([], [], 'k-', linewidth=1.6,  label = 'path at center index')
        plt.xlabel(prop)
        plt.ylabel('Altitude (1000 ft)')
        plt.legend()
        plt.savefig(f'{root}/{prop}_t{trunc}_{casecounter}_cases_pathsgen_QUEST.png', bbox_inches="tight", dpi=500)
        plt.clf()




# find out how many variables were given in QUEST

# naming convention
# KL_{prop}_{num}