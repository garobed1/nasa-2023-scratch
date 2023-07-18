import numpy as np
import os, sys
import json
import argparse
from rand_field import preprocess_data, get_kl_coefficients, truncated_karhunen_loeve_expansion
"""
Treat the original data as Monte Carlo realizations and write them to the 
case directories. This is done entirely outside of quest, so some care
will need to be taken. Get to the point where we can run c3d_parallel_runner.pl.
Options to either use all of the corresponding profiles simultaneously,
or one at a time (e.g. TEMP only, rest default).

Follow this convention: 
# Var_Name Var_Value
A0 T0
A1 T1

I think as long as the cases/case.??.????? directories ands templates with
mk_presb.csh exist, we can use quest2sboom.csh as a shortcut; just generate
the same text files as before
"""


# cmd line options
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true') 
# acceptable args are anything in proplist or 'ALL', which does all of them
parser.add_argument('-m', '--modelprop', action='store', default='TEMP')
# parser.add_argument('-p', '--plot', action='store_true') 
parser.add_argument('-d', '--datafile', action='store')

args = parser.parse_args()
verbose = args.verbose
proparg = args.modelprop
datafile = args.datafile

fulllist = ['TEMP', 'HUMIDITY', 'PRESSURE', 'WINDX', 'WINDY', 'TEMPHUMID']
if proparg == 'ALL':
    proplist = fulllist
elif proparg == 'TEMPHUMID':
    proplist = ['TEMP', 'HUMIDITY']
elif proparg not in fulllist:
    Exception("Valid property/option not given!")
else:
    proplist = [proparg]

# file naming
root = os.getcwd()
# quest_file = 'QUEST.dat'
data_file = 'fulldata.json' # assume its in the directory
if datafile is not None:
    data_file = datafile
case_dir = f'{root}/cases'


# first, generate the KL expansion
if verbose:
    print(f"Root dir: {root}")
    print(f"Preprocessing data ...")

# begin loop over properties
for prop in proplist:
    trunc = None # different properties may call for different numbers of vars

    # NOTE: Not processing or interpolating data here
    with open(data_file) as fj:
        fulldata = json.load(fj)

    datat = np.array(fulldata[prop]).T
    altitudes = np.array(fulldata['altitude']).T

    Ndat = datat.shape[1]
    N = datat.shape[0]

    # since we're not dealing with Quest here, make sure the cases dirs exist
    if not os.path.exists(case_dir):
        os.makedirs(case_dir)
    for i in range(Ndat):
        cnumstr = str(i)
        cnumstr = ((5-len(cnumstr))*'0') + cnumstr
        dirnamei = case_dir + '/case.00.' + cnumstr
        if not os.path.exists(dirnamei):
            os.makedirs(dirnamei)

    if verbose:
        print(f"Processing {prop}, {N} altitudes")

    casecounter = 0
    for case in os.listdir(case_dir):
        casecounter += 1
        if os.path.isdir(case_dir + '/' + case):
            if verbose:
                print(f"Processing {prop} profile for {case}")
            casenum = int(case.split('.')[-1])
            
            # write to file
            with open(case_dir + '/' + case + '/' + f'{prop}_profile.txt' , 'w') as wf:
                for i in range(N):
                    wf.write(f'{altitudes[i,casenum]} {datat[i,casenum]}\n')
            







# find out how many variables were given in QUEST

# naming convention
# KL_{prop}_{num}