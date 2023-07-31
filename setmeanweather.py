#!/usr/bin/env python3

import numpy as np
import os, sys
import argparse
from rand_field import preprocess_data
from humidity_transform import hum2dew, dew2hum

"""
If there is no weather uncertainty being considered, use the 
mean values of the given data for every case instead.
"""

"""
Follow this convention: 
# Var_Name Var_Value
TEMP_path_0 T0
TEMP_path_1 T1

HUMIDITY_path_0 H0
HUMIDITY_path_1 H1
"""

# before anything, first estimate the covariance of the data
# and find its eigenpairs

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true') 
parser.add_argument('-d', '--datadir', action='store')
# if given, don't convert humidity to dew point temperature, create KL, then convert back 
# not recommended, better to convert to dew point using corresponding temp
parser.add_argument('-hd', '--donthumiditytodew', action='store_true')

# exclude a number of data points at the upper atmosphere 
# last 4 points constitute data above ~61000 feet
parser.add_argument('-e', '--exclude_upper', action='store', type=int, default=0)

parser.add_argument('-n', '--numpoints', default=0) # if 0, use mean of altitudes as grid

args = parser.parse_args()
verbose = args.verbose
Ngrid = args.numpoints
datadir = args.datadir
exclude = args.exclude_upper
hd = args.donthumiditytodew
ws = args.writeoriginalstats
# mhflag = args.modelhumidity

proplist = ['TEMP', 'HUMIDITY']

root = os.getcwd()
quest_file = 'QUEST.dat'

data_file = 'fulldata.json'
if datadir is not None:
    data_file = datadir
case_dir = f'{root}/cases'
Ngen = len(os.listdir(case_dir))


if verbose:
    print(f"Root dir: {root}")
    # print(f" o Working in {db}")

# first, generate the KL expansion
if verbose:
    print(f"Preprocessing data (MEAN ONLY) ...")

# begin loop over properties
for prop in proplist:

    altitudes, datat, means, stdvs, name = preprocess_data(data_file, prop, Ngrid)

    # if exclude > 0, remove upper indices
    if exclude > 0:
        altitudes = altitudes[:-exclude]
        datat = datat[:-exclude, :]
        means = means[:-exclude]
        stdvs = stdvs[:-exclude]

    Ndat = datat.shape[1]
    N = datat.shape[0]

    if verbose:
        print(f"Processing {prop}, {N} altitudes")


    if verbose:
        print(f"Generating {prop} (MEAN) profiles for all cases")

    casecounter = 0
    for case in os.listdir(case_dir):
        casecounter += 1

    casestrlist = []
    casenumlist = []

    case_list = os.listdir(case_dir)
    case_list = [x for x in case_list if not x.startswith('CFD_')]
    for case in case_list:
        if os.path.isdir(case_dir + '/' + case):
            with open(case_dir + '/' + case + '/' + f'{prop}_profile.txt' , 'w') as wf:
                for i in range(N):
                    wf.write(f'{altitudes[i]} {means[i]}\n')

