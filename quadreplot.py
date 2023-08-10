#!/usr/bin/env python3

import numpy as np
import os, sys
import argparse
import matplotlib.pyplot as plt
    
from rand_field import preprocess_data

plt.rcParams['font.size'] = 16

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--datadir', action='store')
parser.add_argument('-e', '--exclude_upper', action='store', type=int, default=0)

args = parser.parse_args()
datadir = args.datadir
exclude = args.exclude_upper

proplist = ['TEMP', 'HUMIDITY', 'PRESSURE', 'WINDX', 'WINDY']

pprop2label = {
    "TEMP": "Temperature (F)",
    "HUMIDITY": "Relative Humidity (%)",
    "DEWPOINT": "Dew Point Temperature (F)",
}

root = os.getcwd()
data_file = 'fulldata.json'
if datadir is not None:
    data_file = datadir
case_dir = f'{root}/cases'
Ngen = len(os.listdir(case_dir))

for prop in proplist:
    # check if property exists
    try:
        cur_file = case_dir + '/case.00.00000/' + prop + '_profile.txt'
        with open(cur_file, 'r') as cf:
            pass
    except:
        continue

    altitudes, datat, means, stdvs, name = preprocess_data(data_file, prop, None)

    if exclude > 0:
        altitudes = altitudes[:-exclude]
        datat = datat[:-exclude, :]
        means = means[:-exclude]
        stdvs = stdvs[:-exclude]

    Ndat = datat.shape[1]
    N = datat.shape[0]

    # grab profiles from files
    case_list = os.listdir(case_dir)
    case_list = [x for x in case_list if not x.startswith('CFD_')]
    casecounter = len(case_list)
    datap = np.zeros([N, casecounter])

    for case in case_list:
        if os.path.isdir(case_dir + '/' + case):
            casenum = int(case.split('.')[-1])
            cur_file = case_dir + '/' + case + '/' + prop + '_profile.txt'
            with open(cur_file, 'r') as cf:
                n = 0
                for line in cf:
                    data = line.split(' ')
                    val = float(data[1])
                    datap[n,casenum] = val
                    n += 1
    # get mean and mean +- 1sigma from QUEST output
    wmp = []
    with open(f'{root}/KL{prop.upper()}.Altitude (1000 ft).mean') as f:
        for line in f:
            if not line.startswith('#'):
                wmp.append(float(line.split()[1]))
    wmpsp = []
    with open(f'{root}/KL{prop.upper()}.Altitude (1000 ft).meansigma_p') as f:
        for line in f:
            if not line.startswith('#'):
                wmpsp.append(float(line.split()[1]))
    wmpsm = []
    with open(f'{root}/KL{prop.upper()}.Altitude (1000 ft).meansigma_m') as f:
        for line in f:
            if not line.startswith('#'):
                wmpsm.append(float(line.split()[1]))

    # plot
    mpm1s = np.array([[means+stdvs], [means-stdvs]]).T.reshape([N, 2])
    mpm1sgen = np.array([wmpsp, wmpsm]).T.reshape([N,2])
    plt.plot(datap, altitudes,  '-', linewidth=1.3, color = '0.3', alpha = 0.2, solid_capstyle='projecting')
    plt.plot([], [], '-', linewidth=1.0,  label = 'Synthetic Data (QUEST)')
    # plt.plot(pdatat, paltitudes,  '-', linewidth=1.0)
    # plt.plot([], [], '-', linewidth=1.0,  label = 'Original Data')
    # import pdb; pdb.set_trace()
    plt.plot(means, altitudes, 'k-', linewidth=1.6)
    plt.plot([], [], 'k-', linewidth=1.6, label = r'$\mu$ (Original)')
    plt.plot(mpm1s, altitudes, 'k--', linewidth=1.6)
    plt.plot([], [], 'k--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Original)')
    plt.plot(wmp, altitudes, 'b-', linewidth=1.6)
    plt.plot([], [], 'b-', linewidth=1.6, label = r'$\mu$ (Synthetic)')
    plt.plot(mpm1sgen, altitudes, 'b--', linewidth=1.6)
    plt.plot([], [], 'b--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Synthetic)')
    
    # plt.plot([], [], 'k-', linewidth=1.6,  label = 'path at center index')
    # if prop =='HUMIDITY' and not hd:
    #     pname = 'DEW_POINT_TEMP'
    plt.xlabel(pprop2label[prop])
    plt.ylabel('Altitude (1000 ft)')
    plt.legend(fontsize=13)
    plt.savefig(f'{root}/{prop}_replot_{casecounter}_cases_pathsgen_QUEST.png', bbox_inches="tight", dpi=500)
    plt.clf()