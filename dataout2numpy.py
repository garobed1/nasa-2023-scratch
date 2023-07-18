import os, sys
import numpy as np
import json
import argparse

"""
Convert file I/O sboom output data from .dat files to numpy arrays

For the sake of data_post.py
"""

# cmd line options
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true') 
# read in signatures as well. # NOTE: read_signatures not implemented!
parser.add_argument('-s', '--readsigs', action='store_true')

args = parser.parse_args()
verbose = args.verbose
readsigs = args.readsigs

# file naming
root = os.getcwd()    
casedir = root + '/cases'
if verbose:
    print(f"Root dir: {root}")
# count number of cases
attcounter = 0
for case in os.listdir(casedir):
    attcounter += 1

print(f'{attcounter} outputs')

# find the properties modeled
fulllist = ['TEMP', 'HUMIDITY', 'PRESSURE', 'WINDX', 'WINDY']
proplist = []
for file in os.listdir(casedir + '/case.00.00000'):
    namecand = file.split('_')[0]
    if namecand in fulllist:
        proplist.append(namecand)

# find the loudness properties
fulldata = {}
fulldata['n'] = attcounter
fulldata['Carpet Width (km)'] = np.zeros(attcounter)
fulldata['Negative Side (km)'] = np.zeros(attcounter)
fulldata['Negative Side Angle (deg.)'] = np.zeros(attcounter)
fulldata['Positive Side (km)'] = np.zeros(attcounter)
fulldata['Positive Side Angle (deg.)'] = np.zeros(attcounter)
fulldata['Objective'] = np.zeros(attcounter)
fulldata['PL'] = np.zeros(attcounter)
fulldata['ISBAP'] = np.zeros(attcounter)
fulldata['SELa[0]'] = np.zeros(attcounter)
fulldata['SELb'] = np.zeros(attcounter)
fulldata['SELc'] = np.zeros(attcounter)
fulldata['SELd'] = np.zeros(attcounter)
fulldata['SELe'] = np.zeros(attcounter)
fulldata['Propagation Time (s)'] = np.zeros(attcounter)

# begin looping through cases
gfc = 0
for case in os.listdir(casedir):
    # get case number
    cnumstr = case.split('.')[-1]
    casenum = int(cnumstr)

    # go through the relevant output files
    with open(root + '/cases/' + case + '/sboom/loud.dat') as otf:
        for line in otf:
            if line.strip():
                qoi = line.split()[0]
                if qoi in fulldata.keys():
                    if qoi == 'Objective':
                        val = float(line.split()[4])
                    else:
                        val = float(line.split()[2])
                    fulldata[qoi][casenum] = val

    with open(root + '/cases/' + case + '/sboom/output.out') as otf:
        for line in otf:
            if line.strip():
                qoi = line.split(' ')[0]
                if qoi == 'Propagation':
                    fulldata['Propagation Time (s)'][casenum] = float(line.split()[2])

        # NOTE: Carpet width excluded here, will need to incorporate
    

fulldata_list = {}
for k,v in fulldata.items():
    if isinstance(v, np.ndarray):
        fulldata_list[k] = v.tolist() 
    else:
        fulldata_list[k] = v
with open(f'{root}/{proplist}_data_loud.json', 'w') as fj:
    json.dump(fulldata_list, fj)




     