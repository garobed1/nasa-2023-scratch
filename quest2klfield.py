import numpy as np
import os, sys
import argparse
from rand_field import preprocess_data, get_kl_coefficients, truncated_karhunen_loeve_expansion
from scipy.stats import gaussian_kde


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

# before anything, first estimate the covariance of the data
# and find its eigenpairs

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true') 
parser.add_argument('-p', '--plot', action='store_true') 
parser.add_argument('-d', '--datadir', action='store')
# if given, don't convert humidity to dew point temperature, create KL, then convert back 
# not recommended, better to convert to dew point using corresponding temp
parser.add_argument('-hd', '--donthumiditytodew', action='store_true')
# exclude a number of data points at the upper atmosphere 
# last 4 points constitute data above ~61000 feet
parser.add_argument('-e', '--exclude_upper', action='store', type=int, default=0)
# if true, measure mean, stdv, and PDF error (at each altitude)
parser.add_argument('-m', '--measure_error', action='store_true')
# MODEL HUMIDITY FROM TEMPERATURE/DEW POINT
# parser.add_argument('-H', '--modelhumidity', action='store_true')
parser.add_argument('-n', '--numpoints', default=0) # if 0, use mean of altitudes 
args = parser.parse_args()
verbose = args.verbose
Ngrid = args.numpoints
pflag = args.plot
datadir = args.datadir
exclude = args.exclude_upper
hd = args.donthumiditytodew
me = args.measure_error
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

    # EDF is Edwards AFB    
    altitudes, datat, means, stdvs, name = preprocess_data(data_file, prop, Ngrid)
    
    # if prop is humidity, convert to dew point, but get temp first
    if prop =='HUMIDITY' and not hd:
        # from Extract_statistics.py @author lwhite7 
        # hard clip negative humidities to minimum abs val
        minclip = np.min(abs(datat))
        for i in range(N):
            for j in range(Ndat):
                if datat[i,j] < 0.:
                    datat[i,j] = minclip
        #Get temp data back
        AT, DTT, MT, ST, NT = preprocess_data(data_file, 'TEMP', Ngrid)
        #Convert temp from Fahrenheit to Celsius
        temp_c = (DTT - 32)/1.8 
        #Calculate dew point using equation found here: https://en.wikipedia.org/wiki/Dew_point
        dew_pt = (257.14*(np.log(datat/100*np.exp((18.678-temp_c/234.5)*(temp_c/(257.14+temp_c))))))/(18.678-np.log(datat/100*np.exp((18.678-temp_c/234.5)*(temp_c/(257.14+temp_c)))))
        #convert dew point to Fahrenheit
        dew_pt = 1.8*dew_pt+32
        # Now sub in dew_pt in datat
        datat = dew_pt
        means = np.mean(datat, axis=1)
        stdvs = np.std(datat, axis=1)

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

    # get kl coefficients
    eigval, eigvec = get_kl_coefficients(datat, norm=False)

    if verbose:
        print(f"Generating {prop} profiles for all cases")

    casecounter = 0
    for case in os.listdir(case_dir):
        casecounter += 1

    if me:
        pathgent = np.zeros([N, casecounter])

    casestrlist = []
    casenumlist = []
    for case in os.listdir(case_dir):
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
            
            if pflag and trunc:
                if all(datag == 0.):
                    plt.plot(pathgen, altitudes, 'k-', linewidth=1.6, zorder=999999)
                else:
                    plt.plot(pathgen, altitudes,  '-', linewidth=1.0)

            # hold on to the path if we're measuring error
            if me:
                pathgent[:, casenum] = pathgen[:,0]

            # now convert dew point and temp back to humidity, and clip values over 100
            # actually, what temperature would we use to convert back?
            if prop =='HUMIDITY' and not hd:
                # from Extract_statistics.py @author lwhite7 
                # If we're modeling temperature as well, convert using the corresponding temp values
                if os.path.isfile(case_dir + '/' + case +'/TEMP_profile.txt'):
                    temp_k = np.zeros(N)
                    with open(case_dir + '/' + case + '/' + f'TEMP_profile.txt' , 'r') as wf:
                        pc = 0
                        for line in wf:
                            temp_f = float(line.split()[1])
                            temp_k[pc] = 5/9*(temp_f+459.67)
                            pc += 1

                # Otherwise, use the mean temperature MT, and write MT to every TEMP_profile.txt
                else:
                    with open(case_dir + '/' + case + '/' + f'TEMP_profile.txt' , 'w') as wf:
                        for i in range(N):
                            wf.write(f'{altitudes[i]} {MT[i]}\n')
                    temp_k = 5/9*(MT+459.67)
                #Convert pathsgen dew temp to humidity
                dew_pt_K = (5/9*(pathgen+459.67)).flatten() # convert to kelvin
                pathgen = 100*(6.112*np.exp(17.67*(dew_pt_K-273.15)/(dew_pt_K-29.65))/(6.112*np.exp(17.67*(temp_k-273.15)/(temp_k-29.65))))
                pathgen = np.atleast_2d(pathgen).T
            # if pflag and trunc:
            #     if all(datag == 0.):
            #         plt.plot(pathgen, altitudes, 'k-', linewidth=1.6, zorder=999999)
            #     else:
            #         plt.plot(pathgen, altitudes,  '-', linewidth=1.0)
            # write to file
            with open(case_dir + '/' + case + '/' + f'{prop}_profile.txt' , 'w') as wf:
                for i in range(N):
                    wf.write(f'{altitudes[i]} {pathgen[i][0]}\n')



    if pflag and trunc:
        pname = prop
        plt.plot([], [], '-', linewidth=1.0,  label = 'Synthetic Data (QUEST)')
        plt.plot([], [], 'k-', linewidth=1.6,  label = 'path at center index')
        if prop =='HUMIDITY' and not hd:
            pname = 'DEW POINT TEMP'
        plt.xlabel(pname)
        plt.ylabel('Altitude (1000 ft)')
        plt.legend()
        plt.savefig(f'{root}/{pname}_t{trunc}_{casecounter}_cases_pathsgen_QUEST.png', bbox_inches="tight", dpi=500)
        plt.clf()

    if me and trunc:
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

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        colourmap = mpl.colormaps['rainbow']
        # plt.figure().set_figheight(9.6)
        plt.figure().set_figheight(4.8)
        # plt.plot(meansm, altitudes, 'b-', linewidth=1.1)
        # plt.plot(meansm + stdvsm, altitudes, 'b--', linewidth=1.1)
        # plt.plot(meansm - stdvsm, altitudes, 'b--', linewidth=1.1)
        altind = 9
        for j in [altind]:#range(N):
            maxer = np.max(kdetv[j,:])
            plt.plot(xs[j,:], kdetv[j,:], 'k-', label='true pdf') #/maxer + altitudes[j]
            for i in range(nsamp - 1):
                # y2 = [altitudes[j], altitudes[j]]
                # y1 = [kdemv[j,i]/maxer+y2[0], kdemv[j,i+1]/maxer+y2[1]]
                y1 = [kdemv[j,i], kdemv[j,i+1]]
                plt.fill_between([xs[j,i], xs[j,i+1]],
                                 y1,
                                #  y2,
                                 color=colourmap(kdemv[j,i]/np.max(kdemv[3,:]))
                                 ,alpha=0.6)
        plt.title(f'{prop} PDF at {altitudes[altind]:.1f} thousand ft')
        plt.ylabel(f'PDF({prop})')
        plt.xlabel(prop)
        plt.savefig(f'{root}/{prop}_t{trunc}_{casecounter}_cases_PDF_comp.png', bbox_inches="tight", dpi=500)




# find out how many variables were given in QUEST

# naming convention
# KL_{prop}_{num}