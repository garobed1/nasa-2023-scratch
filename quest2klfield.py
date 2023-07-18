import numpy as np
import os, sys
import argparse
from rand_field import preprocess_data, get_kl_coefficients, truncated_karhunen_loeve_expansion
from humidity_transform import hum2dew, dew2hum
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

# write mean, sigma, m+1s, m-1s alongside each path (Quest likes this format)
parser.add_argument('-ws', '--writeoriginalstats', action='store_true')

# exclude a number of data points at the upper atmosphere 
# last 4 points constitute data above ~61000 feet
parser.add_argument('-e', '--exclude_upper', action='store', type=int, default=0)
# if true, measure mean, stdv, and PDF error (at each altitude)
parser.add_argument('-m', '--measure_error', action='store_true')

parser.add_argument('-n', '--numpoints', default=0) # if 0, use mean of altitudes as grid
args = parser.parse_args()
verbose = args.verbose
Ngrid = args.numpoints
pflag = args.plot
datadir = args.datadir
exclude = args.exclude_upper
hd = args.donthumiditytodew
me = args.measure_error
ws = args.writeoriginalstats
# mhflag = args.modelhumidity

trueproplist = ['TEMP', 'HUMIDITY', 'PRESSURE', 'WINDX', 'WINDY', 'DEWPOINT']
proplist = ['TEMPHUMID', 'TEMP', 'HUMIDITY', 'PRESSURE', 'WINDX', 'WINDY']
# proplist = ['TEMPHUMID']
# 'TEMPHUMID' simulates both TEMP and HUMIDITY/DEW POINT together, do not use
# with TEMP or HUMIDITY

root = os.getcwd()
quest_file = 'QUEST.dat'

data_file = 'fulldata.json'
if datadir is not None:
    data_file = datadir
case_dir = f'{root}/cases'
Ngen = len(os.listdir(case_dir))


if pflag:
    import matplotlib.pyplot as plt
    # store true, unconverted values in dict format
    ppathgent = {}
    d1,d2,d3,d4,d5 = preprocess_data(data_file, 'HUMIDITY', Ngrid)
    for item in trueproplist:
        ppathgent[item] = np.zeros([d2.shape[0]-exclude, Ngen])

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
    if prop == 'TEMPHUMID':
        altitudes, datat, means, stdvs, name = preprocess_data(data_file, 'HUMIDITY', Ngrid)    
    else:
        altitudes, datat, means, stdvs, name = preprocess_data(data_file, prop, Ngrid)

    # if prop is humidity, convert to dew point, but get temp first
    if ((prop == 'HUMIDITY' and not hd) or prop == 'TEMPHUMID') :
        # from Extract_statistics.py @author lwhite7 
        # hard clip negative humidities to minimum abs val
        N, Ndat = datat.shape
        minclip = np.min(abs(datat))
        for i in range(N):
            for j in range(Ndat):
                if datat[i,j] < 0.:
                    datat[i,j] = minclip
        #Get temp data back
        AT, DTT, MT, ST, NT = preprocess_data(data_file, 'TEMP', Ngrid)

        #Get dew point temps
        datat = hum2dew(datat, DTT)

        means = np.mean(datat, axis=1)
        stdvs = np.std(datat, axis=1)
        
        if 0:
            mpm1s = np.array([[means+stdvs], [means-stdvs]]).T.reshape([N, 2])
            plt.plot(datat, altitudes, '-', linewidth=1.0)
            plt.plot([], [], '-', linewidth=1.0,  label = 'Original Data')
            plt.plot(means, altitudes, 'k-', linewidth=1.6)
            plt.plot([], [], 'k-', linewidth=1.6, label = r'$\mu$')
            plt.plot(mpm1s, altitudes, 'k--', linewidth=1.6)
            plt.plot([], [], 'k--', linewidth=1.6, label = r'$\mu \pm 1\sigma$')
            plt.legend()
            plt.xlabel('DEWPT')
            plt.ylabel('Altitude (1000 ft)')
            plt.savefig(f'{name}_DEWPT_pathsinit.png', bbox_inches="tight", dpi=500)
            plt.clf()
            import pdb; pdb.set_trace()

    # if exclude > 0, remove upper indices
    if exclude > 0:
        altitudes = altitudes[:-exclude]
        datat = datat[:-exclude, :]
        means = means[:-exclude]
        stdvs = stdvs[:-exclude]
        if prop == 'TEMPHUMID':
            AT = AT[:-exclude]
            DTT = DTT[:-exclude, :]
            MT = MT[:-exclude]
            ST = ST[:-exclude]


    # if TEMPHUMID, concatenate temp and dew point data
    if prop == 'TEMPHUMID':
        altitudes = np.append(AT, altitudes, axis=0)
        datat = np.append(DTT, datat, axis=0)
        means = np.append(MT, means, axis=0)
        stdvs = np.append(ST, stdvs, axis=0)

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

    # if me:
    pathgent = np.zeros([N, casecounter]) # these include dew conversions

    casestrlist = []
    casenumlist = []

    badcasecounter = 0
    for case in os.listdir(case_dir):
        if os.path.isdir(case_dir + '/' + case):
            # if verbose:
            #     print(f"Generating {prop} profile for {case}")
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
            
            # if pflag and trunc:
            #     if all(datag == 0.):
            #         plt.plot(pathgen, altitudes, 'k-', linewidth=1.6, zorder=999999)
            #     else:
            #         plt.plot(pathgen, altitudes,  '-', linewidth=1.0)


            # store raw path
            pathgent[:, casenum] = pathgen[:,0]

            # now convert dew point and temp back to humidity, and clip values over 100
            # actually, what temperature would we use to convert back?
            if (prop =='HUMIDITY' and not hd) or prop == 'TEMPHUMID':
                # from Extract_statistics.py @author lwhite7 

                # If combining temp and humidity, get temp from the data we just got
                # Write temp here immediately, don't do it for the other cases
                if prop == 'TEMPHUMID':
                    split_ind = int(N/2)
                    temp_use = pathgen[:split_ind]
                    dew_pt_use = pathgen[split_ind:]
                    with open(case_dir + '/' + case + '/' + f'TEMP_profile.txt' , 'w') as wf:
                        for i in range(split_ind):
                            wf.write(f'{altitudes[i]} {temp_use[i]}\n')
                    ppathgent['TEMP'][:, casenum] = temp_use[:,0]

                # If we're modeling temperature as well, convert using the corresponding temp values
                elif os.path.isfile(case_dir + '/' + case +'/TEMP_profile.txt'):
                    # temp_k = np.zeros(N)
                    temp_use = np.zeros(N)
                    with open(case_dir + '/' + case + '/' + f'TEMP_profile.txt' , 'r') as wf:
                        pc = 0
                        for line in wf:
                            temp_use[pc] = float(line.split()[1])
                            pc += 1
                    dew_pt_use = pathgen
                
                # Otherwise, use the mean temperature MT, and write MT to every TEMP_profile.txt
                else:
                    temp_use = MT # mean + 1 sigma seems to work better e.g. higher temps work better 
                    with open(case_dir + '/' + case + '/' + f'TEMP_profile.txt' , 'w') as wf:
                        for i in range(N):
                            wf.write(f'{altitudes[i]} {temp_use[i]}\n')
                    dew_pt_use = pathgen
                
                ppathgent['DEWPOINT'][:, casenum] = dew_pt_use[:,0]
                
                #Convert pathsgen dew temp to humidity
                # NOTE: This function ensures dew point is always 
                # less than dry bulb temp
                if any((dew_pt_use - temp_use)>0) :
                    badcasecounter += 1

                work = dew2hum(dew_pt_use, temp_use)
                
                # if any(work > 100.):
                #     import pdb; pdb.set_trace()
                ppathgent['HUMIDITY'][:, casenum] = work[:,0]

                if prop != 'TEMPHUMID':
                    pathgen = work
                else:
                    pathgen = np.append(temp_use, work)
                    # import pdb; pdb.set_trace()
                    with open(case_dir + '/' + case + '/' + f'HUMIDITY_profile.txt' , 'w') as wf:
                        for i in range(split_ind):
                            wf.write(f'{altitudes[i]} {work[i]}\n')

            # write to file
            if not ((prop =='HUMIDITY' and not hd) or prop == 'TEMPHUMID'):
                ppathgent[prop][:, casenum] = pathgen[:,0]
                with open(case_dir + '/' + case + '/' + f'{prop}_profile.txt' , 'w') as wf:
                    for i in range(N):
                        wf.write(f'{altitudes[i]} {pathgen[i][0]}\n')


    if verbose:
        print(f"Cases where dew point exceeded temperature: {badcasecounter}")

    if pflag and trunc:
        pname = prop

        pproplist = [prop]
        if prop == 'TEMPHUMID':
            pproplist = ['TEMP', 'HUMIDITY', 'DEWPOINT']

            from matplotlib import cm


            # scatter plot
            cmap = cm.coolwarm
            altind = 10
            temp_alt = ppathgent['TEMP'][altind,:200]
            hum_alt = ppathgent['HUMIDITY'][altind,:200]
            # mat = np.cov(temp_alt, hum_alt)
            # import pdb; pdb.set_trace()
            plt.scatter(temp_alt, hum_alt)
            plt.xlim(20., 110.)
            plt.ylim(0., 100.)
            plt.title(f"Temp and Humidity correlation at {altitudes[altind]:.2f} thousand feet")
            plt.xlabel('Temp (F)')
            plt.ylabel('Humidity (%)')
            plt.savefig(f"{root}/{prop}_t{trunc}_{casecounter}_temp_hum_corr.png", bbox_inches='tight')
            plt.clf()


        for pprop in pproplist:
            #retrieve data in case we lost it somehow
            if pprop == 'DEWPOINT':
                AH, DTH, MH, SH, NH = preprocess_data(data_file, 'HUMIDITY', Ngrid)    
                Nh, Ndath = DTH.shape
                minclip = np.min(abs(DTH))
                for i in range(Nh):
                    for j in range(Ndath):
                        if DTH[i,j] < 0.:
                            DTH[i,j] = minclip
                AT, DTT, MT, ST, NT = preprocess_data(data_file, 'TEMP', Ngrid)


                pdew_pt = hum2dew(DTH, DTT)
                # Now sub in dew_pt in datat
                paltitudes = AH
                pdatat = pdew_pt
                pmeans = np.mean(pdatat, axis=1)
                pstdvs = np.std(pdatat, axis=1)

                if 0:
                    # just a test
                    pdew_pt_M = hum2dew(MH, MT)
                    work = dew2hum(pmeans, MT)
                    import pdb; pdb.set_trace()
            else:
                paltitudes, pdatat, pmeans, pstdvs, pname = preprocess_data(data_file, pprop, Ngrid)
            
            if exclude > 0:
                paltitudes = paltitudes[:-exclude]
                pdatat = pdatat[:-exclude, :]
                pmeans = pmeans[:-exclude]
                pstdvs = pstdvs[:-exclude]
            
            Ntrue = pdatat.shape[0]

            pmeanpaths = np.mean(ppathgent[pprop], axis=1)
            pstdvpaths = np.std(ppathgent[pprop], axis=1)

            mpm1s = np.array([[pmeans+pstdvs], [pmeans-pstdvs]]).T.reshape([Ntrue, 2])
            mpm1spaths = np.array([[pmeanpaths+pstdvpaths], [pmeanpaths-pstdvpaths]]).T.reshape([Ntrue, 2])
            plt.plot(ppathgent[pprop], paltitudes,  '-', linewidth=1.0)
            plt.plot([], [], '-', linewidth=1.0,  label = 'Synthetic Data (QUEST)')
            # plt.plot(pdatat, paltitudes,  '-', linewidth=1.0)
            # plt.plot([], [], '-', linewidth=1.0,  label = 'Original Data')
            # import pdb; pdb.set_trace()
            plt.plot(pmeans, paltitudes, 'k-', linewidth=1.6)
            plt.plot([], [], 'k-', linewidth=1.6, label = r'$\mu$ (Original)')
            plt.plot(mpm1s, paltitudes, 'k--', linewidth=1.6)
            plt.plot([], [], 'k--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Original)')
            plt.plot(pmeanpaths, paltitudes, 'b-', linewidth=1.6)
            plt.plot([], [], 'b-', linewidth=1.6, label = r'$\mu$ (Synthetic)')
            plt.plot(mpm1spaths, paltitudes, 'b--', linewidth=1.6)
            plt.plot([], [], 'b--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Synthetic)')
            
            # plt.plot([], [], 'k-', linewidth=1.6,  label = 'path at center index')
            # if prop =='HUMIDITY' and not hd:
            #     pname = 'DEW_POINT_TEMP'
            plt.xlabel(pprop)
            plt.ylabel('Altitude (1000 ft)')
            plt.legend()
            plt.savefig(f'{root}/{pprop}_t{trunc}_{casecounter}_cases_pathsgen_QUEST.png', bbox_inches="tight", dpi=500)
            plt.clf()

            # plt.plot(ppathgent[pprop], paltitudes,  '-', linewidth=1.0)
            # plt.plot([], [], '-', linewidth=1.0,  label = 'Synthetic Data (QUEST)')
            plt.plot(pdatat, paltitudes,  '-', linewidth=1.0)
            plt.plot([], [], '-', linewidth=1.0,  label = 'Original Data')
            # import pdb; pdb.set_trace()
            plt.plot(pmeans, paltitudes, 'k-', linewidth=1.6)
            plt.plot([], [], 'k-', linewidth=1.6, label = r'$\mu$ (Original)')
            plt.plot(mpm1s, paltitudes, 'k--', linewidth=1.6)
            plt.plot([], [], 'k--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Original)')
            plt.plot(pmeanpaths, paltitudes, 'b-', linewidth=1.6)
            plt.plot([], [], 'b-', linewidth=1.6, label = r'$\mu$ (Synthetic)')
            plt.plot(mpm1spaths, paltitudes, 'b--', linewidth=1.6)
            plt.plot([], [], 'b--', linewidth=1.6, label = r'$\mu \pm 1\sigma$ (Synthetic)')
            
            # plt.plot([], [], 'k-', linewidth=1.6,  label = 'path at center index')
            # if prop =='HUMIDITY' and not hd:
            #     pname = 'DEW_POINT_TEMP'
            plt.xlabel(pprop)
            plt.ylabel('Altitude (1000 ft)')
            plt.legend()
            plt.savefig(f'{root}/{pprop}_t{trunc}_{casecounter}_cases_DATA_QUEST.png', bbox_inches="tight", dpi=500)
            plt.clf()

            # write pmeans and mpm1s for later Quest stuff
            if ws:
                with open(root + '/' + f'{pprop}MEAN_profile.txt' , 'w') as wf:
                    wf.write(f'# {Ntrue}\n')
                    wf.write('# mean$(Original) Altitude$(1000$ft) Temperature$(F)\n')
                    for i in range(Ntrue):
                        wf.write(f'{paltitudes[i]} {pmeans[i]}\n')
                with open(root + '/' + f'{pprop}SIGMA_profile.txt' , 'w') as wf:
                    wf.write(f'# {Ntrue}\n')
                    wf.write('# σ$(Original) Altitude$(1000$ft) Temperature$(F)\n')
                    for i in range(Ntrue):
                        wf.write(f'{paltitudes[i]} {pstdvs[i]}\n')
                with open(root + '/' + f'{pprop}MP1S_profile.txt' , 'w') as wf:
                    wf.write(f'# {Ntrue}\n')
                    wf.write('# mean+1σ$(Original) Altitude$(1000$ft) Temperature$(F)\n')
                    for i in range(Ntrue):
                        wf.write(f'{paltitudes[i]} {mpm1s[i,0]}\n')
                with open(root + '/' + f'{pprop}MM1S_profile.txt' , 'w') as wf:
                    wf.write(f'# {Ntrue}\n')
                    wf.write('# mean-1σ$(Original) Altitude$(1000$ft) Temperature$(F)\n')
                    for i in range(Ntrue):
                        wf.write(f'{paltitudes[i]} {mpm1s[i,1]}\n')


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