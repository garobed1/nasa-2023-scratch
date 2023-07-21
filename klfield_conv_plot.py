import numpy as np
import matplotlib.pyplot as plt
import os, sys

"""
After running quest_post.csh and outputting mean/mean+1sigma into stats/ for each project

Give directories with _kl?l* as arguments (? is a specific number, * is an actual wildcard)

This one handles statistics errors of the input KL fields themselves

Can also potentially handle the combined temp/dewpttemp vector if given
NOTE: We show that dew point temp fields converge, but converted humidity does not

Ex.
cd ~/temphumid_lhs
python3 ~/nasa-2023-scratch/klfield_conv_plot.py sboom_atm_uq_temphumid_lhs_kl5l*
"""
# KL field measure to check for
# NOTE: Looking specifically for KLTEMP, KLHUMIDITY, etc.

# UQ method names
qtypenamedict = {
    'dcc':'Dense CC',
    'spcc':'Sparse CC',
    'dgp':'Dense GP',
    'spgp':'Sparse GP',
    'mc':'Monte Carlo',
    'lhs':'LHS'
}
qtypelist = qtypenamedict.keys()

# file naming
datdict = {}
caselist = sys.argv[1:]
root = os.getcwd()
proplist = ['temp', 'humidity', 'windx', 'windy']

# loop over each project
for case in caselist:

    # retrieve properties from directory names
    name = case.split('_')
    kl = name[-1]
    qtype = name[-2]
    pcat = name[-3]
    klp = kl[2] # KL param number
    kll = kl[4] # grid level

    if qtype not in qtypelist:
        qtype = 'dcc'
        pcat = name[-2]


    prop = [pcat]
    if pcat not in proplist:
        if pcat == 'temphumid':
            props = ['temp', 'humidity']
            # props = ['humidity']
    for prop in props:
        ft = klp + qtype + prop
        if ft not in datdict:
            datdict[ft] = {}
            datdict[ft]['qtype'] = qtypenamedict[qtype]
            datdict[ft]['KL'] = klp
            datdict[ft]['prop'] = prop
            datdict[ft]['s'] = []
            datdict[ft]['mean'] = []
            datdict[ft]['sigma'] = []
            datdict[ft]['pdf'] = []

        # get mean
        # get stdv
        wm = []
        ws = []
        wp = []
        try:
            with open(f'{root}/{case}/{prop.upper()}_stat_errs.txt') as f:
                for line in f:
                    wm.append(float(line.split()[1]))
                    ws.append(float(line.split()[2]))
                    wp.append(float(line.split()[3]))
        except:
            # if it doesn't exist, then take the quest post output instead
            # then compute it here
            wme = []
            with open(f'{root}/{case}/{prop.upper()}MEAN_profile.txt') as f:
                for line in f:
                    if not line.startswith('#'):
                        wme.append(float(line.split()[1]))
            wse = []
            with open(f'{root}/{case}/{prop.upper()}SIGMA_profile.txt') as f:
                for line in f:
                    if not line.startswith('#'):
                        wse.append(float(line.split()[1]))
            wmp = []
            with open(f'{root}/{case}/KL{prop.upper()}.Altitude (1000 ft).mean') as f:
                for line in f:
                    if not line.startswith('#'):
                        wmp.append(float(line.split()[1]))
            wmpsp = []
            with open(f'{root}/{case}/KL{prop.upper()}.Altitude (1000 ft).meansigma_p') as f:
                for line in f:
                    if not line.startswith('#'):
                        wmpsp.append(float(line.split()[1]))

            wm = [abs(wme[i] - wmp[i]) for i in range(len(wme))]
            wsp = [abs(wmpsp[i] - wmp[i]) for i in range(len(wmp))]
            ws = [abs(wsp[i] - wse[i]) for i in range(len(wmp))]

        # get number of samples
        with open(f'{root}/{case}/database') as f:
            lines = f.readlines()
            wc = int(lines[2].split()[0])

        datdict[ft]['s'].append(wc)
        datdict[ft]['mean'].append(np.array(wm))
        datdict[ft]['sigma'].append(np.array(ws))
        datdict[ft]['pdf'].append(np.array(wp))

# sort each dict
for key, cdict in datdict.items():
    indices = np.argsort(cdict['s'])
    # import pdb; pdb.set_trace()
    cdict['s'] = np.array(cdict['s'])[indices]
    cdict['mean'] = np.array(cdict['mean'])[indices]
    cdict['sigma'] = np.array(cdict['sigma'])[indices]

# tm = 78.13603476287041
# ts = 0.1485741352415892


# datd4m = np.array([78.14459829870214, 78.13084323997869, 78.13711406730846])
# datd4s = np.array([0.16621340157402, 0.16386430194021, 0.13821609977163])

# dats4m = np.array([78.17430194911633, 78.13508189087734, 78.06498491114561, 78.16114492405056, 78.12782063257961, 78.16800262823406])
# dats4s = np.array([0.15228149394402, 0.11920698119408, 0.23376874522297, 0.04000032038594, 0.16447480977716, 0.07161710114223])

# nd4s = [81, 625, 6561]
# ns4s = [41, 137, 401, 1105, 2929, 7537]
for prop in props:
    for key, cdict in datdict.items():
        if cdict['prop'] == prop:
            plt.plot(cdict['s'], np.sum(cdict['mean'], axis=1)/cdict['mean'].shape[1], label = f"KL {cdict['KL']}, {cdict['qtype']}")
    # [np.sum(cdict['mean'][i]) for i in range(len(cdict['s']))]
    plt.ylabel(f'KL {prop.upper()} Field Mean Error (L1)')
    plt.xlabel('Num Samples')
    plt.yscale('log')
    plt.legend()
    plt.title('KL expansion mean error sample convergence')
    plt.savefig(f'{root}/pictures/kl{prop.upper()}meanconv.png', bbox_inches='tight')
    plt.clf()


    for key, cdict in datdict.items():
        if cdict['prop'] == prop:
            plt.plot(cdict['s'], np.sum(cdict['sigma']/cdict['mean'].shape[1], axis=1), label = f"KL {cdict['KL']}, {cdict['qtype']}")
    # [np.sum(cdict['sigma'][i]) for i in range(len(cdict['s']))]
    plt.yscale('log')
    plt.ylabel(f'KL {prop.upper()} Field Sigma Error (L1)')
    plt.xlabel('Num Samples')
    plt.legend()
    plt.title('KL expansion sigma error sample convergence')
    plt.savefig(f'{root}/pictures/kl{prop.upper()}stdconv.png', bbox_inches='tight')
    plt.clf()
    # import pdb; pdb.set_trace()