import numpy as np
import matplotlib.pyplot as plt
import os, sys

"""
After running quest_post.csh and outputting mean/mean+1sigma into stats/ for each project

Give directories with _kl?l* as arguments (? is a specific number, * is an actual wildcard)

This one handles post cart3d/sboom realization output statistics errors

Ex.
cd ~/temp_cc
python3 ~/nasa-2023-scratch/output_conv_plot.py sboom_atm_uq_temp_dcc_kl4l* sboom_atm_uq_temp_lhs_kl4l*
"""
# output measure to check for
metric = 'seld'

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
for case in caselist:

    # retrieve properties from directory names
    name = case.split('_')
    kl = name[-1]
    qtype = name[-2]
    klp = kl[2] # KL param number
    kll = kl[4] # grid level

    if qtype not in qtypelist:
        qtype = 'dcc'

    ft = klp + qtype

    if ft not in datdict:
        datdict[ft] = {}
        datdict[ft]['qtype'] = qtypenamedict[qtype]
        datdict[ft]['KL'] = klp
        datdict[ft]['s'] = []
        datdict[ft]['mean'] = []
        datdict[ft]['sigma'] = []

    # get mean
    # get stdv
    with open(f'{root}/{case}/stats/{metric}.Off-Track Angle (deg.).mean') as f:
        lines = f.readlines()
        wm = float(lines[0].split()[-1])
    with open(f'{root}/{case}/stats/{metric}.Off-Track Angle (deg.).meansigma_m') as f:
        lines = f.readlines()
        wmms = float(lines[0].split()[-1])
    ws = wm - wmms
    
    # get number of samples
    with open(f'{root}/{case}/database') as f:
        lines = f.readlines()
        wc = int(lines[2].split()[0])

    datdict[ft]['s'].append(wc)
    datdict[ft]['mean'].append(wm)
    datdict[ft]['sigma'].append(ws)

# sort each dict
for key, cdict in datdict.items():
    indices = np.argsort(cdict['s'])
    # import pdb; pdb.set_trace()
    cdict['s'] = np.array(cdict['s'])[indices]
    cdict['mean'] = np.array(cdict['mean'])[indices]
    cdict['sigma'] = np.array(cdict['sigma'])[indices]

"""
These are "exact" values for comparison. Taken from the full data
TODO: Replace as relevant
"""
tm = 78.13603476287041
ts = 0.1485741352415892


for key, cdict in datdict.items():
    plt.plot(cdict['s'], abs(cdict['mean']-tm), label = f"KL {cdict['KL']}, {cdict['qtype']}")

plt.ylabel(f'{metric} Mean Error')
plt.xlabel('Num Samples')
plt.yscale('log')
plt.legend()
plt.title('KL expansion mean error sample convergence')
plt.savefig(f'{root}/pictures/{metric}meanconv.png', bbox_inches='tight')
plt.clf()


for key, cdict in datdict.items():
    plt.plot(cdict['s'], abs(cdict['sigma']-ts), label = f"KL {cdict['KL']}, {cdict['qtype']}")
plt.yscale('log')
plt.ylabel(f'{metric} Sigma Error')
plt.xlabel('Num Samples')
plt.legend()
plt.title('KL expansion sigma error sample convergence')
plt.savefig(f'{root}/pictures/{metric}stdconv.png', bbox_inches='tight')