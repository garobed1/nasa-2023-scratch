import numpy as np
import matplotlib.pyplot as plt
import os, sys

"""
Give directories with _kl?l* as arguments (? is a specific number, * is an actual wildcard)

This one handles post cart3d/sboom realization output statistics errors
"""
metric = 'seld'
qtypenamedict = {
    'dcc':'Dense CC',
    'spcc':'Sparse CC',
    'dgp':'Dense GP',
    'spgp':'Sparse GP',
    'mc':'Monte Carlo',
    'lhs':'LHS'
}
qtypelist = qtypenamedict.keys()
datdict = {}
caselist = sys.argv[1:]
root = os.getcwd()
for case in caselist:
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

tm = 78.13603476287041
ts = 0.1485741352415892


# datd4m = np.array([78.14459829870214, 78.13084323997869, 78.13711406730846])
# datd4s = np.array([0.16621340157402, 0.16386430194021, 0.13821609977163])

# dats4m = np.array([78.17430194911633, 78.13508189087734, 78.06498491114561, 78.16114492405056, 78.12782063257961, 78.16800262823406])
# dats4s = np.array([0.15228149394402, 0.11920698119408, 0.23376874522297, 0.04000032038594, 0.16447480977716, 0.07161710114223])

# nd4s = [81, 625, 6561]
# ns4s = [41, 137, 401, 1105, 2929, 7537]

for key, cdict in datdict.items():
    plt.plot(cdict['s'], abs(cdict['mean']-tm), label = f"KL {cdict['KL']}, {cdict['qtype']}")

plt.ylabel(f'{metric} Mean Error')
plt.xlabel('Num Samples')
plt.yscale('log')
plt.legend()
plt.title('KL expansion mean error sample convergence')
plt.savefig('meanconv.png', bbox_inches='tight')
plt.clf()


for key, cdict in datdict.items():
    plt.plot(cdict['s'], abs(cdict['sigma']-ts), label = f"KL {cdict['KL']}, {cdict['qtype']}")
plt.yscale('log')
plt.ylabel(f'{metric} Std Error')
plt.xlabel('Num Samples')
plt.legend()
plt.title('KL expansion sigma error sample convergence')
plt.savefig('stdconv.png', bbox_inches='tight')