import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, sys

plt.rcParams['font.size'] = 16


output = 'PL'
thresh = 1275.0

outputlabels = {
    'PL': 'PL (dB)',
    'PL2': 'PL (dB)',
}

# look for its files
root = os.getcwd()
dirlist = os.listdir()
outlist = []
for file in os.listdir():
    if file.startswith(output):
        outlist.append(file)

for file in outlist:
    # PDF
    if file.endswith('.cdfpdf'):
        x = []
        y = []
        c = []
        with open(file, 'r') as pdffile:
            for line in pdffile:
                if line.startswith('0'):
                    continue
                linesp = line.split()
                x.append(float(linesp[0]))
                c.append(float(linesp[1]))
                y.append(float(linesp[2]))
    # MEAN
    if file.endswith('.mean'):
        mean = 0
        with open(file, 'r') as mfile:
            line = mfile.readline()
            mean = float(line.split()[1])

    # MP1S
    if file.endswith('.meansigma_p'):
        mp1s = 0
        with open(file, 'r') as mfile:
            line = mfile.readline()
            mp1s = float(line.split()[1])

    # MM1S
    if file.endswith('.meansigma_m'):
        mm1s = 0
        with open(file, 'r') as mfile:
            line = mfile.readline()
            mm1s = float(line.split()[1])

sigma = mean - mm1s


# normalize
x = np.array(x)
y = np.array(y)
y = y / np.max(y)
nsamp = len(x)
colourmap = mpl.colormaps['rainbow']
# plt.grid()
plt.plot(x, y)
# fill with color gradient
for i in range(nsamp - 1):
    plt.fill_between([x[i], x[i+1]],
                     [y[i], y[i+1]],
                     color=colourmap(y[i]),
                    #  ,alpha=0.6,
                     edgecolor=(0,0,0,0.0))

# C612 mach only stats
momean = 72.00291460069828
momp1s = 72.72188287011573
momm1s = 71.28394633128082

# mean, std, other
plt.gca().axvline(mean, linestyle= '-', color = 'k')
plt.gca().axvline(mp1s, linestyle= '--', color = 'k')
plt.gca().axvline(mm1s, linestyle= '--', color = 'k')
# plt.gca().axvline(momean, linestyle= '-', color = 'g')
# plt.gca().axvline(momp1s, linestyle= '--', color = 'g')
# plt.gca().axvline(momm1s, linestyle= '--', color = 'g')
# plt.gca().axvline(thresh, linestyle= '--', color = 'orange')
plt.plot([],[], '-k', label = rf'$\mu$')
plt.plot([],[], '--k', label = rf'$\mu \pm 1\sigma$')
# plt.plot([],[], '-g', label = rf'$\mu$ (Mach)')
# plt.plot([],[], '--g', label = rf'$\mu \pm 1\sigma$ (Mach)')
# plt.plot([],[], '--r', label = f'{thresh} dB Target')
plt.xlabel(outputlabels[output])
plt.ylabel('Normalized PDF Value')
plt.legend(loc=2,fontsize=10)
plt.savefig(f'{output}_PDF_plot.png', bbox_inches='tight')

cdf = 1
tcomp = x - thresh
for i in range(len(tcomp)):
    if tcomp[i] > 0:
        cdf = c[i]
        break
print(mean)
print(sigma)
print(f'{100*(1.-cdf)} % chance to exceed {thresh} dB target')