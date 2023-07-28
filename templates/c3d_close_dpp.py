#!/usr/bin/env python3
#
#  $Id: c3d_close_dpp.py,v 1.8 2023/06/14 19:46:55 mnemec Exp $
#
'''  ...Script for closing near field presure signals, blends out signal
        and then closes it using either a linear ramp or a quadratic.
         (1) optionally plots orig and closed near field signals
         (2) optionally plots in batch mode
         (3) optionally runs sBOOM
         (4) if sBOOM is run, then also plots ground signal (interactive or batch)
         (5) linearizes the close-out function for sBOOM gradients 

        NOTES:
        o need to set the name of sBOOM executable down in main script,  search "SBOOM"
        o if you have problems with matplotlib, comment out the import and the plotting
          at the bottom of the main script, search "PLOTTING"
'''
# =======================================================
import argparse
import os
import math
import matplotlib.pyplot as plt   #       ...optional -- turn off plotting below
import numpy as np                #       ...optional -- just need for plotting
from   os.path import exists


def closeoutFunction(ys, x, y, cs, ce):
    """Signal closeout function using base curve plus damping function
       Returns the new value of y (dpp) at x
            ys := initial y value at start of closeout y(index(ce))
            x  := current x location
            cs := closeoutStart location (in x)
            ce := closeoutEnd location (in x)
            NOTE: for the damping function used below,
                  (ce-cs) should be ~20 meters for the X-59
    Hint:
      istart = next(i for i in range(len(sensorPts)) if sensorPts[i][0]>cs)
      ys     = sensorPts[istart][1]
    """
    y_new = y
    if x > cs and x < ce:
        t = 1 - (x - cs) / (ce - cs)  # t goes from 1 to 0
        f = math.exp(-400*((t-1)**2)) # fast exponential damping
        y_new = t*ys + f*(y-ys)       # linear    base decay curve
        #y_new = t*t*ys + f*(y-ys)    # quadratic base decay curve
    elif x>= ce:
        y_new = 0.
    return y_new

def lin_closeoutFunction_ys(dy, x, cs, ce):
    dy_new = dy
    if x > cs and x < ce:
        t = 1 - (x - cs) / (ce - cs)
        f = math.exp(-400*((t-1)**2))
        dy_new = (t-f)*dy
    elif x >= ce:
        dy_new = 0.
    return dy_new

def lin_closeoutFunction_y(dy, x, cs, ce):
    dy_new = dy
    if x > cs and x < ce:
        t = 1 - (x - cs) / (ce - cs)
        f = math.exp(-400*((t-1)**2))
        dy_new = f*dy
    elif x >= ce:
        dy_new = 0.
    return dy_new

def findSigStart(sensor, tol = 5e-4):
    """Find where the signal first 'starts'/significantly deviates from 0
        Returns the x value where the signal starts
            sensor := list of tuples containing sensor location and value
        NOTE: Iterate from start and check deviation after a tolerance
    """
    N = len(sensor)
    i = 0
    while i < N - 1:
        cury = sensor[i][1]
        nexty = sensor[i+1][1]
        move = nexty-cury
        if move > tol:
            return sensor[i][0]
        i += 1

    print('WARNING: Could not find start of signal, returning first point instead')
    return sensor[0][0]

def findSigStartMach(infile, sdist):
    """Find relative signal position based on mach angle
        Returns reference x value for signal location
            infile := c3d cntl input file string
            sdist  := orthogonal distance from mach cone origin to sensor
    """
    with open(infile) as file:
        for line in file:
            if line.startswith("Mach"):
                words = line.split()
                mach = float(words[1])

    # mach angle formula
    mth = np.arcsin(1./mach)

    return sdist/np.tan(mth)

def findGroundCross(ground, nsafe = 10):
    """Find zero crossing of ground signal
        Returns the x value where the zero crossing occurs
            ground       := ordered list of tuples for ground signal
            nsafe        := make sure it stays negative after this many points
        NOTE: Iterate from start and check if signal and its slope go negative
    """

    N = len(groundPoints)
    i = 0
    ntol = 1e-10
    while i < N - 1:
        cury = ground[i][1]
        nexty = ground[i+1][1]
        move = nexty-cury
        if move < -ntol and nexty < 0. and cury > 0.:
            # before returning, check if it stays negative for nsafe points
            if ground[i+nsafe][1] < -ntol:
                return ground[i][0]
        i += 1

    print('WARNING: Could not find zero crossing, returning first point instead')
    return ground[0][0]
               


# ==============================
#       main script
# ==============================
'''  ...Script for closing near field presure signals, blends out signal
        and then closes it using either a linear ramp or a quadratic.
         (1) optionally plots orig and closed near field signals
         (2) optionally plots in batch mode
         (3) optionally runs sBOOM
         (4) if sBOOM is run, then also plots ground signal (interactive or batch)
         (5) linearizes the close-out function for sBOOM gradients

        NOTES:
        o Set the name of sBOOM executable down in main script,  search "SBOOM"
        o if you have problems with matplotlib, comment out the import and the plotting
          at the bottom of the main script, search "PLOTTING"
'''
#                                    ...1. Parse Args
parser = argparse.ArgumentParser(description="Preprocess sensor, optionally run sBoom or plot")
parser.add_argument("sensorName", help="name of Cart3D lineSensor file")
parser.add_argument("-cs", "--closeoutStart", type=float,
                    help="distance along sensor where closeout of signal will start")
parser.add_argument("-ce", "--closeoutEnd",   type=float,
                    help="location where closeout of signal will end <default:=closeoutStart+15>")
parser.add_argument("-ap", "--aftPadDist",    type=float, default=5.,
                    help="pad distance after sensor & closeout <default:=5>")
parser.add_argument("-nc", "--closeoutNpts",  type=int, default=10,
                    help="number of points in aft pad")
parser.add_argument("-df", "--frontPadDist",  type=float,
                    help="pad distance ahead of sensor <default:= noFrontPad>")
parser.add_argument("-nf", "--frontPadNpts",  type=int, default=10,
                    help="number of points in pad ahead of sensor <default:=10> if padding front")
parser.add_argument("-sv", "--setValue",      type=float,
                    help="option to set value at start of closeout (Use Caution!)")
parser.add_argument("-st", "--shiftAlignTol", type=float, default = 5e-4,
                    help="set tolerance of signal start detection")
parser.add_argument("-sd", "--sensorDist", type=float, 
                    help="distance of sensor from mach cone origin")
parser.add_argument("-sb", "--run_sBOOM",
                    help="run sBOOM with modSig.dat (requires existing presb.input)",    action="store_true")
parser.add_argument("-sf",  "--shiftAlignFront",
                    help="align all signals by where they start referring to nominal",   action="store_true")
parser.add_argument("-sg",  "--shiftAlignGround",
                    help="align all ground signals at the 0 crossing (for varying mach_",action="store_true")
parser.add_argument("-p",  "--plotLive",
                    help="plot (interactive) original & modified near field signals",    action="store_true")
parser.add_argument("-pb", "--plotBatch",
                    help="plot to *pdf (batch) original & ,modified near field signals", action="store_true")
parser.add_argument("-fd", "--finiteDifferenceCheck",
                    help="check closeout linearization with finite differences", action="store_true")
args = parser.parse_args()

# if not args.frontPadDistance == None and args.front

sensorFileName = args.sensorName
print(' o Working on ' + sensorFileName)

if None == args.closeoutStart:
    print('ERROR: Need to provide closeoutStart (--cs) on cmd line\n       Exiting(1)')
    exit(1)

#                             ....2. Read in sensor from Cart3D lineSensor file"
sensorPoints = []
with open(sensorFileName) as file:
    for line in file:
        if not line.startswith("#"):
            words = line.split()
            sensorPoints.append((float(words[3]), float(words[4])))

numSensorPts = len(sensorPoints)
origSensor   = sensorPoints.copy()

#                             ....2.5 Translate signals such that their starts line up
# By default, use case.00.00000 as the reference
if args.shiftAlignFront:
    print('Aligning near-field signals based on nominal mach angle')
    # get nominal case
    sensorPoints0 = []
    # sensorFileName0 = '../../case.00.00000/' + sensorFileName
    sensorFileName0 = '../../case.00.00000/CFD/BEST/FLOW/' + sensorFileName.split('/')[-1]
    with open(sensorFileName0) as file:
        for line in file:
            if not line.startswith("#"):
                words = line.split()
                sensorPoints0.append((float(words[3]), float(words[4])))
    origSensor0 = sensorPoints0.copy()

    # find starts of nominal case and current case
    tol = args.shiftAlignTol
    sigstart0 = findSigStart(origSensor0, tol=tol)
    sigstart = findSigStart(origSensor, tol=tol)

    # alternatively, try computing based on mach angle
    inputFileName0 = '../../case.00.00000/CFD/input.cntl'
    inputFileName = '../CFD/input.cntl'

    sdist = args.sensorDist
    sigstart0 = findSigStartMach(inputFileName0, sdist)
    sigstart = findSigStartMach(inputFileName, sdist)
    

    sigshift = sigstart0 - sigstart

    # shift current case signal and continue
    origSensorShift = []
    for i in range(len(sensorPoints)): origSensorShift.append((origSensor[i][0] + sigshift, origSensor[i][1])) 
    origSensor = origSensorShift
    sensorPoints = origSensorShift.copy()

#           ...know we have a filename and closout location, start working
#                    finish processing inputs, and setup the distances for
aftFlat = args.aftPadDist  #                           tweaking the signal
closeoutSart = args.closeoutStart
#if None == args.closeoutEnd:
#    closeoutEnd = closeoutStart + 15  # default is 15 units after closeout start

closeoutEnd = closeoutSart + 15  if None == args.closeoutEnd else args.closeoutEnd

# if we shifted the signal, adjust the closeout distances as well
if args.shiftAlignFront:
    closeoutSart += sigshift
    closeoutEnd += sigshift
    aftFlat -= sigshift

backPadDist = max(0., closeoutEnd + aftFlat - sensorPoints[-1][0])
backPadNpts = args.closeoutNpts

frontPadDist = args.frontPadDist
frontPadNpts = args.frontPadNpts

if args.shiftAlignFront:
    frontPadDist += sigshift

#                             ....3. Pad front and back if desired
if args.frontPadDist: #               <-- pad front
    xstart = sensorPoints[0][0]
    for i in range(frontPadNpts):
        x = float(i+1)/float(frontPadNpts)
        sensorPoints.insert(0, (xstart - x * x * float(frontPadDist), 0.0))

#                   ...Pad aft of sensor if required for closeout and aft flat
if backPadNpts > 0:
    xend = sensorPoints[-1][0]
    for i in range(backPadNpts):
        x = float(i+1)/float(backPadNpts)
        sensorPoints.append((xend + x * float(backPadDist), 0.0))

#                           ....4. Apply closeout, fading signal to base closeout
cs = closeoutSart #                               curve using gaussian multiplier
ce = closeoutEnd

#             ...set closeout start location and (optionally) start value
istart = next(i for i in range(len(sensorPoints)) if sensorPoints[i][0]>cs)
ys     = sensorPoints[istart][1] if None == args.setValue else args.setValue
#print("==> set y-value at closeout to " + str(ys))

for i in range(istart, len(sensorPoints)):
    pt = sensorPoints[i]
    newPt = (pt[0], closeoutFunction(ys, pt[0], pt[1], cs, ce))
    sensorPoints[i] = newPt

#                          ...write out nearfield file (modSig.dta) for sBoom
if exists("modSig.dat"):
        print('   found existing modSig.dat, moving to "old_modSig.dat"')
        os.system('\mv -f modSig.dat old_modSig.dat')
with open("modSig.dat", "w") as file:
    file.write("# Number of Waveforms = 1\n")
    file.write("# Signature from " + sensorFileName + " ( " + os.getcwd() + " )\n")
    file.write("# Number of Points = " + str(len(sensorPoints)) + "\n")
    for pt in sensorPoints:
        if pt[0] < 0:
            file.write("{x:.9E}".format(x=pt[0]))
        else:
            file.write(" {x:.9E}".format(x=pt[0]))
        if pt[1] < 0:
            file.write("  {p:.9E}\n".format(p=pt[1]) )
        else:
            file.write("   {p:.9E}\n".format(p=pt[1]) )
print("   wrote modSig.dat with modified nearfield signature ("+ str(len(sensorPoints))+ " points)",)

#                                            ... -- optionally run_sBOOM
# cmd = '/usr/bin/time -f "mem=%M elapsed=%E" sboom_2.86 -v -d -e -U -F ' + str(stepFactor)
if args.run_sBOOM:
    if exists("SBground.sig"):
            print('==> sBOOM Found existing SBground.sig --> moving to "old_SBground.sig"')
            os.system('\mv -f SBground.sig old_SBground.sig')
    if exists("loud.dat"):
            os.system('\mv -f loud.dat old_loud.dat')
    #                                    ...optionally run sBOOM
    cmd = 'sboom -v | tee sboom.OUT' #   <-- SET SBOOM command here (v2.9.0)a
    print("running sBOOM with:\n   % " + cmd)
    os.system(cmd)
    if exists("loud.dat"):
        os.system('cat loud.dat')
    else:
        print('sBOOM didnt produce a loud.dat file -- verify that it ran')

# if mach number varies, need to align signals at ground
if args.shiftAlignGround:
    print('Aligning signal with zero crossing of nominal case')
    # get nominal case
    groundFileName0 = '../../case.00.00000/sboom/SBground.sig'    
    groundFileName = 'SBground.sig'
    if exists(groundFileName0) and exists(groundFileName):
        groundPoints0 = []
        with open(groundFileName0) as file:
            for line in file:
                if not line.startswith("#"):
                    words = line.split()
                    groundPoints0.append((float(words[0]), float(words[1])))
        # get current case
        groundPoints = [] 
        with open(groundFileName) as file:
            for line in file:
                if not line.startswith("#"):
                    words = line.split()
                    groundPoints.append((float(words[0]), float(words[1])))
        # find nominal and current zero crossings
        groundCross0 = findGroundCross(groundPoints0)
        groundCross = findGroundCross(groundPoints)

        groundShift = groundCross0 - groundCross
        
        with open('SBground.sig', 'w') as file:
            for dat in groundPoints:
                file.write(f'{(dat[0]+groundShift):.16e} {dat[1]:.16e}\n')

    else:
        print('Case ground solutions not found, no alignment performed')


# post-process gradient array

from pathlib import Path
import re

gradient = []
origG    = [] 

gfile = Path('Gradient.plt')
top   = re.compile('^[#a-zA-Z]')

if gfile.is_file():
    print(' o Detected ' + str(gfile))
    print('   linearizing close-out function, preparing dJdP.dat')
    with gfile.open('r') as file:
        for line in file:
            if not top.match(line):
                words = line.split()
                gradient.append(float(words[0]))

    if len(gradient) != len(sensorPoints):
        print('ERROR: length of gradient ({0:d}) inconsistent with line sensor ({1:d})'
              .format(len(gradient), len(sensorPoints))) 

    origG = gradient.copy()

    for i in range(istart, len(sensorPoints)):
        pt = sensorPoints[i]
        gradient[i] = lin_closeoutFunction_y(origG[i], pt[0], cs, ce)
    
    # full linearization, uncomment this section
#    if args.setValue:
#        print('   constant start-value closeout')
#        for i in range(istart, len(sensorPoints)):
#            pt = sensorPoints[i]
#            # gradient[i] = 0
#            gradient[i] = lin_closeoutFunction_y(origG[i], pt[0], cs, ce)
#    else:
#        print('   linearizing closeoutFunction')
#        gradient[istart] = 0
#        for i in range(istart, len(sensorPoints)):
#            pt = sensorPoints[i]
#            gradient[istart] += lin_closeoutFunction_ys(origG[i], pt[0], cs, ce)
#
#        for i in range(istart, len(sensorPoints)):
#            pt = sensorPoints[i]
#            if i == istart:
#                gradient[i] += lin_closeoutFunction_y(origG[i], pt[0], cs, ce)
#            else:
#                gradient[i] = lin_closeoutFunction_y(origG[i], pt[0], cs, ce)

    # finite-difference debug, save before removing padding
    fdg = gradient.copy()
                
    # undo padding        
    if args.frontPadDist:
        print('   removing front pad ', frontPadNpts)
        gradient[0:frontPadNpts]=[]
        
    if backPadNpts > 0:
        print('   removing back pad ', backPadNpts)
        gradient[-backPadNpts:]=[]

    if len(gradient) != len(origSensor):
        print('ERROR: length of gradient ({0:d}) inconsistent with line sensor ({1:d})'
              .format(len(gradient), len(origSensor))) 
        
    with open('dJdP.dat', 'w') as file:
        for dp in gradient:
            if dp < 0:
                file.write("{:.16e}\n".format(dp))
            else:
                file.write(" {:.16e}\n".format(dp))

    print("   wrote dJdP.dat")

    if args.finiteDifferenceCheck:
        # check closeout linearization with finite differences
        fdeta = 1.e-5
        for i in range(istart, len(sensorPoints)):
            pt          = sensorPoints[i]
            hold        = pt[1]
            sensorPoints[i] = pt[0], hold + fdeta
            stp         = sensorPoints[i][1] - hold

            yP = [0]*len(sensorPoints)
            yM = [0]*len(sensorPoints)
            
            ys = sensorPoints[istart][1] if None == args.setValue else args.setValue
            
            for j in range(istart, len(sensorPoints)):
                yP[j] = closeoutFunction(ys, sensorPoints[j][0], sensorPoints[j][1],
                                         cs, ce)
            
            sensorPoints[i] = pt[0], hold - stp
            
            ys = sensorPoints[istart][1] if None == args.setValue else args.setValue

            for j in range(istart, len(sensorPoints)):
                yM[j] = closeoutFunction(ys, sensorPoints[j][0], sensorPoints[j][1],
                                         cs, ce)

            fdg[i] = 0
            for j in range(istart, len(sensorPoints)):
                fdg[i] += (yP[j] - yM[j])/(2*stp)*origG[j]

            # restore
            sensorPoints[i] = pt

#                            === PLOTTING ===
#                        exit here if no matplotlib
#               (commont out this whole block if you dont have matplotlib)
#                            ================

# colors
lightYinMin = '#6C84B1'  # is yinmin blue blended with white
deepskyblue = '#0d75f8'  # xkcd colors
darkpink    = '#cb416b'

if args.plotLive or args.plotBatch:
    print(' o Plotting results')
    
    x_new  = np.zeros(len(sensorPoints))  # ...create the numpy vectors
    y_new  = np.zeros(len(sensorPoints))
    x_orig = np.zeros(len(origSensor))
    y_orig = np.zeros(len(origSensor))

    for i in range(len(sensorPoints)):
      x_new[i]  = sensorPoints[i][0]
      y_new[i]  = sensorPoints[i][1]

    for i in range(len(origSensor)):
      x_orig[i] = origSensor[i][0]
      y_orig[i] = origSensor[i][1]

    #                             ....actually make the plots
    fig0, ax = plt.subplots(figsize=(10,5.5))
    ax.set_xlim(x_new[0],(x_new[-1]-aftFlat/2.))
    ax.set_xlabel('Distance Along Sensor', fontsize = 14)
    ax.set_ylabel('$\Delta p/p_\infty$',   fontsize = 14) # fancy symbols
    ax.yaxis.get_ticklocs(minor=True)
    ax.minorticks_on()
    ax.grid()
    ax.plot(x_orig,y_orig,linewidth=0.9, color = darkpink, label = 'rawSig.dat')
    ax.plot(x_new,y_new,linewidth=1.4, color = deepskyblue, label = 'modSig.dat')
    if args.shiftAlignFront:
        ax.axvline(sigstart0, linewidth=0.9, color = darkpink, linestyle = '--', label = 'rawSig start')
        ax.axvline(sigstart, linewidth=1.4, color = deepskyblue, linestyle = '--', label = 'modSig start')

    plt.legend(fontsize = 14)
    plt.savefig('dpp_compare.png', dpi=300, bbox_inches='tight')
    print('   plotted raw and modified nearfield signals in dpp_compare.png')
    if args.plotLive :
        plt.show(block = False)  # ...nonblocking plot

if (args.run_sBOOM and (args.plotLive or args.plotBatch) & exists("SBground.sig")):
    fig1, ax = plt.subplots(figsize=(8,5.5))
    groundSig = np.loadtxt("SBground.sig", skiprows = 3)
    xmin = groundSig[0,0]
    xmax = groundSig[-1,0]
    xstart = 0
#    xstart = groundSig[np.argmax(groundSig[:,1] > 0.001),0] # skip leader (look for press > 0.001)
#    xstart += -0.1*(xstart - xmin)  # back by 1/10 of the dist to the start

    ax.set_xlim(xstart,xmax - 0.1*(xmax - xstart))
    ax.set_xlabel('Time, [ms]', fontsize = 14)
    ax.set_ylabel('Pressure, [psf]',   fontsize = 14)
    ax.plot(groundSig[:,0],groundSig[:,1],linewidth=1.4, color = deepskyblue, label = 'SBground.dat')
    ax.yaxis.get_ticklocs(minor=True)
    ax.minorticks_on()
    ax.grid()
    plt.legend(fontsize = 14)
    plt.savefig('groundSig.png', dpi=300, bbox_inches='tight')
    print('   plotted ground signature in groundSig.png')
    if args.plotLive:
        if gfile.is_file():
            plt.show(block = False)
        else:
            plt.show()

if (gfile.is_file() and (args.plotLive or args.plotBatch)):
    x_orig = np.zeros(len(sensorPoints)) # includes padding and closeout
    y_p    = np.zeros(len(sensorPoints))
    y_orig = np.zeros(len(origG))        # linearization computed by sBOOM
    x_nf   = np.zeros(len(origSensor))
    y_pnf  = np.zeros(len(origSensor))
    y_nf   = np.zeros(len(gradient))
    
    for i in range(len(sensorPoints)):
      x_orig[i] = sensorPoints[i][0]
      y_p[i]    = sensorPoints[i][1]*4000+100
      y_orig[i] = origG[i]

    for i in range(len(origSensor)):
      x_nf[i]  = origSensor[i][0]
      y_pnf[i] = origSensor[i][1]*4000+100
      y_nf[i]  = gradient[i]
    
    fig2, ax = plt.subplots(figsize=(10,5.5))
    ax.set_xlim(x_orig[0], x_orig[-1])
    ax.set_ylim(-100, 150)  # hand tune
    ax.set_xlabel('Distance Along Sensor', fontsize = 14)
    ax.set_ylabel('$\partial J/p\'$', fontsize = 14) # fancy symbols
    ax.yaxis.get_ticklocs(minor=True)
    ax.minorticks_on()
    ax.grid()
    ax.plot(x_orig,y_p, linewidth=1.0, linestyle=':', color = 'black', label = 'sBOOM Signature')
    ax.plot(x_nf,y_pnf, linewidth=1.0, linestyle='--', color = 'black', label = 'Cart3D Signature')
    ax.plot(x_orig,y_orig,linewidth=1.0, color = deepskyblue, label = 'sBOOM Gradient')
    ax.plot(x_nf,y_nf,linewidth=1.0, color = darkpink, label = 'Cart3D Gradient')
    
    plt.legend(fontsize = 14)
    plt.savefig('gradplot.png', dpi=300, bbox_inches='tight')
    print('   plotted gradients in gradplot.png')
    if  args.plotLive:
        if args.finiteDifferenceCheck:
            plt.show(block = False)
        else:
            plt.show()
        
    if args.finiteDifferenceCheck:
        y_fd = np.zeros(len(fdg))     
        for i in range(len(fdg)):
            y_fd[i] = fdg[i]
      
        fig3, ax = plt.subplots(figsize=(10,5.5))
        ax.set_xlim(x_orig[0], x_orig[-1])
        ax.set_xlabel('Distance Along Sensor', fontsize = 14)
        ax.set_ylabel('$\partial J/p\'$', fontsize = 14) # fancy symbols
        ax.yaxis.get_ticklocs(minor=True)
        ax.minorticks_on()
        ax.grid()
        ax.plot(x_nf,y_nf, linewidth=1.0, linestyle='-', marker="o", color = 'red', label = 'G')
        ax.plot(x_orig,y_fd, linewidth=1.0, linestyle=':', marker="v", color = 'black', label = 'FD')

        plt.legend(fontsize = 14)
        plt.savefig('gradfd.png', dpi=300, bbox_inches='tight')
        print('   plotted finite difference checks in gradfd.png')
        if  args.plotLive:
            plt.show()
        
exit(0)
