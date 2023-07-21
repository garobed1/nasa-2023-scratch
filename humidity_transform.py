import numpy as np
from scipy.stats import boxcox, yeojohnson
from scipy.special import inv_boxcox#, inv_yeojohnson

from scipy.optimize import minimize_scalar

"""
Contain all possible transformations of humidity in here, to clean up quest2klfield

So far we have:
sqrt transform
dew point transform (requires a dry bulb temperature)

Hard to say which might get better results

Hum is given in percent, temps should be given in F

Property format:
Hum.shape = (number_of_path_points, number_of_paths)
aka each column corresponds to a path along altitudes
"""

selector = 0

# global parameter that holds the optimal lambda for the box cox transform
lam_opt = np.zeros(1)

def hum2dew(hum, temp, arg_units='F', return_units='F'):

    if selector == 0:
        temp_c = temp
        if arg_units == 'F':
            #Convert temp from Fahrenheit to Celsius
            temp_c = (temp - 32)/1.8 

        #Calculate dew point using equation found here: https://en.wikipedia.org/wiki/Dew_point
        dew_pt_c = (257.14*(np.log(hum/100*np.exp((18.678-temp_c/234.5)*(temp_c/(257.14+temp_c))))))/(18.678-np.log(hum/100*np.exp((18.678-temp_c/234.5)*(temp_c/(257.14+temp_c)))))

        dew_pt = dew_pt_c
        if return_units == 'F':
            #convert dew point to Fahrenheit
            dew_pt = 1.8*dew_pt_c+32
    
    if selector == 1:
        dew_pt = np.sqrt(hum)
    
    if selector == 2:
        dew_pt = np.log(hum)
    
    if selector == 3:
        N = hum.shape[0]
        lam_opt.resize(N)
        dew_pt = np.zeros_like(hum)
        for i in range(N):
            dew_pt[i], lam_opt[i] = boxcox(hum[i,:], optimizer=optimizer)
            # dew_pt[i], lam_opt[i] = yeojohnson(hum[i,:])#, optimizer=optimizer)

    return dew_pt

def dew2hum(dew, temp, arg_units='F', dew_units='F'):

    # Ensure that dew point temperature is clipped to be less than its corresponding
    # dry bulb temperature
    if selector == 0:
        if any((dew - temp)>0) :
            dew = np.clip(dew, a_min=None, a_max=temp)

        temp_c = temp
        dew_c = dew
        if arg_units == 'F':
            temp_c = f2c(temp)
        if dew_units == 'F':
            dew_c = f2c(dew)

        # work = 100*(6.112*np.exp(17.67*(dew_c)/(dew_c+243.5))/(6.112*np.exp(17.67*(temp_c)/(temp_c+243.5))))
        b = 18.678
        c = 257.14
        work = 100*np.exp(b*dew_c/(dew_c+c) - b*temp_c/(temp_c+c))
    if selector == 1:
        work = dew**2
    
    if selector == 2:
        work = np.exp(dew)
    
    if selector == 3:
        N = dew.shape[0]
        work = np.zeros_like(dew)
        for i in range(N):
            work[i] = inv_boxcox(dew[i,:], lam_opt[i])
            # work[i] = inv_yeojohnson(dew[i,:], lam_opt[i])

    # final pass
    if any((work - 100.)>0.) :
        work = np.clip(work, a_min=None, a_max=100.)

    return work


def f2c(temp_f):
    return (temp_f - 32)/1.8 

def c2f(temp_c):
    return 1.8*temp_c+32

options = {'xatol': 1e-12}  # absolute tolerance on `x`
def optimizer(fun):
    return minimize_scalar(fun, bounds=(0.2, 4.),
                                    method="bounded", options=options)
