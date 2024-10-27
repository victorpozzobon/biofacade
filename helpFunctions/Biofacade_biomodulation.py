# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
## --- Clear all --- ###
# import sys
# sys.modules[__name__].__dict__.clear()
## --- --------- --- ###

import numpy as np

# Metabolism modulation based on temperature
def muOfT(Sorokiniana, T):
    
    mu_max = 1 # dummy value to get a fraction between [0,1]
    Tmin = 0
    Topt = 0
    Tmax = 0
    
    if Sorokiniana:
        # Sorokiniana
        # Source : Bernard et Rémond 2012
        Tmin = 5.2
        Topt = 38.7
        Tmax = 45.8
    else:
        # Param Mayo data that you fitted
        # min opt  max => -22.0, 31.3, 45.6
        Tmin = -22.0
        Topt = 31.3
        Tmax = 45.6
   
    f_T = (Topt - Tmin) * (T - Topt)
    g_T = (Topt - Tmax) * (Topt + Tmin - 2 * T)
    mu_T = mu_max * (T - Tmax) * (T - Tmin)**2 / ((Topt - Tmin) * (f_T - g_T)) # Fraction of mu_max, with NEGATIVE VALUES
    mu_T_culled = mu_T * (T > Tmin) * (T < Tmax) # Fraction of mu_max, culled at 0
    
    return mu_T_culled

# Pigment content based on averaged illumination
def PigOfI(A, b, a, I):
    return A * np.exp( - I / (a * I + b))

def PigEquilibrium(pigmentLetter, I):
    
    EqContent = 0
    # Chlorophyll a
    if pigmentLetter == "a":
        EqContent = PigOfI(25.2, 56.7, 0.361, I)
        
    # Chlorophyll b
    elif pigmentLetter == "b":
        EqContent = PigOfI(12.0, 59.4, 0.364, I)
        
    # Lutein    
    elif pigmentLetter == "l":
        EqContent = PigOfI(6.53, 17.6, 0.315, I)
        
    else:
        print('Problem on pigment selection')
    
    return EqContent