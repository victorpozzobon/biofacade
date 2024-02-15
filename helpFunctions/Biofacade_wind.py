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

# Convection coefficient on the facade
def facadeConvectionCoefficient(U_10, IncidenceAngle, r_Building, r_Meteo, H_Building, H_Meteo):
    # Output
    h = 0 # W/m²/K
    
    # Velocity on the facade, m/s
    U_rescaled = U_10 * (r_Building/r_Meteo)**0.0706 * \
        np.log( (H_Building + r_Building) / r_Building) / \
            np.log( (H_Meteo + r_Meteo) / r_Meteo)

    # Incidence angles
    i_below = int(np.floor(IncidenceAngle/30.))
    i_above = i_below + 1
    theta_span = IncidenceAngle - i_below * 30.
    
    # Convection coefficient
    #possibleTheta = np.array([0,    30,     60,   90,  120,  150,  180,  210,  240,  270,  300,  330,  360])
    A_Theta       = np.array([4.90, 4.63, 4.25, 2.78, 1.44, 1.85, 2.25, 1.85, 1.44, 2.78, 4.25, 4.63, 4.90])
    B_Theta       = np.array([0.86, 0.87, 0.88, 0.87, 0.83, 0.84, 0.84, 0.84, 0.83, 0.87, 0.88, .087, 0.86])
    h_below = A_Theta[i_below] * U_rescaled ** B_Theta[i_below]
    h_above = A_Theta[i_above] * U_rescaled ** B_Theta[i_above]
    h = h_below + (h_above-h_below) * theta_span / 30.
    
    return h

# dirWind_new => 0 - comes from the North, 90° comes from the East
# velWind_new in m/s 