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

def SnellRefraction(theta_i, n_i, n_r):
    # compute sin theta_r
    sin_theta_r = n_i / n_r * np.sin(theta_i)
    theta_r = -1
    
    # verifying and computing refraction
    if sin_theta_r < 0:
        theta_r = -2
        print("Huge problem in Snell's law")
    elif sin_theta_r <= 1:
        theta_r = np.arcsin(sin_theta_r)
        
    return theta_r # -2 => Huge problem, -1 => total reflection, [0-pi/2] => refraction

def reflectivity_one_interface(theta_i, theta_r):
    reflectivity = 1.0
    # print("Angles")
    # print(theta_i)
    # print(theta_r)
    # compute reflectivity only if theta_r is valid
    if (theta_r >= 0) and (theta_r <= np.pi/2) and (theta_i >= 0) and (theta_i <= np.pi/2) and (theta_i + theta_r > 0):
        reflectivity = 0.5 *( np.tan(theta_i - theta_r) **2 / np.tan(theta_i + theta_r) **2 +
                              np.sin(theta_i - theta_r) **2 / np.sin(theta_i + theta_r) **2 )
    elif (theta_i + theta_r == 0):
        reflectivity = 0.
    # print("Reflectivity")
    # print(reflectivity)
    
    return reflectivity

def transmission_multiple_interface(theta_i, n_s):
    # n_s is an numpy array of the successive refraction indexes
    n_interface = len(n_s) - 1
    transmission = 1.0
    theta_i_current = theta_i
    for i in range(0, n_interface):
        theta_r = SnellRefraction(theta_i_current, n_s[i], n_s[i+1])
        transmission = transmission * (1.0 - reflectivity_one_interface(theta_i_current, theta_r))
        if transmission < 1e-15:
            transmission = 0
            break
        theta_i_current = theta_r
        # print(transmission)
    return transmission

def transmission_multiple_interface_cos_distribution(n_s):
    # compute the transmitivity for cos distibuted ray
    dtheta = 0.01
    transmission_sum = 0
    for theta_i in np.arange(0, 90, dtheta):
        # print(theta_i)
        # print(transmission_sum)
        angle_rad = (theta_i + dtheta / 2 ) / 180 * np.pi
        transmission_sum += np.cos(angle_rad) * transmission_multiple_interface(angle_rad, n_s)  * (dtheta / 180. * np.pi)

    return transmission_sum

# print(transmission_multiple_interface(45 / 180. * np.pi, [1.00028276, 1.5082, 1.00028276, 1.5082, 1.3390]))
# print(transmission_multiple_interface_cos_distribution([1.0, 1.3, 1.0, 1.5, 1.3]))
# print(transmission_multiple_interface_cos_distribution([1.3, 1.5, 1.0, 1.3, 1.0]))
# print(transmission_multiple_interface_cos_distribution([1.00028276, 1.5082, 1.00028276, 1.5082, 1.339]))
# print(transmission_multiple_interface_cos_distribution([1.339, 1.5082, 1.00028276]))
# print(transmission_multiple_interface_cos_distribution([1.339, 1.5082, 1.00028276, 1.5082, 1.00028276]))
# print(transmission_multiple_interface_cos_distribution([1.00028276, 1.5082, 1.339]))
# print(transmission_multiple_interface(0, [1.00028276, 1.5082, 1.00028276, 1.5082, 1.339]))
