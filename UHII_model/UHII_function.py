#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:29:50 2023

@author: victor
"""

def UHII_t(UHII_Daily_night_to_come, UHII_Daily_before, t_current, t_sundawn, basal_level = 0.5):
   
    # UHII model based on: 
    # Montávez, J. P., Rodríguez, A., & Jiménez, J. I. (2000). A study of the Urban Heat Island of Granada. International Journal of Climatology, 20(8), 899‑911. https://doi.org/10.1002/1097-0088(20000630)20:8<899::AID-JOC433>3.0.CO;2-I
    # Lai, J., Zhan, W., Huang, F., Voogt, J., Bechtel, B., Allen, M., Peng, S., Hong, F., Liu, Y., & Du, P. (2018). Identification of typical diurnal patterns for clear-sky climatology of surface urban heat islands. Remote Sensing of Environment, 217, 203‑220. https://doi.org/10.1016/j.rse.2018.08.021
    
    import numpy as np
    
    # Idea: a straight line between two sinusoids
    hour_since_dawn = t_current - t_sundawn # h
    UHII = -100 # K unlikeky to happen, to detect problem
    # Correction to have a smooth transition between the two modes
    if UHII_Daily_before < basal_level : 
        UHII_Daily_before = basal_level
        
        
    # Determine UHII value as a funtion of time
    if UHII_Daily_night_to_come >= basal_level:
        if hour_since_dawn >= 0 and hour_since_dawn <= 10:
            UHII = basal_level + (UHII_Daily_before - basal_level) * (np.cos(np.pi * (hour_since_dawn - 0) / 10) +1)/2
        elif hour_since_dawn > 10 and hour_since_dawn <= 15:
            UHII = basal_level + (UHII_Daily_night_to_come - basal_level) * (-np.cos(np.pi * (hour_since_dawn - 10) / 5) +1)/2
        elif hour_since_dawn > 15:
            UHII = UHII_Daily_night_to_come 
        else:
            print("Daily hour problem ...")
    
    # Determine UCII value as a funtion of time
    if UHII_Daily_night_to_come < basal_level: 
        if hour_since_dawn >= 0 and hour_since_dawn <= 10:
            UHII = UHII_Daily_night_to_come + (UHII_Daily_before - UHII_Daily_night_to_come) * (np.cos(np.pi * (hour_since_dawn - 0) / 10) +1)/2
        elif hour_since_dawn > 10 and hour_since_dawn <= 15:
            UHII = UHII_Daily_night_to_come + (basal_level - UHII_Daily_night_to_come) * (-np.cos(np.pi * (hour_since_dawn - 10) / 5) +1)/2
        elif hour_since_dawn > 15:
            UHII = basal_level 
        else:
            print("Daily hour problem ...")
    
    return UHII
