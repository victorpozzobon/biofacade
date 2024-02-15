# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 08:34:15 2018

@author: victor
"""
############## ATTENTION ##############
#                                     #
#      Checked on 2023 Febr. 01       #
#               WORKING               #
#                                     #
############## ATTENTION ##############

# Based: Recommended Practice for the Calculation of Daylight Availability, Journal of the Illuminating Engineering Society, vol. 13, nᵒ 4, p. 381‑392, juill. 1984, doi: 10.1080/00994480.1984.10748791.

import numpy as np

def solar_input(julian_date, time_h, latitude_deg, longitude_deg, azimuth_deg, meteo_mode):
    
    # TIME
    J = julian_date # Julian Date => number of the day in the year
    ts = time_h # hours
    
    # METEO
    mode = meteo_mode # 0 clear sky, 1 partly cloudy, 2 cloudy
    # Nébulosité totale = Sky cover method. For meteo mode : 0 to 3 tenths -> mode 0, 4 to 7 tenths -> mode 1, 8 to 10 tenths -> mode 2
    
    # POSITION
    latitude = latitude_deg * np.pi / 180 # in rad
    longitude = longitude_deg * np.pi / 180 # in rad
    ae = azimuth_deg / 180 * np.pi
    
    # PHOTO-CONVERSION
    #1 W/m2 = 126.7 lux
    #1 lux = 0.0185 µmolPhoton/m2/s
    #1 W/m2 = 126.7*0.0185 = 2.3495 µmolPhoton/m2/s 
    # Rendement PAR
    eta = 0.423 # https://www.fondriest.com/environmental-measurements/parameters/weather/photosynthetically-active-radiation/
    
    # CALCULATION
    l = latitude
    L = longitude
    SM = np.floor(longitude * 180 / np.pi / 15) * 15 /180*np.pi # Standard merdian
    ET = 0.170 * np.sin(4 * np.pi / 373 * (J - 80)) - 0.129 * np.sin(2 * np.pi / 355 * (J - 8))
    td = ts + 1 # daylight time in hour
    t = ts + ET + 12 * (SM - L) / np.pi # solar time
    
    delta = 0.4093 * np.sin(2 * np.pi / 368 * (J - 81)) # solar declination in rad
    at = np.arcsin(np.sin(l)*np.sin(delta) - np.cos(l)*np.cos(delta)*np.cos(np.pi * t / 12)) # solar altitude
    ai = -1
    
    # Can we see the sun? 
    if at > 0:
        as_ = np.arctan(np.cos(delta)*np.sin(np.pi * t / 12) / (np.cos(l)*np.sin(delta) + np.sin(l)*np.cos(delta)*np.cos(np.pi * t / 12))) # solar azimuth
        az = as_ - ae # Solar elevation azimuth
        ai = np.arccos(np.cos(at) * np.cos(az)) # incidence angle
        ap = np.arctan(np.sin(at) / np.cos(ai)) # profile angle
        
        # SOLAR FLUX
        Ext = 127.5 * (1 + 0.034 * np.cos(2 * np.pi / 365 * (J - 2))) # Flux solaire extraterrestre en klux
        m = 1 / np.sin(at) # air mass
        
        # OUTPUTS - SKYLIGHT (NOT DIRECT SOLAR)
        Ekv = -1 # Diffuse on a vertical surface ! in klux
        Ekh = -1 # Diffuse on a horizontal surface ! in klux
        c = -1
        # clear sky
        if mode == 0:
            Ekv = (4.0 * at ** 1.3 + 12.0 * np.sin(at)**0.3*np.cos(at)**1.3) * ( (2 + np.cos(az)) / (3 - np.cos(az)))
            Ekh = 8.2 * np.sin(at)**0.5 + 6.9*np.sin(at)*np.cos(at)*np.cos(az)
            c = 0.21
        # partly cloudy
        elif mode == 1:
            Ekv = (12.0 * at + 30.2 * np.sin(at) ** 0.8 * np.cos(at)) * ((1 + np.cos(az)) / (3 - np.cos(az)))
            Ekh = 22.7 * np.sin(at) + 14.1*np.sin(at)**1.3*np.cos(at)*np.cos(az)
            c = 0.8
        # Cloudy sky
        elif mode == 2:
            Ekv = 8.5 * np.sin(at)
            Ekh = 10.7 * np.sin(at)
            c = 1e9
        else:
            print("Lol !!!")
    
        # OUTPUTS - DIRECT SOLAR
        Edn = Ext * np.exp(- c * m) # DNI
        Edh = Edn * np.sin(at) # DNI on an horizontal surface !
        Edv = Edn * np.cos(ai) # DNI on a vertical surface ! 
    else:
        Edn, Edh, Edv, Ekh, Ekv, ai = 0, 0, 0, 0, 0, -1

    return [Edn, Edh, Edv, Ekh, Ekv, at, ai] # in klux / rad for at

def test_routine():
    # TIME
    J = 101 # Julian Date => number of the day in the year
    ts = 10 # hours
    
    # METEO
    mode = 0 # 0 clear sky, 1 partly cloudy, 2 cloudy
    # Nébulosité totale = Sky cover method. For meteo mode : 0 to 3 tenths -> mode 0, 4 to 7 tenths -> mode 1, 8 to 10 tenths -> mode 2
    
    # POSITION
    # Test case from: Recommended Practice for the Calculation of Daylight Availability, Journal of the Illuminating Engineering Society, vol. 13, nᵒ 4, p. 381‑392, juill. 1984, doi: 10.1080/00994480.1984.10748791.
    latitude = 0.68/np.pi*180 # °
    longitude = 1.34/np.pi*180 # °
    azimuth = -90 # azimuth par rapport au sud de la surface verticale, plein sud = 0, est = -90, ouest + 90

    # Test 
    [Edn, Edh, Edv, Ekh, Ekv, at, ai] = solar_input(J, ts, latitude, longitude, azimuth, mode)
    print("A. DNI: {:.1f} klux".format(Edn))
    print("B. Indicent on horizontal: {:.1f} klux".format(Edh))
    print("C. Indicent on vertical: {:.1f} klux".format(Edv))
    print("(B² + C² is not equal to A², as C depends on orientation, while B does not)")
    print("Sky on horizontal: {:.1f} klux".format(Ekh))
    print("Sky on vertical: {:.1f} klux".format(Ekv))

# test_routine()