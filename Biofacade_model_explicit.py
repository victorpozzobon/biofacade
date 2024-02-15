# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
## --- Clear all --- ###
import sys
sys.modules[__name__].__dict__.clear()
## --- --------- --- ###

## --- Importing Python modules --- ###
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os.path

## --- Importing case-specific module --- ###
from UHII_model.UHII_function import UHII_t
from helpFunctions.Biofacade_preprocessing import extracting_and_subsampling
from helpFunctions.Biofacade_reflectivity import transmission_multiple_interface, transmission_multiple_interface_cos_distribution
from helpFunctions.Biofacade_wind import facadeConvectionCoefficient

print('Starting execution ...')
startTime = time.time()
plt.close('all')

def span(array):
    return np.max(array)-np.min(array)

'''
The core function which compute reservoir temperature and averaged illumination over time
'''
def biofacade(x=[0, 0.6, 0.9, 0.9], 
              y=[1, 0, 0.5, 0.05, 1, 0, 0, 0], 
              supDataMeteo = None, 
              supDataLight = None, 
              runID = 0, # to be used for parallel computing
              stationID = -1):
    
    ### --- Overall parameters --- ###
    # Selection station based on its ID
    # Default, stationID = -1 => "Marignane"
    if not stationID > 0:
        station_name = "Marignane" # station_name = "Marignane" #"Orly" #"Prunay" #"Marignane"
    else:
        station_name = str(stationID)

    # Time associated parameter
    transient_day = 0 # number of days to be discarded at the beginning of the process to avoid a bias by a transient period
    dt = 0.5 # h, timestep for differential equation integration
    
    ### --- Processing uncertain parameters (x vector) --- ###
    # Sobol mode - We are operating in pseudo-Sobol mode to account for uncertainty on parameters
    Sobol = bool(x[0]) # Sobol mode
    if Sobol:
        U_in = 0.1 # inner air velocity
        eps_In = x[1] # -, emissivity of the building indoor
        F_sky = 0.5 # -, sky view factor
        eps_MC = x[2] # -, emissivity of the microalgae culture
        eps_Sur = x[3] # -, emissivity of the surrounding buildings
        # T_in = x[6] # -, indoor temperature
        
    ### --- Processing operating modes parameters (y vector) --- ###
    # y = glazing, radiative film, concentration (alpha), thickness, in/out air, boiler, strain, orientation
    
    # 1st - single or double glazing
    n_glaz = int(y[0]) # -, number of glazing for 1. single or 2. double-glazing, or more
        
    # 2nd - radiative film
    film_mode = int(y[1]) # 0 - no film, 1 - greenhouse, 2 - heat management

    # If there is no film, then transmittances are equal to 1
    Tfilm_UV = 1
    Tfilm_Vis = 1
    Tfilm_NIR = 1
    Tfilm_MFIR = 1
    if film_mode == 1: # Based on Nwoba et al., 2019
        Tfilm_UV = 0.2622 # not used
        Tfilm_Vis = 0.8915
        Tfilm_NIR = 0.8315
        Tfilm_MFIR = 0.5114
    if film_mode == 2: # Based on Waché et al., 2020
        Tfilm_UV = 1
        Tfilm_Vis = 0.764
        Tfilm_NIR = 0.429
        Tfilm_MFIR = 0.429
    
    
    # 3rd - culture concentration
    alpha = y[2] # -, transmittend fraction of green light 

    # 4th - culture compartment thickness
    e_MC = y[3] # m, culture reservoir tickness
    
    # 5th - sparged gas from outside or inside
    Tgas_from_building = bool(y[4]) # if True Tgas_in = T_in
    
    # 6th - sparged gas variable temperature
    boilerWinter = bool(y[5]) # if True Tgas_in += 60 °C during winter
    T_boiler = 153.62 # °C, based on Serrano et al., 2013
    FumeDilution = 4 # - 
    
    # 7th - choice of strain
    Sorokiniana = bool(y[6]) # if True go for thermophyllic one.
    
    # Chlorella vulgaris thermal affinity parameters
    T_lower = 17.7 # °C
    T_upper = 39.9 # °C
    # Chlorella sorokiniana thermal affinity parameters
    if Sorokiniana:
        T_lower = 30.6 # °C
        T_upper = 43.2 # °C
    
    # 8th - orientation
    azimuth = y[7] # °, basically, biofacade orientation


    # %%% --- Loading the preprocessed data ---  %%% #
    # Filenames declaration
    folder = './'
    # File containing the illumination information based on the Julian time and position
    illumination_filename = folder + "data/illuminations_vector_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth) + ".csv"
    # File containing the weather monitoring information
    meteo_filename = folder + "data/meteo_vector_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth) + ".csv"
    
    # If files do not exist, preprocess then load
    if (not os.path.isfile(illumination_filename)) and (not os.path.isfile(meteo_filename)) and (supDataMeteo is None):
        print('Creating database ...')
        extracting_and_subsampling(station_name = station_name, dt = dt, azimuth = azimuth, transient_day = transient_day)
    
    # Actual data loading of weather data
    if not Sobol: print('Loading weather database ...') # Suppressed in Sobol mode to avoid chaotic consol spam
    if supDataMeteo is None:
        dataWeather = np.genfromtxt(meteo_filename, delimiter = " ")
    else:
        dataWeather = supDataMeteo # in some case, e.g., parallel run, the weather data can be supplied to avoid mutiple loading
        
    # Dispatching weather related data
    if not Sobol: print('Sorting meteo fields') # Suppressed in Sobol mode to avoid chaotic consol spam
    t_new = dataWeather[:, 0]
    temperature_new = dataWeather[:, 1]
    cloudCover_new = dataWeather[:, 2]
    clear_Sky_Heat_Flux = dataWeather[:, 3]
    UHII_new = dataWeather[:, 4]
    pressure_new = dataWeather[:, 5]
    velWind_new = dataWeather[:, 6]
    dirWind_new= dataWeather[:, 7] 
    if not Sobol: print('Meteo data loaded')
    
    # Actual data loading of illumination data
    if not Sobol: print('Loading illumination database ...') # Suppressed in Sobol mode to avoid chaotic consol spam
    if supDataLight is None:
        illuminations_vector = np.genfromtxt(illumination_filename, delimiter = " ")
    else:
        illuminations_vector = supDataLight # in some case, e.g., parallel run, the illumination data can be supplied to avoid mutiple loading
    
    # Dispatching illumination related data
    if not Sobol: print('Sorting illumination fields')
    at = illuminations_vector[:, 5] # solar altitude in rad
    ai = illuminations_vector[:, 6] # solar incidence angle in rad
    
    if not Sobol: print('Illumination data loaded')
    
    if not Sobol: print('Database loaded')
    if not Sobol: print('Starting to compute biofacade temperature ...')
    if not Sobol: print('Initilization of parameters ...')
    
    ### --- Years identification --- ###
    # it is used to latter dispatch the results between the different years
    first_year = 2013
    current_year = first_year
    n_year = int(t_new[-1] / 365 / 24) # A bit rough but functionnal
    n_year = max(n_year, 1)
    year_time_limit = np.zeros(n_year+1)
    year_time_limit[0] = transient_day * 24
    for i in range(0, n_year):
        year_time_limit[i+1] = year_time_limit[i] + 365 * 24
        # Managing leap year -- 2020 was a leap year
        if (current_year - 2020) % 4 == 0:
            year_time_limit[i+1] += 24
        current_year += 1
    
    ### --- Winters identification --- ###
    # it is used to activate the boiler fumes recirculation (in case the strategy is applied)
    winterTime = np.ones(len(t_new))
    for i in range(int(transient_day * 24 / dt), len(t_new)):
        for j in range(1, len(year_time_limit)):
            if (t_new[i] > year_time_limit[j] -274*24) and (t_new[i] < year_time_limit[j] -60*24): #1er avril -> 1er novembre, no heating
                # Not winter time
                winterTime[i] = 0

    # %%% --- Physical & Operational parameters ---  %%% #
    
    ### Biofacade geometrical parameters & operating conditions ###
    # Geometry
    # e_MC = 0.05 # m, culture reservoir tickness - set as parameter
    H_In = 4 # m, biofacade height
    W_MC = 1 # m, biofacade width
    H_ref = 2.08 # m, reference height for indoor free convection correlation
    L_ref = (4+1)/2. # m, characteristica llegth for indoor forcer convection
    e_air = 0.015 #m, double glazing air layer thickness
    e_PMMA = 0.015 #m, PMMA glazing thickness
    r_Building, r_Meteo, E_Meteo = 0.3, 0.3, 10 # m, wind velocity corrections
    E = 20 # m, biofacade elevation on the host building
    
    # Operation
    T_in = 22 # °C, building indoor temperature
    T_gas = 22 # °C inlet temperature of the sparged gaz
    # alpha = 0.5 # -, transmitted green light - set as parameter
    F = 0.2 # vvm, areation flowrate -> To be converted to per second
    
    ### Biological properties ###
    eta_PS = 0.0434 # -, photosynthetic efficiency
    
    ### Phisical properties ###
    sigma_boltz = 5.67e-8 # W/m2/k^4
    k_air = 0.026 # W/m/K, air thermal conductivity 
    nu_air = 1.85e-5 #m²/s, air kinematic viscosity
    Cp_air = 1004 # J/kg/K, air specific heat
    Pr = 0.71 #-, air Prandtl number
    if not Sobol: U_in = 0.1 # m/s, indoor air velocity - set as parameter
    k_PMMA = 0.18987 # W/m/K, PMMA thermal conductivity 
    rho_water = 1000 # kg/m3
    Cp_water = 4183 # J/kg/K
    rhoCpV = rho_water * Cp_water * H_In * e_MC * W_MC
    
    # Radiative properties 
    f_vis_sun = 0.487 # -, visible fraction of the solar spectrum
    if not Sobol: F_sky = 0.5 # -, sky view factor - set as parameter
    F_sur = 1 - F_sky # -, surrounding view factor
    if not Sobol: eps_MC = 0.9 # -, emissivity of the microalgae culture - set as parameter
    if not Sobol: eps_Sur = 0.9 # -, emissivity of the surrounding buildings - set as parameter
    if not Sobol: eps_In = 0.6 # -, emissivity of the building indoor - set as parameter
    
    # Refraction indexes @ 400 nm
    n_air = 1.00028276 # https://refractiveindex.info/?shelf=other&book=air&page=Ciddor
    n_water = 1.3390 # https://refractiveindex.info/?shelf=main&book=H2O&page=Hale
    n_PMMA = 1.5082 # https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Zhang-Mitsubishi
    
    n_s_indoor = [n_air, n_PMMA, n_water] # Successive interfaces for the indoor
    # Successive interfaces for the outdoor
    n_s_outdoor = []
    # Adding glazing 
    for i in range(0, n_glaz):
        n_s_outdoor.append(n_air)
        n_s_outdoor.append(n_PMMA)
    n_s_outdoor.append(n_water)
    
    # Computing fixed transmittance
    # Inward flux
    tau_sun = transmission_multiple_interface(0, n_s_outdoor)
    tau_sky_abs = transmission_multiple_interface_cos_distribution(n_s_outdoor)
    tau_sur_abs = transmission_multiple_interface_cos_distribution(n_s_outdoor)
    tau_in_abs = transmission_multiple_interface_cos_distribution(n_s_indoor)
    # Outward flux
    n_s_indoor.reverse()
    n_s_outdoor.reverse()
    tau_sky_emi = transmission_multiple_interface_cos_distribution(n_s_outdoor)
    tau_sur_emi = transmission_multiple_interface_cos_distribution(n_s_outdoor)
    tau_in_emi = transmission_multiple_interface_cos_distribution(n_s_indoor)
    n_s_outdoor.reverse() # Reset to correct value
    
    # %%% --- Numerical allocation & Initialization ---  %%% #
    # Some declarion and preallocation
    Results_Matrix = np.zeros([len(t_new), 22])
    T_bio = np.zeros(len(t_new))
    UHII_storage = np.zeros(len(t_new))
    
    # Initialization
    T_bio[0] = T_in
    dayTime = False # flag to classify the timesteps
    UHII_Daily_night_to_come, UHII_Daily_before, t_sundawn = 0, 0, 0 
    
    # Some variable to use resistance in series to find 'skin' temperatures of outermost glazing
    TPMMAIn = T_in
    TPMMAInOld = TPMMAIn + 1 
    TPMMAOut = temperature_new[0]-273.15
    TPMMAOutOld = TPMMAOut + 1
    h_conv_save = 0 # maximum outdoor convection coefficient
    
  
    # %%% --- Computing ---  %%% #
    if not Sobol: print('Time loop ...')
    for i in range(0, len(t_new)-1):
        # Some display
        if (runID == 0) and (i % 10000 == 0): print("Iteration {:d} / {:d}".format(i, len(t_new)))
        
        # To avoid to many memory access and too long formulas
        CC = cloudCover_new[i]
        Tair_Out = temperature_new[i]
        
        # Strategies on gas temperature
        if Tgas_from_building:
            T_gas = T_in
        else:
            T_gas = Tair_Out - 273.15
            
        if boilerWinter and winterTime[i] > 0:
            T_gas = ((FumeDilution - 1)*T_gas + T_boiler) / FumeDilution # Double validated formula, based on mass and energy balance
        
        # Direct solar 
        Phi_Sun_Tot = illuminations_vector[i, 2] * 1000 / 126.7 # In W/m²
        Phi_Sun_Vis = (3-2*alpha**2-alpha)/3. * (1-eta_PS) * f_vis_sun * Phi_Sun_Tot
        Phi_Sun_IR = (1 - f_vis_sun) * Phi_Sun_Tot
        if ai[i] >= 0:
            tau_sun = transmission_multiple_interface(ai[i], n_s_outdoor)
        else:
            tau_sun = 0
        
        Phi_Sun_Abs = tau_sun * (Tfilm_Vis*Phi_Sun_Vis + Tfilm_NIR*Phi_Sun_IR)
        
        # Sky radiative exchange
        Phi_Sky_Vis = illuminations_vector[i, 4] * 1000 / 126.7 # In W/m²
        Phi_Sky_Tot = 0
        Eps_skyT_sky4 = (9.36578e-6 * (1-CC) * Tair_Out**6 
                + Tair_Out**4 * CC * 
                                ( 
                                (1-0.84*CC) * (0.527 + 0.161 * np.exp(8.45 * ( 1 - 273/Tair_Out )))
                                + 0.84 * CC
                                )
                                )
        Phi_Sky_Tot = sigma_boltz * Eps_skyT_sky4
        Phi_Sky_Abs = F_sky * (0.4 * tau_sun + 0.6 * tau_sky_abs) * ((3-2*alpha**2-alpha)/3. * (1-eta_PS) * Tfilm_Vis * Phi_Sky_Vis + Tfilm_MFIR * (Phi_Sky_Tot-Phi_Sky_Vis))
        Phi_Sky_Emi = Tfilm_MFIR * F_sky * tau_sky_emi * eps_MC * sigma_boltz * (T_bio[i] + 273.15)**4
        
        # Surrounding radiative exchange
        if not dayTime and at[i] >= 0:
            dayTime = True
            UHII_Daily_before = UHII_Daily_night_to_come
            UHII_Daily_night_to_come = UHII_new[i]
            t_sundawn = t_new[i]
        if dayTime and at[i] < 0:
            dayTime = False
        
        UHII_ti = UHII_t(UHII_Daily_night_to_come, UHII_Daily_before, t_new[i], t_sundawn, basal_level = 0.5)
        UHII_storage[i] = UHII_ti
        Trur = (1-CC) * (2.82 + 1.15 * (Tair_Out-273.15)) + CC * (1.33 + 1.00 * (Tair_Out-273.15))
        Tsur = Trur + UHII_ti
        
        
        Phi_Sur_Abs = Tfilm_MFIR * F_sur * tau_sur_abs * eps_Sur * sigma_boltz * (Tsur + 273.15)**4
        Phi_Sur_Emi = Tfilm_MFIR * F_sur * tau_sur_emi * eps_MC * sigma_boltz * (T_bio[i] + 273.15)**4
        
        # Indoor radiation
        Phi_In_Rad_Abs = Tfilm_MFIR * tau_in_abs * eps_In * sigma_boltz * (T_in + 273.15)**4
        Phi_In_Rad_Emi = Tfilm_MFIR * tau_in_emi * eps_MC * sigma_boltz * (T_bio[i] + 273.15)**4
        
        # Indoor convection
        # Loop to find TPMMAIn (indoor PMMA skin temperature)
        while (abs(TPMMAInOld - TPMMAIn) > 0.01):
            h_in_conv_free = 2.04 * (H_In/H_ref * abs(TPMMAIn - T_in))**0.23
            
            Re_Ref = U_in * L_ref / nu_air
            h_in_conv_forced = k_air / L_ref * 0.664 * Re_Ref **0.5 * Pr ** (1./3.)
            
            Phi_conv_net_in = (T_in - T_bio[i]) / ( \
                                                  1/np.max([h_in_conv_forced, h_in_conv_free, 1e-6]) + 
                                                  e_PMMA / k_PMMA)
            TPMMAInOld = TPMMAIn
            TPMMAIn = max(min(T_in - Phi_conv_net_in * (e_PMMA / k_PMMA), TPMMAInOld +1), TPMMAInOld -1)
        TPMMAInOld = TPMMAIn + 1. # to initiate for the next round
        
        # Outdoor convection
        # Loop to find TPMMAOut (outdoor PMMA skin temperature)
        while (abs(TPMMAOutOld - TPMMAOut) > 0.01):
            rho_air = pressure_new[i] * 0.029 / 8.315 / (Tair_Out)
            # Differentiating to get local value
            RaH_up = 9.81 * (1/(Tair_Out) * abs(TPMMAOut - (Tair_Out-273.15)) * (E+H_In/2.)**3) / nu_air / (k_air / Cp_air / rho_air)
            h_out_conv_free_up = k_air / (E+H_In/2.) * (0.825 + 0.387 * RaH_up**(1./6.) / (1 + (0.492*RaH_up)**(9./16.))**(8./27.) )
            RaH_down = 9.81 * (1/(Tair_Out) * abs(TPMMAOut - (Tair_Out-273.15)) * (E-H_In/2.)**3) / nu_air / (k_air / Cp_air / rho_air)
            h_out_conv_free_down = k_air / (E-H_In/2.) * (0.825 + 0.387 * RaH_down**(1./6.) / (1 + (0.492*RaH_down)**(9./16.))**(8./27.) )
            h_out_conv_free = ((E+H_In/2.)*h_out_conv_free_up - (E-H_In/2.)*h_out_conv_free_down) / H_In
            
            
            U_Meteo = velWind_new[i] # in m/s
            # azimuth = facade orientation # azimuth par rapport au sud de la surface verticale, plein sud = 0, est = -90, ouest + 90
            # dirWind_new = wind origin # in degree. Wind 0° => coming from the North, then, clock rotation.
            WindIncidenceAngle = abs(dirWind_new[i] - (azimuth+180.))
            h_out_conv_forced = facadeConvectionCoefficient(U_Meteo, WindIncidenceAngle, r_Building, r_Meteo, E, E_Meteo)
            if h_conv_save < h_out_conv_forced: h_conv_save = h_out_conv_forced
            Phi_conv_net_out = ((Tair_Out-273.15) - T_bio[i]) / ( \
                                                  1/np.max([h_out_conv_forced, h_out_conv_free, 1e-6]) + 
                                                  n_glaz * e_PMMA / k_PMMA + (n_glaz-1) * e_air / k_air)
            TPMMAOutOld = TPMMAOut
            TPMMAOut = max(min(T_bio[i] - Phi_conv_net_out * (n_glaz * e_PMMA / k_PMMA + (n_glaz-1) * e_air / k_air), TPMMAOutOld +1), TPMMAOutOld -1)
        TPMMAOutOld = TPMMAOut + 1. # to initiate for the next round
        
        # Gas flow power
        P_gas = F / 60 * e_MC * H_In * W_MC  * Cp_air * (T_gas - T_bio[i])
        

        # Integration    
        T_bio[i+1] = T_bio[i] + dt * 3600 / rhoCpV * ( \
                                                      H_In * W_MC * ( \
                                                               Phi_conv_net_in + Phi_conv_net_out + Phi_In_Rad_Abs -
                                                               Phi_In_Rad_Emi + Phi_Sky_Abs - Phi_Sky_Emi +
                                                               Phi_Sur_Abs - Phi_Sur_Emi + Phi_Sun_Abs) +
                                                          P_gas )
            
        
        # Individual ligth fluxes
        visible_light = (tau_sun * illuminations_vector[i, 2] + (0.4 * tau_sun + 0.6 * tau_sky_abs) * illuminations_vector[i, 4]) * 0.0185 * 1000 # in µmolPhoton/m²/s
        visible_light = visible_light * Tfilm_Vis
        visible_light_averaged = visible_light / np.log(alpha) * (alpha + alpha**2 - 2) / 3
        
        # Storing results
        '''
        0  - Time, h
        1  - Outdoor air temperature, °C
        2  - Reservoir temperture, °C
        3  - Net convective heat flux with the indoor, W/m2
        4  - Net convective heat flux with the outdoor, W/m2
        5  - Absorbed radiative heat flux from the indoor, W/m2
        6  - Emitted radiative heat flux toward the indoor, W/m2
        7  - Absorbed radiative heat flux from the sky, W/m2
        8  - Emitted radiative heat flux toward the sky, W/m2
        9  - Absorbed radiative heat flux from the surrounding, W/m2
        10 - Emitted radiative heat flux toward the surrounding, W/m2
        11 - Absorbed radiative heat flux from the sun, W/m2
        12 - Net power exhanged with sparging gas, W
        13 - Skin temperature of PMMA, indoor side, °C
        14 - Skin temperature of PMMA, outdoor side, °C
        15 - Clear sky heat flux, W/m2 (not used)
        16 - Incident visible light, µmolPhoton/m²/s
        17 - Volume averaged visible light, µmolPhoton/m²/s
        18 - Urban Heat Insland Intensity, K
        19 - Solar altitude, rad
        20 - Solar incidence, rad
        21 - Surrounding apparent temperature, °C
        '''

        Results_Matrix[i, : ] = [t_new[i], Tair_Out-273.15, T_bio[i], Phi_conv_net_in, Phi_conv_net_out, \
            Phi_In_Rad_Abs, Phi_In_Rad_Emi, Phi_Sky_Abs, Phi_Sky_Emi, Phi_Sur_Abs, Phi_Sur_Emi, Phi_Sun_Abs, P_gas / (H_In * W_MC), \
                TPMMAIn, TPMMAOut, clear_Sky_Heat_Flux[i], visible_light, visible_light_averaged, UHII_ti, at[i], ai[i], Tsur]

    # Clean up    
    Results_Matrix = np.array(Results_Matrix)   
    
    # %%% --- Plotting the results ---  %%% #
    if not Sobol:
        ### --- External & biofacade temperature --- ###
        ax = plt.subplot(2, 2, 1)
        ax.plot(t_new, temperature_new-273.15, alpha = 0.5, label = "Outdoor temperature")
        ax.plot(t_new, T_bio, alpha = 0.5, label = "Biofacade temperature")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Temperature (°C)")
        ax.set_ylim([-10, 50])
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
        
        ### --- UHII --- ###
        ax = plt.subplot(2, 2, 2)
        ax.plot(t_new, Results_Matrix[:, 18], label = "UHII")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("UHII (°C)")
        ax.set_ylim([-5, +5])
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
        
        ### --- Clear sky heat flux on horizontal surface -> determine sundawn --- ###
        ax = plt.subplot(2, 2, 3)
        ax.plot(t_new, clear_Sky_Heat_Flux, label = "Inicident direct heat flux on an horizontal surface")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Heat flux (W/m²)")
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
        
        ### --- Absorbed direct sun flux on the biofacade --- ###
        ax = plt.subplot(2, 2, 4)
        ax.plot(t_new, Results_Matrix[:, 7], label = "Absorbed heat flux from the sun")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Heat flux (W/m²)")
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
    
    # %%% --- Deeper analysis ---  %%% #
    
    ### --- Plotting each quantity of interest --- ###
    name = ["Time", "Tair_Out", "T_bio", "Phi_conv_net_in", "Phi_conv_net_out", \
        "Phi_In_Rad_Abs", "Phi_In_Rad_Emi", "Phi_Sky_Abs", "Phi_Sky_Emi", "Phi_Sur_Abs", "Phi_Sur_Emi", "Phi_Sun_Abs", "P_gas",
        "TPMMAIn", "TPMMAOut", "clear_Sky_Heat_Flux", "visible_light", "visible_light_averaged", "UHII_t", "at", "ai", "Tsur"]
        
    # %%% --- Plotting the metrics ---  %%% #
    
    ### --- Converting results to a dataframe to compute performances --- ###
    df_results_whole = pd.DataFrame(data = Results_Matrix, columns = name)
    harvest = np.zeros([n_year, 6+6+1])
    harvest_names = ["LNok_Tlow_h", "LNok_Tok_h", "LNok_Thigh_h", "Lok_Tlow_h", "Lok_Tok_h", "Lok_Thigh_h", 
                     "T_below_0_h", "TotalProductible_h", "10DecValue_C", "10DecSpan_C", "1DecValue_C", "1DecSpan_C",
                     "year"]
    
    for yearIndex in range(0, n_year):
        df_results = df_results_whole[(df_results_whole["Time"] > year_time_limit[yearIndex]) 
                                      & (df_results_whole["Time"] <= year_time_limit[yearIndex+1])].copy(deep=True)
        time_below_0 = np.sum((df_results["T_bio"] < 0)) * dt

        ### Extracting daytime -> 
        df_productible = df_results[(df_results["at"] > 0)]
        if not Sobol: print("Number of timestep with at > 0: {:d}".format(len(df_productible)))
        ### Differentiating between light-sufficient and light deficient ->
        df_light_ok = df_productible[df_productible["visible_light_averaged"] >= 150]
        df_light_Nok = df_productible[df_productible["visible_light_averaged"] < 150]
        
        ### Differentiating between the three temperature range ->
        # T_lower = 17.7 # °C
        # T_upper = 39.9 # °C
        
        harvest[yearIndex, :6] = dt * np.array([np.sum((df_light_Nok["T_bio"] < T_lower)), 
                                            np.sum((df_light_Nok["T_bio"] > T_lower) & (df_light_Nok["T_bio"] <= T_upper)), 
                                            np.sum((df_light_Nok["T_bio"] > T_upper)),
                                            np.sum((df_light_ok["T_bio"] < T_lower)),
                                            np.sum((df_light_ok["T_bio"] > T_lower) & (df_light_ok["T_bio"] <= T_upper)),
                                            np.sum((df_light_ok["T_bio"] > T_upper))])
        
        ### Computing percentages ->
        Total_productible = len(df_productible) * dt
        harvest[yearIndex, 6] = time_below_0
        harvest[yearIndex, 7] = Total_productible
        
        ### Decile analysis ->
        df_results['Decile_rank'] = pd.qcut(df_results['T_bio'], 10, labels = False)
        # Decile analysis -- too hot
        harvest[yearIndex, 8] = np.average(df_results[df_results["Decile_rank"] == 9]['T_bio'])
        harvest[yearIndex, 9] = span(df_results[df_results["Decile_rank"] == 9]['T_bio'])
        # Decile analysis -- too cold
        harvest[yearIndex, 10] = np.average(df_results[df_results["Decile_rank"] == 0]['T_bio'])
        harvest[yearIndex, 11] = span(df_results[df_results["Decile_rank"] == 0]['T_bio'])
        
        harvest[yearIndex, 12] = yearIndex + first_year
        
        if not Sobol: print("Time below 0°C: {:.2f} h". format(time_below_0))
        
        if not Sobol:
            np.save("results_vector_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth), harvest)
            np.save("results_vector_extended_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth), Results_Matrix)
        
        if not Sobol: print("Maximum outdoor convection coefficient: {:.2f} W/m²/K".format(h_conv_save))
        if not Sobol: print("Average temperature: {:.2f} °C".format(np.average(T_bio)))
    
    return harvest


    if not Sobol: print('Normal end of execution.\nExecution time: ' + format(time.time() - startTime, '.2f') + ' s')

'''
x vector
0 - sobol mode: 0 - disabled, 1 - enabled
1 - emissivity of the building indoor: [0, 1]
2 - emissivity of the microalgae culture: [0, 1]
3 - emissivity of the surrounding buildings: [0, 1]

y vector
0 - nb. of front glazing: 1 (single glazing), 2 (double glazing), or more
1 - film_mode: 0 - no film, 1 - greenhouse, 2 - heat management
2 - culture concentration, as transmittend fraction of green light: [0, 1]
3 - culture compartment thickness (m), >0
4 - sparged gas from outside or inside: 0 - temperature from the outdoor, 1 - temperature from the building
5 - sparged gas variable temperature (K): 0 - no change over the year, 1 - +60°C in winter
6 - choice of strain: 0 - Chlorella vulgaris (mesophyllic), 1 - Chlorella sorokiniana (thermophyllic)
7 - azimuth (°), basically, biofacade orientation: [-90, 90]
'''
resultsArray = biofacade(x=[0, 0.6, 0.9, 0.9], 
              y=[1, 0, 0.5, 0.05, 1, 0, 0, 0], 
              supDataMeteo = None, 
              supDataLight = None, 
              runID = 0, 
              stationID = -1)

for year in range(0, resultsArray.shape[0]):
    current_year = resultsArray[year, :]
    print("--- Year {:d} ---".format(int(current_year[-1])))
    print("Light Nok T low: {:.1f} h".format(current_year[0]))
    print("Light Nok T ok: {:.1f} h".format(current_year[1]))
    print("Light Nok T high: {:.1f} h".format(current_year[2]))
    print("Light ok T low: {:.1f} h".format(current_year[3]))
    print("Light ok T ok: {:.1f} h".format(current_year[4]))
    print("Light Nk T high: {:.1f} h".format(current_year[5]))
    print("Time below 0 °C: {:.1f} h".format(current_year[6]))
    print("Total productible time: {:.1f} h".format(current_year[7]))
    print("Average of the 9th decile: {:.1f} °C".format(current_year[8]))
    print("Span of the 9th decile: {:.1f} °C".format(current_year[9]))
    print("Average of the 1st decile: {:.1f} °C".format(current_year[10]))
    print("Span of the 1st decile: {:.1f} °C".format(current_year[11]))
    print("")
