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
from helpFunctions.Biofacade_biomodulation import muOfT, PigEquilibrium

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
              z=[0.95, 3, 1],
              supDataMeteo = None, 
              supDataLight = None, 
              runID = 0, # to be used for parallel computing
              stationID = -1,
              transient_day = 41,
              dt = 0.5):
    
    '''###############################
    ### --- Overall parameters --- ###
    ###############################'''
    # Selection station based on its ID
    # Default, stationID = -1 => "Marignane"
    if not stationID > 0:
        station_name = "Marignane" # station_name = "Marignane" #"Orly" #"Prunay" #"Marignane"
    else:
        station_name = str(stationID)

    # Time associated parameter
    # With 41 for year 2023, Starts on Monday, 21 November 2022, 1am
    #transient_day = 41 # number of days to be discarded at the beginning of the process to avoid a bias by a transient period
    #dt = 0.5 # h, timestep for differential equation integration
    
    '''#######################################################
    ### --- Processing uncertain parameters (x vector) --- ###
    #######################################################'''
    
    # Sobol mode - We are operating in pseudo-Sobol mode to account for uncertainty on parameters
    Sobol = bool(x[0]) # Sobol mode
    if Sobol:
        U_in = 0.1 # inner air velocity
        eps_In = x[1] # -, emissivity of the building indoor
        F_sky = 0.5 # -, sky view factor
        eps_MC = x[2] # -, emissivity of the microalgae culture
        eps_Sur = x[3] # -, emissivity of the surrounding buildings
        # T_in = x[6] # -, indoor temperature
        
    '''#############################################################
    ### --- Processing operating modes parameters (y vector) --- ###
    #############################################################'''
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

    '''########################################################
    ### --- Processing bioprocess parameters (z vector) --- ###
    ########################################################'''
    # z = dilution, 
    
    # 1st - fraction of the reservoir volume replaced when the transmitted light reaches the limit value
    dilution = z[0] # -, [0, 1[
    
    # 2sd - maintenance mode. 0 - null, 1 - constant light respiration, 2 - constant dark respiration, 3 - time varying
    maintenanceMode = int(z[1])
    
    # 3rd - temperature mode. 0 - no temperature dependency, 1 - biological processes are modulated by temperature
    temperatureMode = int(z[2])
    
    # 4th - retroaction of the pigment content on cross section. 0 - variation but no retroaction on cross section, 1 - variation and retroaction
    chlCrossSection = int(z[3])
    
    # 5th - cross sections wavelength band dependency. 0 - no, 1 - yes
    wavelengthBandMode = int(z[4])
    
    '''###############################################
    # %%% --- Loading the preprocessed data ---  %%% #
    ###############################################'''
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
    
    '''#################################
    ### --- Years identification --- ###
    #################################'''
    # it is used to latter dispatch the results between the different years
    first_year = 2023
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
    
    '''###################################
    ### --- Winters identification --- ###
    ###################################'''
    # it is used to activate the boiler fumes recirculation (in case the strategy is applied)
    winterTime = np.ones(len(t_new))
    for i in range(int(transient_day * 24 / dt), len(t_new)):
        for j in range(1, len(year_time_limit)):
            if (t_new[i] > year_time_limit[j] -274*24) and (t_new[i] < year_time_limit[j] -60*24): #1er avril -> 1er novembre, no heating
                # Not winter time
                winterTime[i] = 0

    '''###################################################
    # %%% --- Physical & Operational parameters ---  %%% #
    ###################################################'''
    
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
    I_building = 323 * 4.61/5.98 # lux, illumination coming from indoor lighting - 4.61/5.98 is to adapt it from cold white to solar spectrum
    
    ### Biological properties ###
    eta_PS_0 = 0.0434 # -, maximum photosynthetic efficiency
    eta_PS = 0.0 # -, photosynthetic efficiency at time 0 (in the dark)
    sigma_abs = 692 # m2/kg, absorption cross-section
    sigma_abs_red = 450 # m2/kg, absorption cross-section
    sigma_abs_green = 100 # m2/kg, absorption cross-section
    sigma_abs_blue = 600 # m2/kg, absorption cross-section
    # Matrix storing values:    Whole   Blue    Green   Red
    sigma_matrix = np.array([  [253.,   480.,   110.,   185.  ],  # Dry matter basis
                               [16.5e3, 31.5e3, 7.31e3, 12.0e3]]) # Total chloropĥyll basis
    mu_max = 2.17/24 # 1/h, maximum growth rate
    I_ref = 48.3 # µmolPhoton/m²/s, scaling of photosinthetic efficiency
    HHV = 20.51e6 # J/kg -- Chlorella vulgaris HHV Oliver et al. 
    me_light = 3.88e-2/24 # light respisration h-1, from carbon capture publication
    me_dark = 7.42e-3/24 # dark respiration h-1, from carbon capture publication
    tau_pig = 8.96 # h, characteristic time of pigment content modulation
    pigmentLetterCode = ["a", "b", "l"] # -, coded letter to identify pigments
    
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
    
    '''#######################################################
    # %%% --- Numerical allocation & Initialization ---  %%% #
    #######################################################'''
    # Some declarion and preallocation
    Results_Matrix = np.zeros([len(t_new), 33])
    T_bio = np.zeros(len(t_new))
    X = np.zeros(len(t_new))
    Pigments = np.zeros([len(t_new), 3])
    UHII_storage = np.zeros(len(t_new))
    
    # Initialization
    T_bio[0] = T_in
    X[0] = 0.001 # g/L <=> kg/m3
    # TODO justify values
    Pigments[0, :] = np.array([15, 5, 2]) # mg/g, Chlorophyll a, b, and lutein 
    dayTime = False # flag to classify the timesteps
    UHII_Daily_night_to_come, UHII_Daily_before, t_sundawn = 0, 0, 0 
    
    # Some variable to use resistance in series to find 'skin' temperatures of outermost glazing
    TPMMAIn = T_in
    TPMMAInOld = TPMMAIn + 1 
    TPMMAOut = temperature_new[0]-273.15
    TPMMAOutOld = TPMMAOut + 1
    h_conv_save = 0 # maximum outdoor convection coefficient
    
    '''###########################
    # %%% --- Computing ---  %%% #
    ###########################'''
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
        
        # Computing sun transmitivity
        if ai[i] >= 0:
            tau_sun = transmission_multiple_interface(ai[i], n_s_outdoor)
        else:
            tau_sun = 0
            
        # Updating cross section
        # Please note that when all the cross sections have the same value
        # it means that there is NO wavelength band dependency
        
        if chlCrossSection == 0: # Dry matter basis
            if wavelengthBandMode == 0: # No spectral dependency. One value for the whole visible spectrum
                sigma_abs_red = sigma_matrix[0, 0]
                sigma_abs_green = sigma_abs_red
                sigma_abs_blue = sigma_abs_red
                
            elif wavelengthBandMode == 1: # Based on red, green, blue bands
                sigma_abs_red = sigma_matrix[0, 3]
                sigma_abs_green = sigma_matrix[0, 2]
                sigma_abs_blue = sigma_matrix[0, 1]
                
        elif chlCrossSection == 1: # Total chlorophyll basis
            chlTot = Pigments[i, 0] + Pigments[i, 1] # Total chlorophyll content, mg/g
            
            if wavelengthBandMode == 0: # No spectral dependency. One value for the whole visible spectrum
                sigma_abs_red = sigma_matrix[1, 0] * chlTot / 1000 # because cross section are in m²/kg of chl. 
                sigma_abs_green = sigma_abs_red
                sigma_abs_blue = sigma_abs_red
                
            elif wavelengthBandMode == 1: # Based on red, green, blue bands
                sigma_abs_red = sigma_matrix[1, 3] * chlTot / 1000. # because cross section are in m²/kg of chl. 
                sigma_abs_green = sigma_matrix[1, 2] * chlTot / 1000.
                sigma_abs_blue = sigma_matrix[1, 1] * chlTot / 1000.
            
        # Illumination -> to get eta_PS
        visible_light = (tau_sun * illuminations_vector[i, 2] + (0.4 * tau_sun + 0.6 * tau_sky_abs) * illuminations_vector[i, 4]) # in klux
        visible_light = visible_light * Tfilm_Vis
        
        # Indoor illumination
        # FROM 7 am to 21 pm, 5 days a week, starts on Monday, 21 November 2022, 1am
        if ((t_new[i] +1)%24 > 7) and ((t_new[i] +1)%24 < 21) and (int((t_new[i] +1)/24)%7 < 5): 
            visible_light = visible_light + I_building/1000 # in klux
        
        
        P_Abs = H_In * W_MC * visible_light * (3 - np.exp(- sigma_abs_red * X[i] * e_MC)
                                               - np.exp(- sigma_abs_green * X[i] * e_MC)
                                               - np.exp(- sigma_abs_blue * X[i] * e_MC)) / 3 # in klux
        P_Abs = P_Abs * 1000 / 126.7 # In W
        
        # I_average, volume-averaged illumination
        I_average = visible_light * ( (1 - np.exp(- sigma_abs_red * X[i] * e_MC)) / (sigma_abs_red * X[i] * e_MC)
                                      +(1 - np.exp(- sigma_abs_green * X[i] * e_MC)) / (sigma_abs_green * X[i] * e_MC)
                                      +(1 - np.exp(- sigma_abs_blue * X[i] * e_MC)) / (sigma_abs_blue * X[i] * e_MC)) / 3 # in klux
        I_average = I_average * 1000 * 0.0185 # in µmolPhoton/m²/s
        
        if I_average > 0:
            eta_PS = eta_PS_0 * 4 * I_ref / I_average * ( 1 / ( 1 + np.exp(- I_average / I_ref) ) -0.5)
            regulated_Fraction = 1 - 4 * I_ref / I_average * ( 1 / ( 1 + np.exp(- I_average / I_ref) ) -0.5)
        else:
            eta_PS = 0
            regulated_Fraction = 0

        # Direct solar 
        Phi_Sun_Tot = illuminations_vector[i, 2] * 1000 / 126.7 # In W/m²
        Phi_Sun_Vis = (1-eta_PS) * f_vis_sun * Phi_Sun_Tot * (3 - np.exp(- sigma_abs_red * X[i] * e_MC)
                                                              - np.exp(- sigma_abs_green * X[i] * e_MC)
                                                              - np.exp(- sigma_abs_blue * X[i] * e_MC)) / 3 # In W/m²
        Phi_Sun_IR = (1 - f_vis_sun) * Phi_Sun_Tot # In W/m²
        Phi_Sun_Abs = tau_sun * (Tfilm_Vis*Phi_Sun_Vis + Tfilm_NIR*Phi_Sun_IR) # In W/m²
        
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
        
        ### Temperature ###
        # Integration
        T_bio[i+1] = T_bio[i] + dt * 3600 / rhoCpV * ( \
                                                      H_In * W_MC * ( \
                                                               Phi_conv_net_in + Phi_conv_net_out + Phi_In_Rad_Abs -
                                                               Phi_In_Rad_Emi + Phi_Sky_Abs - Phi_Sky_Emi +
                                                               Phi_Sur_Abs - Phi_Sur_Emi + Phi_Sun_Abs) +
                                                          P_gas )
        
        ### Biomass ###
        
        
        # Tempreraure effect
        if temperatureMode == 1:
            metabolism_scaling = muOfT(Sorokiniana=Sorokiniana, T=T_bio[i])
        else:
            metabolism_scaling = 1.

        # Maintenance
        if maintenanceMode == 0:
            me = 0
        elif maintenanceMode == 1:
            me = me_light
        elif maintenanceMode == 2:
            me = me_dark
        elif maintenanceMode == 3:
            if at[i] > 0: # To avoid indoor to induce "always on" (fake day) culture
                me = me_light
            else:
                me = me_dark
        else:
            print("Problem with maintenance mode. Value = " + str(maintenanceMode))   
        
        
        # Intergration of biomass
        X[i+1] = X[i] + dt * metabolism_scaling * (
            3600 * P_Abs / HHV * eta_PS / (H_In * W_MC * e_MC) 
            - me * X[i])
        
        # Integration of pigments
        # Chl a, b and lutein
        for pig in range(0, len(pigmentLetterCode)):
            Pigments[i+1, pig] = Pigments[i, pig] + dt * metabolism_scaling * 1/tau_pig * (PigEquilibrium(pigmentLetterCode[pig], I_average) - Pigments[i, pig])
            
            
        # Dilution
        outGreenLightFraction = np.exp(- sigma_abs_green * X[i+1] * e_MC) * 100
        if (outGreenLightFraction < alpha * 100) and (at[i] > 0):
            DV = dilution * H_In * W_MC * e_MC
            prod_MC = dilution * H_In * W_MC * e_MC * X[i+1]
        else:
            DV = 0
            prod_MC = 0 
        X[i+1] = (H_In * W_MC * e_MC - DV) / (H_In * W_MC * e_MC) * X[i+1]
            
        
        # DECREPATED - ONLY IN BUILDING MODE
        # KEPT FOR RETROCOMPATIBILITY REASON (ALLEDGLY ...)
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
        22 - Power absorbed by the culutre, W
        23 - Average light intensity within the culture, µmolPhoton/m²/s
        24 - Photosynthetic efficiency, -
        25 - Biomass concentration, g/l
        26 - Fraction of green light transmitted, %
        27 - Metabolism scalling, %
        28 - Microalgae production, kg
        29 - Chlorophyll a content, mg/g
        30 - Chlorophyll b content, mg/g
        31 - Lutein content, mg/g
        32 - Deviation from maximal photosynthetic efficiency, -
        '''

        Results_Matrix[i, : ] = [t_new[i], Tair_Out-273.15, T_bio[i], Phi_conv_net_in, Phi_conv_net_out, \
            Phi_In_Rad_Abs, Phi_In_Rad_Emi, Phi_Sky_Abs, Phi_Sky_Emi, Phi_Sur_Abs, Phi_Sur_Emi, Phi_Sun_Abs, P_gas / (H_In * W_MC), \
                TPMMAIn, TPMMAOut, clear_Sky_Heat_Flux[i], visible_light, visible_light_averaged, UHII_ti, at[i], ai[i], Tsur,
                P_Abs, I_average, eta_PS, X[i], outGreenLightFraction, metabolism_scaling*100, prod_MC, 
                Pigments[i, 0], Pigments[i, 1], Pigments[i, 2], regulated_Fraction]

    # Clean up
    Results_Matrix[-1, : ] = Results_Matrix[-2, : ]
    Results_Matrix = np.array(Results_Matrix)   
    
    '''######################################
    # %%% --- Plotting the results ---  %%% #
    ######################################'''
    if not Sobol:
        idx2023 = int(transient_day*24/dt)
        ### --- Absorbed power --- ###
        ax = plt.subplot(2, 2, 1)
        ax.plot(t_new[idx2023:], Results_Matrix[idx2023:, 22], label = "Absorbed power")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Absorbed power (W)")
        # ax.set_ylim([-10, 50])
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
        
        ### --- X --- ###
        ax = plt.subplot(2, 2, 2)
        ax.plot(t_new[idx2023:], Results_Matrix[idx2023:, 25], label = "X")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("X (g/L)")
        # ax.set_ylim([0, 15])
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
        
        ### --- eta_PS --- ###
        ax = plt.subplot(2, 2, 3)
        ax.plot(t_new[idx2023:], Results_Matrix[idx2023:, 27], label = "Metabolism scaling")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Metabolism scaling (%)")
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
        
        ### --- I_average --- ###
        ax = plt.subplot(2, 2, 4)
        ax.plot(t_new[idx2023:], Results_Matrix[idx2023:, 29], label = "Chlorophyll a")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Chlorophyll a (mg/g)")
        plt.tight_layout(pad=0.5)
        plt.legend()
        plt.show()
    
    '''#################################
    # %%% --- Deeper analysis ---  %%% #
    #################################'''
    
    ### --- Plotting each quantity of interest --- ###
    name = ["Time", "Tair_Out", "T_bio", "Phi_conv_net_in", "Phi_conv_net_out", \
        "Phi_In_Rad_Abs", "Phi_In_Rad_Emi", "Phi_Sky_Abs", "Phi_Sky_Emi", "Phi_Sur_Abs", "Phi_Sur_Emi", "Phi_Sun_Abs", "P_gas",
        "TPMMAIn", "TPMMAOut", "clear_Sky_Heat_Flux", "visible_light", "visible_light_averaged", "UHII_t", "at", "ai", "Tsur",
        "P_Abs", "I_average", "eta_PS", "X", "outGreenLightFraction", "metabolism_scaling", "prod_MC", "Chl_a", "Chl_b", "Lutein",
        "regulated_Fraction"]
        
    # %%% --- Plotting the metrics ---  %%% #
    
    ### --- Converting results to a dataframe to compute performances --- ###
    df_results_whole = pd.DataFrame(data = Results_Matrix, columns = name)
    harvest = np.zeros([n_year, 14])
    
    for yearIndex in range(0, n_year):
        df_results = df_results_whole[(df_results_whole["Time"] > year_time_limit[yearIndex]) 
                                      & (df_results_whole["Time"] <= year_time_limit[yearIndex+1])].copy(deep=True)
        harvest[yearIndex, 0] = yearIndex + first_year
        
        # Computing some results
        proc_1year = np.sum(df_results["prod_MC"])
        harvest[yearIndex, 1] = proc_1year
        time_below_0 = np.sum((df_results["T_bio"] < 0)) * dt
        harvest[yearIndex, 2] = time_below_0

        # Extracting daytime
        df_productible = df_results[(df_results["at"] > 0)]
        Total_productible = len(df_productible) * dt
        harvest[yearIndex, 3] = Total_productible
        
        # Pigment productions
        k = 0
        for pigName in ["Chl_a", "Chl_b", "Lutein"]:
            proc_1year = np.sum(df_results["prod_MC"] * df_results[pigName])
            harvest[yearIndex, 4+k] = proc_1year
            k = k + 1
        
        # Absorbed energy
        harvest[yearIndex, 7] = np.sum(df_results["P_Abs"]) * 3600 * dt
        harvest[yearIndex, 8] = np.average(df_results[(df_results["at"] > 0)]["eta_PS"])
        harvest[yearIndex, 9] = np.average(df_results[(df_results["at"] > 0)]["regulated_Fraction"])
        harvest[yearIndex, 10] = np.average(df_results[(df_results["at"] > 0)]["I_average"])
        
        # Pigment content
        harvest[yearIndex, 11] = np.average(df_results["Chl_a"])
        harvest[yearIndex, 12] = np.average(df_results["Chl_b"])
        harvest[yearIndex, 13] = np.average(df_results["Lutein"])
        
        
        if not Sobol:
            np.save("results_vector_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth), harvest)
            np.save("results_vector_extended_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth), Results_Matrix)
    
    return harvest    

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

z vector
0 - fraction of the reservoir volume replaced when the transmitted light reaches the limit value: [0, 1[
1 - maintenance mode. 0 - null, 1 - constant light respiration, 2 - constant dark respiration, 3 - time varying
2 - temperature mode. 0 - no temperature dependency, 1 - biological processes are modulated by temperature
3 - retroaction of the pigment content on cross section. 0 - variation but no retroaction on cross section, 1 - variation and retroaction
4 - cross sections wavelength band dependency. 0 - no, 1 - yes

resultsArray
0 - current year, -
1 - mass of microalgae produced, kg
2 - duration below 0 °C, h
3 - total potential production duration (at > 0), h
4 - chlorophyll a produced, g
5 - chlorophyll b produced, g
6 - lutein produced, g
7 - absorbed energy, GJ
8 - average photoconversion efficiency (during daytime), % absolute
9 - average deviation to maximum photoconversion efficiency (during daytime), % maximum producticle
10- average illumination (during daytime), µmolPhoton/m²/s
11- average cell chlorophyll a content, mg/g
12- average cell chlorophyll b content, mg/g
13- average cell lutein content, mg/g

'''
resultsArray = biofacade(x=[0, 0.6, 0.9, 0.9], 
              y=[2, 1, 0.75, 0.08, 1, 1, 0, 0],
              z=[0.20, 3, 1, 1, 1],
              supDataMeteo = None, 
              supDataLight = None, 
              runID = 0, 
              stationID = -1)

print("")

HHV = 20.51e6 # J/kg -- Chlorella vulgaris HHV Oliver et al. 
eta_PS_0 = 0.0434 # -, maximum photosynthetic efficiency

for year in range(0, resultsArray.shape[0]):
    current_year = resultsArray[year, :]
    print("--- Year {:d} ---".format(int(current_year[0])))
    print("Produced biomass: {:.3f} kg".format(current_year[1]))
    print("Produced chl. a:  {:.3f} g".format(current_year[4]))
    print("Produced chl. b:  {:.3f} g".format(current_year[5]))
    print("Produced lutein:  {:.3f} g".format(current_year[6]))
    print("Absorbed energy:  {:.3f} GJ".format(current_year[7]/1e9))
    print("    Biomass eq.:  {:.3f} kg".format(current_year[7]/HHV*eta_PS_0))
    print("Avg. PCE:         {:.3f} %".format(current_year[8]*100))
    print("Avg. deviation:   {:.3f} %".format(current_year[9]*100))
    print("Avg. I avg.:      {:.3f} µmolPhoton/m²/s".format(current_year[10]))
    print("Avg. cell chl. a: {:.3f} g".format(current_year[11]))
    print("Avg. cell chl. b: {:.3f} g".format(current_year[12]))
    print("Avg. cell lutein: {:.3f} g".format(current_year[13]))
    print("Time below 0 °C:  {:.1f} h".format(current_year[2]))
    print("Total productible time: {:.1f} h".format(current_year[3]))
    print("")

print('Normal end of execution.\nExecution time: ' + format(time.time() - startTime, '.2f') + ' s')