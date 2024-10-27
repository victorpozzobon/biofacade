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
import pandas as pd
import matplotlib.pyplot as plt
import time
import helpFunctions.Solar_Time as st

def extracting_and_subsampling(station_name = "Marignane", dt = 0.1, azimuth = 0, transient_day = 7):
    # azimuth = 0 # azimuth par rapport au sud de la surface verticale, plein sud = 0, est = -90, ouest + 90
    startTime = time.time()
    # Pandas option
    pd.set_option('display.max_columns', None)
    
    # station_name = "Marignane"
    # dt = 0.5
    # azimuth = 0
    # transient_day = 0
    
    # Loading the data
    print('Starting to load the database ...')
    currentTime = time.time()
    dataset = pd.read_csv("donnees-synop-essentielles-omm.csv", delimiter = ';')
    print('Database loaded in {:3.3g} s'.format(time.time() - currentTime))
    
    # To print columns names:
    # dataset.columns
    
    # Restricting dataset
    dataset_short = dataset.loc[: , ("ID OMM station", "Date", 'Direction du vent moyen 10 mn', 'Vitesse du vent moyen 10 mn', 'Température', 'Nebulosité totale', 'Latitude',
    'Longitude', 'communes (name)', 'Point de rosée', 'Pression station', 'Humidité')]
    
    # Replacing Nan nebulosity by 0
    print('Number of Nan before correcting nebulosity: {:d}'.format(dataset_short.isna().sum().sum()))
    dataset_short['Nebulosité totale'].fillna(0, inplace=True)
    print('Number of Nan after correcting nebulosity: {:d}'.format(dataset_short.isna().sum().sum()))
    
    # Checking Nan left
    dataset_short_na = dataset_short[dataset_short.isna().any(axis=1)]
    plt.figure()
    plt.hist(dataset_short_na['ID OMM station'])
    plt.title("Number of Na per station")
    plt.title("ID OMM station")
    plt.show()
    
    for station in dataset_short_na['ID OMM station'].unique():
        data_station = dataset_short[dataset_short["ID OMM station"] == station]
        data_station_technique = data_station.loc[: , ("Date", 'Direction du vent moyen 10 mn', 'Vitesse du vent moyen 10 mn', 'Température', 'Nebulosité totale', 'Point de rosée', 'Pression station')]
        nb_nan = data_station_technique.isna().sum().sum()
        if nb_nan > 0:
            nb_tot = len(data_station)
            name = data_station['communes (name)'].unique()[-1]
            long = data_station['Longitude'].unique()[-1]
            lat = data_station['Latitude'].unique()[-1]
            print("Station {:s}, long {:f}, lat {:f}, nb nan {:d}, nb tot {:d}, {:3.3g} %".format(str(name), long, lat, nb_nan, nb_tot, nb_nan/nb_tot * 100))
            
        # To print the nan surrounding:
        # np.where(data_station_technique.isna())[0]
        # data_station_technique.iloc[36320:36330, :]
        
    # Isolating "Marignane" by default
    currentTime = time.time()
    print("\nStarting to process the test station ...")
    
    print("Station name: " + station_name)
    # station_name = "Marignane"
    if not station_name.isnumeric():
        data_meteo = dataset_short[dataset_short['communes (name)'] == station_name].dropna().reset_index()
    else:
        data_meteo = dataset_short[dataset_short['ID OMM station'] == int(station_name)].dropna().reset_index()
    # Orly
    # data_meteo = dataset_short[dataset_short['ID OMM station'] == 7149].dropna().reset_index()
    # Marignane
    # data_meteo = dataset_short[dataset_short['ID OMM station'] == 7650].dropna().reset_index()
    
    # Sorting the times and computing meteo
    # Preallocating
    data_meteo["year"] = -1
    data_meteo["julianDay"] = -1
    data_meteo["hour"] = -1
    data_meteo["timeStamp"] = -1
    data_meteo["meteoMode"] = -1
    data_meteo["cloudCover"] = -1
    data_meteo["UHII"] = -1
    
    # Dealing with each rwo
    for i in range(0, len(data_meteo)):
        # Extracting time
        extractedDate = data_meteo.loc[i, 'Date']
        year = int(extractedDate.split('-')[0])
        month = int(extractedDate.split('-')[1])
        day = int(extractedDate.split('-')[2].split('T')[0])
        hour = int(extractedDate.split('T')[1].split(':')[0]) - (int(extractedDate.split('+')[1].split(':')[0])==2)*1 # There is a "-1" in the Solar_Time librairy
    
        # Allocating
        data_meteo.loc[i, 'timeStamp'] = int('{:d}{:02d}{:02d}{:02d}'.format(year, month, day, hour))
        data_meteo.loc[i, 'year'] = year
        ts = pd.Timestamp(year = year, month = month, day = day, hour = 6, second = 0, tz = 'Europe/Berlin')
        ts_ref = pd.Timestamp(year = year, month = 1, day = 1, hour = 6, second = 0, tz = 'Europe/Berlin')
        data_meteo.loc[i, 'julianDay'] = int(ts.to_julian_date() - ts_ref.to_julian_date() + 1)
        data_meteo.loc[i, 'hour'] = hour
        
        # Meteo
        # Nébulosité totale = Sky cover method. For meteo mode : 0 to 3 tenths -> mode 0, 4 to 7 tenths -> mode 1, 8 to 10 tenths -> mode 2
        nebulosity = int(data_meteo.loc[i, 'Nebulosité totale'] / 10)
        meteoMode = - 1
        if nebulosity >=0 and nebulosity <= 3:
            meteoMode = 0
        elif nebulosity >=4 and nebulosity <= 7:
            meteoMode = 1
        elif nebulosity >=8 and nebulosity <= 10:
            meteoMode = 2
        else:
            print("Problem on nebulosity for index {:d}".format(i))
        data_meteo.loc[i, 'meteoMode'] = meteoMode
        data_meteo.loc[i, "cloudCover"] = data_meteo.loc[i, 'Nebulosité totale'] / 100.
    
    # Caculating UHII of the day
    # Selecting each full day
    for currentYear in data_meteo["year"].unique():
        for currentDay in data_meteo["julianDay"].unique():
            # Computing the average and UHII
            U_avg = np.average(data_meteo[(data_meteo["year"] == currentYear) & (data_meteo["julianDay"] == currentDay)]["Vitesse du vent moyen 10 mn"].values)
            CC_avg = np.average(data_meteo[(data_meteo["year"] == currentYear) & (data_meteo["julianDay"] == currentDay)]["cloudCover"].values)
            RH_avg = np.average(data_meteo[(data_meteo["year"] == currentYear) & (data_meteo["julianDay"] == currentDay)]["Humidité"].values)
            UHII = -0.54 * U_avg - 1.48 * CC_avg - 0.039 * RH_avg + 7.63
            data_meteo.loc[(data_meteo["year"] == currentYear) & (data_meteo["julianDay"] == currentDay), "UHII"] = UHII
            
    # Sorting by record time
    data_meteo.sort_values(by='timeStamp', ascending=True, inplace=True)
    data_meteo = data_meteo.reset_index()
    
    
    # %%% Theoretical values %%% #
    # Contructing the time vector
    longitude = data_meteo["Longitude"].iloc[0]
    latitude = data_meteo["Latitude"].iloc[0]
    # azimuth = 0 # azimuth par rapport au sud de la surface verticale, plein sud = 0, est = -90, ouest + 90
    meteo_mode = 0
    
    day_time_vector = []
    illuminations_vector = []
    for day in range(1, 365+1):
        for hours in np.arange(0, 24, 0.25):
            day_time_vector.append(np.array([day, hours, latitude, longitude, azimuth, meteo_mode]))
            illuminations_vector.append(st.solar_input(
                day, hours, latitude, longitude, azimuth, meteo_mode))
    day_time_vector = np.array(day_time_vector)
    illuminations_vector = np.array(illuminations_vector)
    
    # PHOTO-CONVERSION
    # 1 W/m2 = 126.7 lux
    # 1 lux = 0.0185 µmolPhoton/m2/s 
    # 1 W/m2 = 126.7*0.0185 = 2.3495 µmolPhoton/m2/s 
    # OR 1 lux = 0.0185 µmolPhotonPAR/m2/s  https://hortione.com/2021/08/lux-to-ppfd-convertor-ppfd-to-lux/
    # The above is more coherent
    # Rendement PAR
    # In case: eta = 0.423 # https://www.fondriest.com/environmental-measurements/parameters/weather/photosynthetically-active-radiation/
    
    def WattToPAR(x):
        return 2.3495 * x
    def PARtoWatt(x):
        return 1 / 2.3495 * x
    
    plt.figure(figsize=(6, 12))
    ax = plt.subplot(1, 2, 1)
    ax.plot(day_time_vector[:, 0] * 24 + day_time_vector[:, 1], illuminations_vector[:, 0] * 1000 /126.7, 
             label = "DNI", alpha = 0.5)
    ax.plot(day_time_vector[:, 0] * 24 + day_time_vector[:, 1], illuminations_vector[:, 2] * 1000 /126.7, 
             label = "Vertical surface", alpha = 0.5)
    ax.plot(day_time_vector[:, 0] * 24 + day_time_vector[:, 1], illuminations_vector[:, 1] * 1000 /126.7, 
             label = "Horizontal surface", alpha = 0.5)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("DNI (W/m²)")
    ax.set_ylim([0, 1000])
    secax_y = ax.secondary_yaxis(
        'right', functions=(WattToPAR, PARtoWatt))
    secax_y.set_ylabel('µmolPhotonPAR/m²/s')
    
    plt.legend()
    plt.show()
    
    
    # %%% Meteo modulated values %%% #
    # Year selection
    transient_day = 7
    selected_years = np.arange(2013, 2022+1, 1)
    data_meteo_1year = data_meteo[
        ((data_meteo['year'] == selected_years[0]-1) & (data_meteo['julianDay'] >= 365 - transient_day)) 
        | ((data_meteo['year'] >= selected_years[0]) & (data_meteo['year'] <= selected_years[-1])) 
        | ((data_meteo['year'] == selected_years[-1]+1) & (data_meteo['julianDay'] == 1))]
    
    # EMERGENCY MODE
    # 1 year mode
    # data_meteo_1year = data_meteo[(data_meteo['year'] == 2022)]
    
    # Just for testing
    # data_meteo_1year.loc[:, "meteoMode"] = 0
    
    # Subsampling
    t_old = data_meteo_1year.loc[:, 'julianDay'].values * 24 + data_meteo_1year.loc[:, 'hour'].values + data_meteo_1year['year'].values * 365 * 24
    t_old = t_old - np.min(t_old) +1 # Because the first recording of the day is at time 01:00 am
    t_new = np.arange(1, np.max(t_old)+dt, dt) # h
    meteoMode_new = np.round(np.interp(t_new, t_old, data_meteo_1year.loc[:, 'meteoMode'].values)).astype('int')
    dirWind_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'Direction du vent moyen 10 mn'].values)
    velWind_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'Vitesse du vent moyen 10 mn'].values)
    temperature_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'Température'].values)
    pressure_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'Pression station'].values)
    julianDay_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'julianDay'].values)
    hour_new = np.zeros(len(t_new))
    hour_new[0] = 1 # Because the first recording of the day is at time 01:00 am
    newDay = False
    for i in range(1, len(hour_new)):
        # Dealing with hours
        hour_new[i] = hour_new[i-1] + dt
        if hour_new[i] >= 24: 
            hour_new[i] = hour_new[i] - 24
            newDay = True
        # Dealing with Julian days
        if not newDay:
            julianDay_new[i] = julianDay_new[i-1]
        else:
            julianDay_new[i] = julianDay_new[int(i+12/dt)]
            newDay = False
            
    nebulosite_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'Nebulosité totale'].values)
    cloudCover_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'cloudCover'].values)
    UHII_new = np.interp(t_new, t_old, data_meteo_1year.loc[:, 'UHII'].values)
    
    day_time_vector = []
    illuminations_vector = []
    clear_Sky_Heat_Flux = []
    # for i in data_meteo_1year.index:
    #     day_time_vector.append(np.array([data_meteo_1year.loc[i, 'julianDay'], data_meteo_1year.loc[i, 'hour'], latitude, longitude, azimuth, data_meteo_1year.loc[i, 'meteoMode']]))
    #     illuminations_vector.append(st.solar_input(
    #         data_meteo_1year.loc[i, 'julianDay'], data_meteo_1year.loc[i, 'hour'], latitude, longitude, azimuth, data_meteo_1year.loc[i, 'meteoMode']))
    
    for i in range(0, len(t_new)):
        day_time_vector.append(np.array([julianDay_new[i], hour_new[i], 
                                         latitude, longitude, azimuth, meteoMode_new[i]]))
        illuminations_vector.append(st.solar_input(
            julianDay_new[i], hour_new[i], latitude, longitude, azimuth, meteoMode_new[i]))
        clear_Sky_Heat_Flux.append(st.solar_input(
            julianDay_new[i], hour_new[i], latitude, longitude, azimuth, 0) [1] * 1000 / 126.7) # DNI on horizontal surface, in W/m²
    
    
    day_time_vector = np.array(day_time_vector)
    illuminations_vector = np.array(illuminations_vector)
    clear_Sky_Heat_Flux = np.array(clear_Sky_Heat_Flux)
    
    ax = plt.subplot(1, 2, 2)
    ax.plot(day_time_vector[:, 0] * 24 + day_time_vector[:, 1], illuminations_vector[:, 0] * 1000 /126.7, 
             label = "DNI", alpha = 0.5)
    ax.plot(day_time_vector[:, 0] * 24 + day_time_vector[:, 1], illuminations_vector[:, 2] * 1000 /126.7, 
             label = "Vertical surface", alpha = 0.5)
    ax.plot(day_time_vector[:, 0] * 24 + day_time_vector[:, 1], illuminations_vector[:, 1] * 1000 /126.7, 
             label = "Horizontal surface", alpha = 0.5)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("DNI (W/m²)")
    ax.set_ylim([0, 1000])
    secax_y = ax.secondary_yaxis(
        'right', functions=(WattToPAR, PARtoWatt))
    secax_y.set_ylabel('µmolPhotonPAR/m²/s')
    plt.tight_layout(pad=0.5)
    plt.legend()
    plt.show()
    
    
    ### --- Storing data --- ###
    data_to_write = np.c_[t_new, temperature_new, cloudCover_new, clear_Sky_Heat_Flux, UHII_new, pressure_new, velWind_new, dirWind_new]
    
    np.savetxt("data/meteo_vector_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth) + ".csv", data_to_write)
    np.savetxt("data/illuminations_vector_" + station_name + "_dt_" + str(dt) + "_azi_" + str(azimuth) + ".csv", illuminations_vector)


print("Preprocessing done ...")

