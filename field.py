"""
Created on Tues Jun 6 

@author: judyjinn

Dog Tracking GPS Data



"""




'''                         #######################
#------------------------   ## ---   Set up  --- ##     ------------------------
                            #######################
'''
"""import every package you'll need"""
import os 
import scipy.io as scio #scipy
import numpy as np #to use arrays

# TkAGG requiredfor tkinter to properly import and allow simulataneous 
# matplotlib figures
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import matplotlib.animation as animation #animate plots, unused
import matplotlib.cm as cm
from matplotlib import gridspec
import h5py 
from datetime import datetime as dt
import sympy as sp  #solves equations
import pandas as pd 
from shapely.geometry import Point
from shapely.geometry import Polygon
import geopandas    
import random
import math
import tkinter as Tk #opening files
from tkinter import filedialog
import time
import re
import pickle
import sys
import scipy.stats as stats
import statsmodels.api as sm
from scipy import spatial
from statsmodels.formula.api import ols
import statsmodels.formula.api as smf
import gpxpy as gpx
from datetime import datetime as dt


os.getcwd() 
# auto-detect width of terminal for displaying dataframes
pd.set_option('display.max_columns',0) 
plt.ion() 
# turn off the chaining warning from pandas
pd.options.mode.chained_assignment = None


'''                         #######################
#------------------------   ## ---   Notes   --- ##     ------------------------
                            #######################
'''
'''

Head up is 1, head down is 0

Convert the GPS tracker CSV file using http://www.gpsvisualizer.com/ 
to a gpx file to put it into Google Maps

TODO: fix calculation for area between trails

'''


'''                          #######################
#-------------------------   ## --- Functions --- ##     -----------------------
                             #######################
'''

def convert_GPS(data, header, Hz):
    ''' Converts GPS from NMEA to lat/long format

    Args:
        data:       pandas data frame; contains raw GPS data
        header:     list of str; header to be put as colum names
        Hz:         recording rate of GPS

    Returns: 
        data:       pandas data frame; containing lat/long GPS coordinates
    '''
    
    data.columns = header
    data = data.drop('rm1',1)
    data = data.drop('rm2',1)
    
    # Convert to dec degrees
    # If using this code in the future, make sure lat/long NMEA continues to have 2 digits and 3 digits respectively for min/sec conversion
    def lat_NMEA2decdeg(data):
        lat_dec = (data['lat_NMEA']//100.0)+(data['lat_NMEA']%100.0/60.0)
        if (data['lat_dir'] =='S'):
            lat_dec *= (-1)
        return lat_dec
    
    def long_NMEA2decdeg(data):
        long_dec = (data['long_NMEA']//100.0)+(data['long_NMEA']%100.0/60.0)
        if (data['long_dir'] =='W'):
            long_dec *= (-1)
        return long_dec
        
    data['lat_dec'] = data.apply(lat_NMEA2decdeg, axis=1)
    data['long_dec'] = data.apply(long_NMEA2decdeg, axis=1)

    
    return data
    
def weather_station():
    ''' Concatenates all the mobile weather station data from the two experiments
        returns them as a pandas data frame

    Args:
        None

    Returns: 
        WS_data:   pandas data frame; containing lat/long GPS coordinates
                    from trail aging experiments as well as time of day expts
    '''
    
    TofD_WS_file = path+'/TofD_WS.csv'
    TofD_WS = pd.read_csv(TofD_WS_file, sep=',')
    trail_age_file = path+'/trail_age_WS.csv'
    trail_age_WS = pd.read_csv(trail_age_file, sep=',')
    WS_data = pd.concat([TofD_WS, trail_age_WS])
    WS_data = WS_wind(WS_data)
    return WS_data
    
    
def clean_trail(raw_GPS, trail_name, header, Hz, trail_path):
    ''' Transforming data for anlyses. Adds a time column, calculates centroids
        of GPS data to downsample from 5 Hz to 1 Hz.

    Args:
        raw_GPS:        pandas data frame; GPS data converted to lat/long
        trail_name:     str; file name for GPS data
        header:         list of str; contains column names
        Hz:             int; speed of GPS recording
        trail_path:     str; path to file

    Returns: 
        raw_GPS:        pandas data frame; GPS with lat/long and added time col
        centroid_1Hz:   pandas data frame; downsampled GPS data
    '''
    
    # Modify Time of Day trail by finding average GPS points every second
    raw_GPS = convert_GPS(raw_GPS, header, Hz)
    
    raw_GPS['time'] = raw_GPS['hr'].apply(lambda x: str(x)) + ':' + \
        raw_GPS['min'].apply(lambda x: str(x))+ ':' + \
        raw_GPS['sec'].apply(lambda x: str(x).zfill(2))

    if (os.path.isfile(
        trail_path + '/' + trail_name + '_centroid_1Hz.csv')==False
        ):
        
        print('Calculating centroids')
        # find num of observations per time (should be about 5 because of 5Hz)
        # Use list(t_grp) to view
        t_grp =  raw_GPS[['lat_dec', 'long_dec']].groupby( 
            raw_GPS['time'], sort=False
            )
        t_grp_size = raw_GPS[['lat_dec', 'long_dec']].groupby( 
            raw_GPS['time'], sort=False).size().reset_index()
        t_grp_size.columns=['time', 'XYgrp_size']
        t_grp_lat = raw_GPS['lat_dec'].groupby( 
            raw_GPS['time'], sort=False).sum()
        t_grp_long = raw_GPS['long_dec'].groupby( 
            raw_GPS['time'], sort=False
            ).sum()
        centroid_latlong = pd.concat([t_grp_lat,t_grp_long], 
            axis=1).reset_index()
        centroid_latlong = pd.merge(centroid_latlong,t_grp_size, on='time')
        centroid_latlong['center_lat'] = (
            centroid_latlong['lat_dec']/centroid_latlong['XYgrp_size']
            )
        centroid_latlong['center_long'] = (
            centroid_latlong['long_dec']/centroid_latlong['XYgrp_size']
            )

        centroid_1Hz = centroid_latlong[['center_lat','center_long']]
        centroid_1Hz.columns = ['lat','long']
        centroid_1Hz.to_csv(trail_path + '/' + trail_name + '_centroid_1Hz.csv', 
            sep=',',index=False
            )
    else:
        centroid_1Hz = pd.read_csv(trail_path + '/' + trail_name + \
        '_centroid_1Hz.csv'
        )
            
    return raw_GPS, centroid_1Hz
    
def graph_trail(dog_name, dog, dog_tilt, dog_WS, trail, condition, save_path):
    ''' Graphs human trail, dog trail, dog head up/down, and wind directions
        Saves graph as png

    Args:
        dog_name:   str; name of dog
        dog:        pd df; GPS trail of dog
        dog_tilt:   pd df; head up/down of dog
        dog_WS:     pd df; mobile weather station readings for dog's trail
        trail:      pd df; GPS of human trail
        condition:  str; name of experiment
        save_path:  str; path to save location

    Returns: 
        None
    '''
    
    fig = plt.figure(figsize = (15,5))
    ax = fig.gca()
    plt.style.use('ggplot')
    ax.grid(color='lightgray', linestyle='-', linewidth=1)
    plt.axis('off')
    # plt.ticklabel_format(useOffset=False) # turn off scientific plotting
    fig.suptitle(dog_name+ ' ' + condition, fontsize = 20)

    ax.plot(trail['long'], trail['lat'], color = 'crimson', linewidth=3, 
        label='Trail'
        )  
    ax.plot(dog['long'],dog['lat'], color = 'limegreen', linewidth=2, 
        label=dog_name
        )        
    ax.scatter(dog_tilt['long_dec'], dog_tilt['lat_dec'], color = 'darkgreen', 
        s = 40, zorder = 6, alpha = 0.5, label = (dog_name + 'Head down')
        )    
    ax.scatter(dog_WS['wp_long'], dog_WS['wp_lat'], color = 'crimson', s = 40, 
        marker = "s", zorder = 5, label = 'Weather Station Points'
        )    
    ax.quiver(dog_WS['wp_long'], dog_WS['wp_lat'], 
        dog_WS['wind_cos(rad)']*dog_WS['Wind Speed'], # scale by wind speed
        dog_WS['wind_sin(rad)']*dog_WS['Wind Speed'], 
        zorder = 1, width = 0.002, color = 'crimson',alpha=0.5,
        label = 'Wind Direction'
        )
        
    ax.set_facecolor('white')
    ax.axis('off')
    ax.legend(loc='lower left').draggable()
    plt.subplots_adjust(bottom = None, top = 0.9)
    fig.savefig((save_path+'/'+dog_name+'_'+condition+'.png'))
    plt.close()
    
    return
    
def WS_wind(WS_file):
    ''' Finding vert/horiziontal components of wind direction for graphing

    Args:
        WS_file:    pd df; all mobile weather station data

    Returns: 
        WS_file:    pd df; containing two new cols for wind direction components
    '''
    
    def wind_cos(data):
        data['wind_cos(rad)'] = math.cos(data['Wind Direction (deg)'] * 
            math.pi/180.0)
        return data
    def wind_sin(data):
        data['wind_sin(rad)'] = math.sin(data['Wind Direction (deg)'] * 
            math.pi/180.0)
        return data
    WS_file['Wind Direction (deg)']=WS_file['Wind Direction (deg)']-90
    WS_file = WS_file.apply(wind_cos, axis=1)
    WS_file = WS_file.apply(wind_sin, axis=1)

    return WS_file

def find_dog_file(files_path, dog_name, condition):
    ''' Pulls correct CSV file for dog trail being analyzed

    Args:
        files_path: str; path to file
        dog_name:   str; name of dog
        condition:  str; name of experiment

    Returns: 
        dog_track:  pd df; data of dog GPS and head up/down
    '''
    
    tracks_files = os.listdir(files_path)
    for file_name in tracks_files:
        if ((dog_name in file_name) & (condition.replace(' ','') in file_name)):
            track_file =  files_path + '/' + file_name
            dog_track = pd.read_csv(track_file, sep=',', header=None)
    return dog_track


def area_btwn_trails(dog_centroid, trail):
    ''' Creates a polygon from the human trail and dog trail.
        Calculates the area of polygon.
        Note: Currently there is some problem with this calculation. Fix later.

    Args:
        dog_centroid:   pd df; downsampled GPS data for dog
        trail:          pd df; GPS data of human trail

    Returns: 
        area:           float; Area of the polygon created
    '''
    
    dog_centroid_list = dog_centroid.values.tolist()
    trail_list=trail.values.tolist()
    polygon_points = []
    
    # Concatenate all points together and close off the beginnig and end
    for xy in dog_centroid_list:
        polygon_points.append(xy)
    
    for xy in trail_list[::-1]:
        polygon_points.append(xy)

    for xy in dog_centroid_list[0:1]:
        polygon_points.append(xy)

    # create a polygon from all points
    polygon = Polygon(polygon_points)
    area = polygon.area
    
    return area

def distance_from_trail(dog_centroid, trail, file_name, dist_path,):
    ''' Uses K-D trees using all points from the dog's trail to find the nearest
        GPS point to the human's trail. Saves as CSV or opens if already exists

    Args:
        dog_centroid:   pd df; downsampled GPS data for dog
        trail:          pd df; GPS data of human trail
        file_name:      str; what to name the saved file
        dist_path:      str; path to save/open the distance CSV file

    Returns: 
        dog2trail_dist: pd df; data containing the distance between points
    '''
    
    if not os.path.exists(dist_path):
        os.makedirs(dist_path)
    
    if (os.path.isfile(
        dist_path + '/' + file_name + '_dist2path.csv')==False
        ):
    
        print('calculating distances... ', file_name)
        trail_X = trail.as_matrix()
        dog_centroid_X =  dog_centroid.as_matrix()
    
        KDTree_nearest = []
        for pt in dog_centroid_X:
            dist, index = spatial.KDTree(trail_X).query(pt)
            lat_x, long_y = trail_X[spatial.KDTree(trail_X).query(pt)[1]]
            KDTree_nearest.append([dist, index, lat_x, long_y])

        tmp = pd.DataFrame(KDTree_nearest)
        dog2trail_dist = pd.concat([dog_centroid, tmp],axis=1)
        dog2trail_dist.columns = ['dog_lat', 'dog_long', 'dist2trail', 
            'trail index', 'trail_lat', 'trail_long'
            ]
        dog2trail_dist.to_csv(dist_path + '/' + file_name + \
            '_dist2path.csv', sep=',',index=False
            )
    else:
        dog2trail_dist = pd.read_csv(dist_path + '/' + file_name + \
        '_dist2path.csv'
        )
    
    return dog2trail_dist

def time_of_day(path, header, Hz, field_info, WS_data):
    ''' Calculates stats about the time of day experiment regarding trail length,
        time to completing the trail, meteorlogical data, etc.

    Args:
        path:       str; main path of script
        header:     list of str; names for coluns
        Hz:         int; rate of GPS recordings
        field_info: pd df; contains information about each experiment
        WS_data:    pd df; all mobile weather station readings

    Returns: 
        ToD_trail_diff: pd df; data containing stats about each dog's trail
    '''
    
    # Time of Day trail
    # Find centroids for raw GPS data if file has not yet been created
    trail_path = (path + '/TrackingData/briones_TofD/'
        'Trail')
    if (os.path.isfile(
        trail_path + '/TofD_trail_centroid_1Hz.csv')==False
        ):
        raw_TofD_trail = pd.read_csv(trail_path+'/2017_06_07_BrionesTrail.TXT',
            sep = ','
            )
        raw_TofD_trail = raw_TofD_trail.iloc[:-1500]
        TofD_trail = clean_trail(raw_TofD_trail,'TofD_trail', header, Hz, 
            trail_path
            )

    else:
        TofD_trail = pd.read_csv(path + '/TrackingData/briones_TofD/' 
            'Trail/Briones_1hr_centroid1Hz.csv', sep=','
            )
    # check for existing graph folder    
    TofD_graph_folder = path+'/TrackingData/briones_TofD/graphs'    
    if not os.path.exists(TofD_graph_folder):
        os.makedirs(TofD_graph_folder)
        
    field_info_TofD = field_info[
        field_info['expt'] == 'Time of Day'].reset_index() 
    list_conditions = field_info_TofD['condition'].unique().tolist()
    
    # create empty dict to store information about each dog
    TofD_stats_dict = {}
    for condition in list_conditions:
        TofD_stats_dict[condition]={}
    
    # Loop through all dog trials and calculate data for each
    TofD_files = path+'/TrackingData/briones_TofD/dogs'
    TofD_dist_folder = path+'/TrackingData/briones_TofD/dogs/dog2trail_dist'
    for trial in range(0, len(field_info_TofD)):
        date = (field_info_TofD.iloc[trial]['date'])
        dog_name = (field_info_TofD.iloc[trial]['dog'])
        trail = (field_info_TofD.iloc[trial]['trail_laid_by'])
        condition = (field_info_TofD.iloc[trial]['condition'])
        
        # Find correct file
        dog_track = find_dog_file(TofD_files, dog_name, condition)

        # Check if centroids folder exists
        TofD_centroids_folder = path+'/TrackingData/briones_TofD/dogs/centroids'
        if not os.path.exists(TofD_centroids_folder):
            os.makedirs(TofD_centroids_folder)
        
        # Downsample GPS data if not already done
        save_name = (dog_name + '_' + condition)
        dog_track, dog_centroid = clean_trail(dog_track, save_name, header, Hz, 
            TofD_centroids_folder)
        
        # Find start/stop time for each trial and then cut GPS data accordingly
        try:
            start_row = dog_track[dog_track['time']==start_time].index[0]
            dog_track=dog_track[start_row:]
        except:
            pass
            
        try:
            end_row = dog_track[dog_track['time']==end_time].index[0]
            dog_track=dog_track[:end_row]
        except:
            pass
        
        # Subset dog head up/down data
        dog_tilt = dog_track[dog_track['tilt']==0]
        # Subset mobile weather station data for the dog and trail
        dog_WS = WS_data[(WS_data['Date']==date) & (WS_data['Dog']==dog_name)]
        
        # # Graph all TofD trials
        # # Comment out if unneeded
        # graph_trail(dog_name, dog_centroid, dog_tilt, dog_WS, TofD_trail,
        #     condition, TofD_graph_folder)
        
        # Find area between trail and dog trail
        traildog_area = area_btwn_trails(dog_centroid, TofD_trail)
        
        # Find distance from trail and dog trail
        file_name = dog_name + '_' + condition.replace(' ','')
        dog_dist2trail = distance_from_trail(dog_centroid, TofD_trail, 
            file_name, TofD_dist_folder
            )
        
        # Stats
        start = dog_track['time'].iloc[0]
        end = dog_track['time'].iloc[-1]
        total_sec, total_time = trial_time(start, end)
        
        
        TofD_stats_dict[condition][dog_name] = {
            'total seconds': total_sec,
            'total time': total_time,
            'total head down': len(dog_tilt),
            'percent head down': len(dog_tilt)/len(dog_track),
            'trail length': field_info_TofD.iloc[trial]['trail length'],
            'dog trail length': field_info_TofD.iloc[trial]['dog trail length'],
            'trial number': field_info_TofD.iloc[trial]['trial num'],
            'start order': field_info_TofD.iloc[trial]['start order'],
            'btwn trail area': traildog_area
            }
    
    # Create a pandas data frame out of the dict of all dog trails    
    TofD_stats_df = pd.DataFrame.from_dict(
        {
            (i,j): TofD_stats_dict[i][j] 
            for i in TofD_stats_dict.keys() 
            for j in TofD_stats_dict[i].keys()
        }, orient='index').reset_index()
        
    TofD_stats_df = TofD_stats_df.rename(index=str, 
        columns={"level_0": "Condition", "level_1": "Dog"}
        )
     
    # Get averages for meteorological data   
    trail_avg_diff = weather_station_stats(WS_data)
    
    ToD_trail_diff = trail_avg_diff[
        trail_avg_diff['Experiment']=='Time of Day'
        ]
    ToD_trail_diff = pd.merge(
        TofD_stats_df, ToD_trail_diff, on=['Dog', 'Condition']
        )
    
    ToD_trail_diff['dog2trail_ratio'] = ToD_trail_diff['dog trail length']/ \
        ToD_trail_diff['trail length']
    
    ToD_trail_diff.to_csv(path+'/TrackingData/briones_TofD/ToD_trail_diff.csv', 
        sep=',', index=False
        )
    
    return ToD_trail_diff 
    
def trail_age(path, header, Hz, field_info, WS_data):
    ''' Calculates stats about the time of day experiment regarding trail length,
        time to completing the trail, meteorlogical data, etc.

    Args:
        path:       str; main path of script
        header:     list of str; names for coluns
        Hz:         int; rate of GPS recordings
        field_info: pd df; contains information about each experiment
        WS_data:    pd df; all mobile weather station readings

    Returns: 
        trail_age_vals: pd df; data containing stats about each dog's trail
    '''
    
    # Time of Day trail
    # Find centroids for raw GPS data if file has not yet been created
    trail_path = (path + '/TrackingData/russell_trailage')
    trail_age_trail = pd.read_csv(trail_path + \
        '/Russell_trailage_trail.csv', sep=','
        )
    
    # check for existing graph folder      
    trail_age_graph_folder = path+'/TrackingData/russell_trailage/graphs'    
    if not os.path.exists(trail_age_graph_folder):
        os.makedirs(trail_age_graph_folder)
        
    field_info_trail_age = field_info[
        field_info['expt'] == 'Trail Age'].reset_index() 
    list_conditions = field_info_trail_age['condition'].unique().tolist()
    
    # create empty dict to store information about each dog
    trail_age_stats_dict = {}
    for condition in list_conditions:
        trail_age_stats_dict[condition]={}
    
    # Loop through all dog trials and calculate data for each
    trail_age_files = path+'/TrackingData/russell_trailage/dogs'
    trail_age_dist_folder = path+'/TrackingData/russell_trailage/dogs/dog2trail_dist'
    for trial in range(0, len(field_info_trail_age)):
        date = (field_info_trail_age.iloc[trial]['date'])
        dog_name = (field_info_trail_age.iloc[trial]['dog'])
        trail = (field_info_trail_age.iloc[trial]['trail_laid_by'])
        condition = (field_info_trail_age.iloc[trial]['condition'])
        
        # Find correct file
        dog_track = find_dog_file(trail_age_files, dog_name, condition)
        
        # Search for existing centroids files
        trail_age_centroids_folder = path+'/TrackingData/russell_trailage/dogs/centroids'
        if not os.path.exists(trail_age_centroids_folder):
            os.makedirs(trail_age_centroids_folder)
        
        # Create new centroids or load them if not preexisting
        save_name = (dog_name + '_' + condition)
        dog_track, dog_centroid = clean_trail(dog_track, save_name, header, Hz, 
            trail_age_centroids_folder)
        
        # Get start/end times for each trial and clip GPS & head data    
        try:
            start_row = dog_track[dog_track['time']==start_time].index[0]
            dog_track=dog_track[start_row:]
        except:
            pass
            
        try:
            end_row = dog_track[dog_track['time']==end_time].index[0]
            dog_track=dog_track[:end_row]
        except:
            pass
        
        # subset dog head up/down data
        dog_tilt = dog_track[dog_track['tilt']==0]
        dog_WS = WS_data[(WS_data['Date']==date) & (WS_data['Dog']==dog_name)]
        
        # # Graph all trail_age trials
        # graph_trail(dog_name, dog_centroid, dog_tilt, dog_WS, trail_age_trail,
        #     condition, trail_age_graph_folder)
        
        # Find area between trail and dog trail
        traildog_area = area_btwn_trails(dog_centroid, trail_age_trail)
        
        # Find distance from trail and dog trail
        file_name = dog_name + '_' + condition
        dog_dist2trail = distance_from_trail(dog_centroid, trail_age_trail, 
            file_name, trail_age_dist_folder
            )
        
        # Stats
        start = dog_track['time'].iloc[0]
        end = dog_track['time'].iloc[-1]
        total_sec, total_time = trial_time(start, end)
        
        
        trail_age_stats_dict[condition][dog_name] = {
            'total seconds': total_sec,
            'total time': total_time,
            'total head down': len(dog_tilt),
            'percent head down': len(dog_tilt)/len(dog_track),
            'trail length': field_info_trail_age.iloc[trial]['trail length'],
            'dog trail length': field_info_trail_age.iloc[trial]['dog trail length'],
            'trial number': field_info_trail_age.iloc[trial]['trial num'],
            'start order': field_info_trail_age.iloc[trial]['start order'],
            'btwn trail area': traildog_area
            }
        
    trail_age_stats_df = pd.DataFrame.from_dict(
        {
            (i,j): trail_age_stats_dict[i][j] 
            for i in trail_age_stats_dict.keys() 
            for j in trail_age_stats_dict[i].keys()
        }, orient='index').reset_index()
        
    trail_age_stats_df = trail_age_stats_df.rename(index=str, 
        columns={"level_0": "Condition", "level_1": "Dog"}
        )
        
    trail_age_vals = weather_station_stats(WS_data)
    
    trail_age_vals = trail_age_vals[
        trail_age_vals['Experiment']=='Trail Age'
        ]
    trail_age_vals = pd.merge(
        trail_age_stats_df, trail_age_vals, on=['Dog', 'Condition']
        )
    
    trail_age_vals['dog2trail_ratio'] = trail_age_vals['dog trail length']/ \
        trail_age_vals['trail length']
    
    trail_age_vals.to_csv(path+'/TrackingData/russell_trailage/trail_age_vals.csv', 
        sep=',', index=False
        )
    
    return trail_age_vals
    
def trial_time(start, end):
    '''
    # Function: get_msec
    # Input: a series of data where values are str in hh:mm:ss.ms format 
    # Outputs: np.array of time converted to milliseconds (int)
    '''

    h,m,s = list(map(int,end.split(':')))
    end_sec = (h*60*60)+(m*60)+s
    
    h,m,s = list(map(int,start.split(':')))
    start_sec = (h*60*60)+(m*60)+s
    
    total_sec = end_sec -start_sec
    
    h = (total_sec // 3600)        
    m = (total_sec % 3600 // 60)   
    s = (total_sec % 3600% 60)
    total_time = [h,m,s]  
    
    
    return total_sec, total_time

def weather_station_stats(WS_data):
    ''' Calculates the average meteorological readings either by day or by
        each dog's trail. Then calculates the change between when the person
        laid the trail and when the dog was trailing.

    Args:
        WS_data:   pd df; all mobile weather station data

    Returns: 
        trail_avg_diff: pd df; the average of meteorogical data
    '''
    
    dogs_WS_data = WS_data[
        (WS_data['Dog']=='Kapo') | (WS_data['Dog']=='Voz') | \
        (WS_data['Dog']=='Zinka') | (WS_data['Dog']=='Gig')
        ]
    
    # Find average of the entire day
    dogs_daily_avg = dogs_WS_data.groupby(
        ['Date', 'Condition']).mean().reset_index() 
    
    # Find average during one dog's run
    dogs_trail_avg = dogs_WS_data.groupby(
        ['Date', 'Experiment', 'Condition', 'Dog']).mean().reset_index()
    
    # Subset out people's trails
    WS_data_people = WS_data[
        (WS_data['Dog']!='Kapo') & (WS_data['Dog']!='Voz') & \
        (WS_data['Dog']!='Zinka') & (WS_data['Dog']!='Gig')
        ]
    
    #  Averages during the person's trail
    people_trail_avg = WS_data_people.groupby(
        ['Date', 'Condition', 'Dog']).mean().reset_index()
    
    # Merge dog trail avg with the people trails and find difference
    trail_avg = pd.merge(dogs_trail_avg, people_trail_avg, on = 'Date')
    col_names = people_trail_avg.columns.tolist()
    
    temp_diff = trail_avg['Temp_y'] - trail_avg['Temp_x']
    trail_avg_diff = pd.DataFrame(temp_diff)
    trail_avg_diff.columns = ['temp_change']
    
    #  Find the difference and make new columns
    for col in col_names[4:8]:
        trail_avg_diff[col+'_change'] = (
        trail_avg[col+'_y'] - trail_avg[col+'_x']
        )
    trail_avg_diff = pd.concat([dogs_trail_avg, trail_avg_diff],axis=1
        )
        
    return trail_avg_diff
    
'''                       #######################
#----------------------   ## ---    MAIN   --- ##     --------------------------
                          #######################
'''

if __name__ == '__main__':
    path = os.getcwd()
    Hz = 5
    header = ['rm1', 'year', 'month', 'day', 'hr', 'min', 'sec', 'ms', 
        'lat_NMEA', 'lat_dir', 'long_NMEA', 'long_dir', 'rm2', 'tilt']

    field_info = pd.read_csv(path+'/field_info.csv', sep=',')
    WS_data = weather_station()
    # ToD_trail_diff = time_of_day(path, header, Hz, field_info, WS_data)
    trail_age_vals = trail_age(path, header, Hz, field_info, WS_data)
    
     
     
     
     
    # Testing area below this line
    # path='/Users/judyjinn/Python/DogTracking/Field/TrackingData/briones_TofD/dogs/dog2trail_dist'
    # dog_dist2trail = pd.read_csv(path+'/Gig_afternoon_dist2path.csv', sep=',')
    # dist_freq = dog_dist2trail.groupby(['trail index', 'trail_lat','trail_long']).count()['dist2trail'].reset_index()
    # dist_freq.columns = ['trail index', 'trail_lat', 'trail_long', 'freq']
    
    # n = 50
    # x = np.random.randn(n)
    # y = x * np.random.randn(n)
    #
    # fig, ax = plt.subplots(2, figsize=(6, 6))
    #
    # ax[0].scatter(x, y, s=50)
    #
    # sizes = (np.random.randn(n) * 8) ** 2
    # ax[1].scatter(x, y, s=sizes)
    #
    # fig.show()
    
        #
    # fig = plt.figure(figsize = (15,5))
    # ax = fig.gca()
    # plt.style.use('ggplot')
    # ax.grid(color='lightgray', linestyle='-', linewidth=1)
    # plt.axis('off')
    # # plt.ticklabel_format(useOffset=False) # turn off scientific plotting
    #
    # plt.scatter(x,y)
    # ax.set_facecolor('white')
    # ax.axis('off')
    # ax.legend(loc='lower right').draggable()
    # plt.subplots_adjust(bottom = None, top = 0.9)
    # plt.show()
    #
