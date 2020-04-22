"""

INFO
####

A script for home country detection for each individual Twitter user.

Representing two methods:
    
    1) Home is a country where user's geotagged tweet count is the greatest
    based on unique days classification
    
    2) Home is a country where user's geotagged tweet count is the greatest
    based on unique weeks classification


USAGE
#####

Run the script with the following command:
    python detect_homes.py -i path/to/tweets.pkl -c column_name -o path/to/output.pkl
    -uh path/to/user_homes.csv -uc path/to/user_count_per_area.csv


NOTE
####

Run post_location_detector.py before running this script!

The same script can be used to detect probable home continent, country or
municipality. Although with smaller areas, the reliability of home detection
diminishes due to idiosyncracies in each user's posting and geotagging
behavior and due to how location data is handled by Twitter.


@author: Samuli Massinen

edited by: Tuomas Väisänen
"""

# Import necessary libraries
import pickle
import pandas as pd
import geopandas as gpd
import sys
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define input file path
ap.add_argument("-i", "--input", required=True,
                help="Path to input pickle containing locations of posts.")

# Define column name for post location names
ap.add_argument("-c", "--column", required=True, default='post_country',
                help="Name of column containing location name of post. Default: "\
                    "post_country")

# Define output file path
ap.add_argument("-o", "--output", required=True,
                help="Path to the output pickle.")

# Define output file path for user home csv
ap.add_argument("-uh", "--userhome", required=True,
                help="Path to output csv for home detection per unique user")

# Define output file path for user home csv
ap.add_argument("-uc", "--ucount", required=True,
                help="Path to output csv for user count per detected home location")

# Parse arguments
args = vars(ap.parse_args())

# Read Twitter data
with open(args['input'], 'rb') as f:
    data = pickle.load(f)

# Add new columns for datetime: one with everything, other with only years, months and days
# Add columns for different home detection methods: unique days and unique weeks. Also, create a column for unique week calculation
data['datetime'] = None
data['datetime'] = pd.to_datetime(data['created_at'])
data['datetime_stripped'] = None
data['datetime_stripped'] = data['datetime'].dt.strftime("%Y-%m-%d")
data['home_unique_days'] = None
data['home_unique_weeks'] = None
data['unique_weeks'] = None

# Group the data by userid in order to define home country easily for each individual
grouped = data.groupby('user_id')

def calculateUniqueWeeks(grouped_gdf):
    
    """
    A Function for calculating a unique week for every individual Twitter post in GeoDataFrame (week number  + year).
    
    Parameter
    ----------
    grouped_gdf: <gpd.GeoDataFrame>
    
        A GeoDataFrame grouped by some unique user ID (e.g. 'userid'). Hence, key refers to users Twitter ID and values
        to every Twitter post an individual has made.
        
    Output
    ------
    <gpd.GeoDataFrame>
        A GeoDataFrame including unique week information for each individual post.
     
    """
    
    y = 1
    # Create an accessory list for contemporary storing of values
    user_list = []
    
    for key, values in grouped_gdf:
        
        print("Processing: ", y, "/" ,len(grouped_gdf))
        y = y + 1
        
        individual = values
        individual = individual.sort_values(by='datetime')
        
        for index, row in individual.iterrows():
            
            date = row['datetime']
            year = date.year
            week_number = date.isocalendar()[1]
            unique_week = str(week_number) + "_" + str(year)
            
            individual.loc[index, 'unique_weeks'] = unique_week
            
        user_list.append(individual)
        print("Unique week added!")
    
    # Convert individuals list to GeoDataFrame
    twitter_home_week = gpd.GeoDataFrame(pd.concat(user_list, ignore_index=True))
    
    return(twitter_home_week)

def detectHomeCountry(grouped_gdf, home_detection_method, colname):
    
    """
    A Function for detecting user's home country.
    
    Parameters
    ----------
    grouped_gdf: <gpd.GeoDataFrame>
    
        A GeoDataFrame grouped by some unique user ID (e.g. 'userid'). Hence, key refers to users Twitter ID and values
        to every Twitter post an individual has made.
        
    home_detection_method: <str>
    
        Input either 'home_unique_days' or 'home_unique_weeks'
        
    Output
    ------
    <gpd.GeoDataFrame>
        A GeoDataFrame including user's home country.
     
    """
    
    y = 1
    # Create accessory lists for contemporary storing of values
    user_list = []
    
    # Loop over grouped df to extract post location durations
    for key, values in grouped_gdf:
        
        print("Processing: ", y, "/" ,len(grouped_gdf))
        y = y + 1
        
        # Extract individuals tweet history and order it by datetime
        individual = values
        individual = individual.sort_values(by='datetime')
        
        # check home detection method
        if home_detection_method == 'home_unique_days':
        
            method = individual.drop_duplicates(subset='datetime_stripped')
            value_list = method[colname].value_counts()
            value_dict = dict(method[colname].value_counts())
            
        elif home_detection_method == 'home_unique_weeks':
            
            method = individual.drop_duplicates(subset='unique_weeks')
            value_list = method[colname].value_counts()
            value_dict = dict(method[colname].value_counts())
        
        else:
            
            print("Specify either 'home_unique_days' or 'home_unique_weeks' as home_detection method!")
            sys.exit()
        
        # Check if c
        if (len(value_list) == 1 or value_list[0] != value_list[1]):
                
            home = next(iter(value_dict))
            
            for index, row in individual.iterrows():
                individual.loc[index, home_detection_method] = home
                    
            user_list.append(individual)   
            print(home, 'added!')
            
        elif value_list[0] == value_list[1]:
    
            home = next(iter(value_dict)) + ' or ' +  list(value_dict.keys())[1]
            for index, row in individual.iterrows():
                individual.loc[index, home_detection_method] = home
                    
            user_list.append(individual)   
            print(home, 'added!')
            print("Problematic! Two are same")
                
        elif value_list[0] == value_list[1] == value_list[3]:
                
            home = next(iter(value_dict)) + ' or ' +  list(value_dict.keys())[1] + ' or ' +  list(value_dict.keys())[2]
            
            for index, row in individual.iterrows():
                individual.loc[index, home_detection_method] = home
                
            user_list.append(individual)
            print(home, 'added!')
            print("Problematic! Three are same")
                
        else:
            home = 'Too many options'
            
            for index, row in individual.iterrows():
                individual.loc[index, home_detection_method] = home
                
            user_list.append(individual)
            print(home, 'added!')
            print("Problematic! Too many options")
            
    # Convert individuals list to GeoDataFrame
    twitter_home = gpd.GeoDataFrame(pd.concat(user_list, ignore_index=True))
    
    return(twitter_home)

# determine home country based on unique days
twitter_home_unique_days = detectHomeCountry(grouped, 'home_unique_days',
                                             args['column'])

# calculate unique weeks accumulatively, with home_unique_days
print('[INFO] - Creating grouped weeks by user id..')
grouped_week = twitter_home_unique_days.groupby('user_id')
print('[INFO] - Calculating unique weeks..')
unique_weeks = calculateUniqueWeeks(grouped_week)

# group unique weeks by user_id
print('[INFO] - Grouping unique weeks by user id..')
grouped_weeks = unique_weeks.groupby('user_id')

# determine home country based on unique weeks
print('[INFO] - Detecting home locations of users based on unique weeks..')
twitter_home_unique_weeks = detectHomeCountry(grouped_weeks,
                                              'home_unique_weeks',
                                              args['column'])

# Save output
print('[INFO] - Saving initial results to pickle..')
twitter_home_unique_weeks.to_pickle(args['output']) 

# generate csv for detected home countries per user
print('[INFO] - Saving csv files of user counts per home location and unique user home locations..')
hom_det = twitter_home_unique_weeks[['user_id','home_unique_weeks']].drop_duplicates()
hom_det.to_csv(args['userhome'], sep=';', encoding='utf-8')

# get user counts per detected country
user_counts_per_country = hom_det.home_unique_weeks.value_counts()
user_counts_per_country.to_csv(args['ucount'], sep=';', encoding='utf-8')

print('[INFO] - ... done!')