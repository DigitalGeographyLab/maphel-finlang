#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

INFO
####

This script divides home detection counts between areas for clear and ambiguous
cases for users in moderate and low diversity classes. If a user received an
ambiguous home detection eg:
    
    "Helsinki or Espoo"
    
Both Helsinki and Espoo would get .5 added to their user home count.


USAGE
#####

Run the script with the following command:
    python area_language_diversity.py -i path/to/tweets.pkl -a path/to/areas.gpkg
    -o path/to/output.gpkg


NOTE
####

Please run detect_homes.py and user_language_diversity.py before using this
script, otherwise it will not work. In the article this script was run only
after detecting the home MUNICIPALITY/REGION.

@author: Tuomas Väisänen
"""

import pandas as pd
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define input file path
ap.add_argument("-i", "--input", required=True,
                help="Path to input pickle containing geotagged tweets.")

# Define user home detection csv path
ap.add_argument("-ho", "--home", required=True,
                help="Path to input file containing detected home locations "\
                    "per unique user.")

# Define path for output csv for moderate diversity users
ap.add_argument("-om", "--outmod", required=True,
                help="Path to the output csv containing parsed home detections"\
                    " of users with moderate diversity classification.")

# Define path for output csv for low diversity users
ap.add_argument("-ol", "--outlow", required=True,
                help="Path to the output csv containing parsed home detections"\
                    " of users with low diversity classification.")

# Parse arguments
args = vars(ap.parse_args())

# read data in
print('[INFO] - Reading data in..')
df = pd.read_pickle(args['input'])

# read municipal detections in
homes = pd.read_csv(args['home'], sep=';', encoding='utf-8')

# drop column
homes = homes.drop(columns=['Unnamed: 0'])

# rename column
homes = homes.rename(columns={'home_unique_weeks':'home_municipality'})

# join dataframes
merged = pd.merge(df, homes, on='user_id')

# extract moderate and lows
print('[INFO] - Separating diversity classes...')
mods = merged[merged['div_class'] =='Moderate']
mods = mods.drop_duplicates(keep='first', subset='user_id')
lows = merged[merged['div_class'] == 'Low']
lows = lows.drop_duplicates(keep='first', subset='user_id')


# function to parse counts
def homeparse(df):
    # separate clear detections from ambiguous home detections
    amb = df[df['home_municipality'].str.contains(' or ')]
    clr = df[~df['home_municipality'].str.contains(' or ')]
    # count clear cases
    clr = clr['home_municipality'].value_counts()
    # create lists of ambiguous cases
    amb['home_municipality'] = amb['home_municipality'].str.split(' or ')
    # create counts
    amb['counts'] = amb['home_municipality'].apply(lambda x: 1 / len(x))
    # explode lists
    amb = amb.explode('home_municipality')
    # count values
    amb = amb.groupby('home_municipality')['counts'].sum()
    # join counts
    counts = clr.add(amb, fill_value=0)
    # correct municipality names
    #counts = counts.rename(index=correct)
    return counts

# get counts
print('[INFO] - Parsing home location counts...')
modcounts = homeparse(mods)
lowcounts = homeparse(lows)

# to csv
print('[INFO] - Saving to csv files..')
modcounts.to_csv(args['outmod'], sep=',', encoding='utf-8')
lowcounts.to_csv(args['outlow'], sep=',', encoding='utf-8')

print('[INFO] - ... done!')