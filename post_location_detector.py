# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:21:47 2019

INFO
####

This script does a spatial join between a point and polygon layer. In case of
Twitter data, this script detects in which polygon (e.g. country, municipality etc)
the posts were made in. Make sure both datasets are in the same coordinate
system before running this.


USAGE
#####

Run the script with the following command:
    python post_location_detector.py -p tweets.gpkg -a areas.gpkg -c name
    -oc post_location -o path/to/output.pkl


NOTE
####

This is the first step in home detection analysis and thus it is recommended
to use administrative borders including maritime borders as this increases
the likelihood of a correct detection.

It is recommended to use country or area id codes in the name column as actual names
can be written differently across various datasets e.g.:
    'Finland' vs 'The republic of Finland'.

@author: Tuomas Väisänen
"""

import geopandas as gpd
import sys
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define input point file path
ap.add_argument("-p", "--points", required=True,
                help="Path to geotagged posts geopackage file.")

# Define input area file path
ap.add_argument("-a", "--areas", required=True,
                help="Path to area polygon file.")

# Define name of column containing area names in area file
ap.add_argument("-c", "--column", required=True,
                help="Location name column in areas polygon.")

# Define column name for location in output file
ap.add_argument("-oc", "--ocolumn", required=False, default='post_country',
                help="Location name column in output file. Default is post_country.")

# Define output file path
ap.add_argument("-o", "--output", required=True,
                help="The path to the output pickle.")

# Parse arguments
args = vars(ap.parse_args())

# read files in
points = gpd.read_file(args['points'])
counties = gpd.read_file(args['areas'])

# check matching CRS
print('[INFO] - Checking whether coordinate systems match between areas and points...')
if points.crs == counties.crs:
    print('[INFO] - Coordinate systems match!')
else:
    sys.exit('[INFO] - Error: Not matching CRS. Make sure the coordinate systems match!')

# drop duplicates
points = points.drop_duplicates(subset='id')

# perform spatial join (points-in-polygon)
print('[INFO] - Performing spatial join...')
joined = gpd.sjoin(points, counties, how='inner', op='intersects')

# rename columns
print('[INFO] - Renaming columns to be more informative')
joined = joined.rename(columns={args['column']:args['ocolumn']})

# saving posts with
print('[INFO] - Saving to pickle..')
joined.to_pickle(args['output'])

# extract points not in any polygon
free_pts = points[~points.index.isin(joined.index)]
print('[INFO] - Amount of points not located within any area polygon: ' + str(len(free_pts)))

print('[INFO] - ... done!')