#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

INFO
####

This script calculates diversity metrics of language use for each area and
saves the output as a geopackage file.


USAGE
#####

Run the script with the following command:
    python area_language_diversity.py -i path/to/tweets.pkl -a path/to/areas.gpkg
    -o path/to/output.gpkg


NOTE
####

This script assumes the top twitter languages in Finland and some column names
are hard-coded. Please alter script if using somewhere else.

@author: Tuomas Väisänen
"""

import pandas as pd
import geopandas as gpd
import skbio.diversity.alpha as sk
import numpy as np
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define input file path
ap.add_argument("-i", "--input", required=True,
                help="Path to input pickle containing geotagged tweets.")

# Define user home detection csv path
ap.add_argument("-a", "--areas", required=True,
                help="Path to input file containing area polygons.")

# Define path for output pickle
ap.add_argument("-o", "--output", required=True,
                help="Path to the output geopackage.")

# Parse arguments
args = vars(ap.parse_args())

# read area polygon layer in
print('[INFO] - Reading data in..')
areas = gpd.read_file(args['areas'])

# read tweets in
tweets = pd.read_pickle(args['input'])

# get unique langs per tweet
tweets['langs'] = tweets['langs'].apply(lambda x: list(set(x)))

# explode tweets with multiple langs to independent rows
tweets = tweets.explode('langs')

# get all unique languages to list
langlist = list(set(tweets['langs'].values.tolist()))

# convert tweets to geodataframe
tweets = gpd.GeoDataFrame(tweets, crs='EPSG:4326',
                          geometry=gpd.points_from_xy(tweets.longitude, tweets.latitude))

# reproject tweets to finnish national coordinate system
print('[INFO] - Projecting to national CRS..')
tweets_tm35 = tweets.to_crs(areas.crs)

# spatial join
print('[INFO] - Performing spatial join..')
sjoin = gpd.sjoin(areas, tweets_tm35)

# count languages in each area
langcounts = pd.pivot_table(sjoin, index='NATCODE', columns='langs', aggfunc={'langs':len}).fillna(0)

# clean up dataframe
langcounts.columns = langcounts.columns.droplevel()

# merge language counts to areas polygon
areas = areas.merge(langcounts, how='left',on='NATCODE')

# top 10 language list
top10 = ['fi','en','et','ru','sv','es','ja','fr','pt','de']

# get post counts and sums for whole of finland
postmean = areas['post_count'].mean()
postsum = areas['post_count'].sum()

# loop over areas
print('[INFO] - Counting language proportions in areas..')
for i, row in areas.iterrows():
    # loop over top 10 languages
    for l in top10:
        # get post count per area
        posts = row['post_count']
        # get post count in area per language
        lposts = row[l]
        # get national language mean
        lpostmean = areas[l].mean()
        # get national language sum
        lpostsum = areas[l].sum()
        # set up column names for language
        colname = l + '_prop' # proportion against post counts
        colname2 = l + '_mean_prop' # proprtion against spatial unit means
        colname3 = l + '_sum_prop' # proportion aghainst overall proportion
        colname4 = l + '_lang_prop' # proportion against same language total
        # insert proportion of posts in current language in area
        areas.at[i, colname] = ((int(lposts) / int(posts)) * 100) - 100
        # insert proportion of posts in current language in area against national average
        areas.at[i, colname2] = ((((int(lposts) / int(posts)) * 100) / ((float(lpostmean) / float(postmean)) * 100)) * 100) - 100
        # insert proportion of posts in current language in area against national average of sums
        areas.at[i, colname3] = ((((int(lposts) / int(posts)) * 100) / ((int(lpostsum) / int(postsum)) * 100)) * 100) - 100
        # insert proportion of posts in current language in area against sum of all posts in current language
        areas.at[i, colname4] = (int(lposts) / int(lpostsum)) * 100

# get dominant language from selected columns
areas['propmax'] = areas[['fi_prop','en_prop','et_prop','ru_prop','sv_prop','es_prop','ja_prop','fr_prop','pt_prop','de_prop']].idxmax(axis=1)
areas['mean_propmax'] = areas[['fi_mean_prop','en_mean_prop','et_mean_prop','ru_mean_prop','sv_mean_prop','es_mean_prop','ja_mean_prop','fr_mean_prop','pt_mean_prop','de_mean_prop']].idxmax(axis=1)
areas['sum_propmax'] = areas[['fi_sum_prop','en_sum_prop','et_sum_prop','ru_sum_prop','sv_sum_prop','es_sum_prop','ja_sum_prop','fr_sum_prop','pt_sum_prop','de_sum_prop']].idxmax(axis=1)

# get all language column names
cols = list(areas[langlist].columns)

# loop over areas
print('[INFO] - Calculating diversity metrics per area..')
for i, row in areas.iterrows():
    # get counts of languages
    otus = list(row[cols])
    # drop zeros
    otus = [i for i in otus if i != 0]
    # calculate diversity metrics
    areas.at[i, 'dominance'] = sk.dominance(otus)
    areas.at[i, 'berger'] = sk.berger_parker_d(otus)
    areas.at[i, 'menhinick'] = sk.menhinick(otus)
    areas.at[i, 'singletons'] = sk.singles(otus)
    areas.at[i, 'shannon'] = np.exp(sk.shannon(otus, base=np.e))
    areas.at[i, 'unique'] = sk.observed_otus(otus)

# save to file
print('[INFO] - Saving output geopackage...')
areas.to_file(args['output'], driver='GPKG')

print('[INFO] - ... done!')
