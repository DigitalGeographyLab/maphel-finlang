#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

INFO
####

This script plots temporal changes in Shannon-Weiner entropy across areas and
saves the output as a pdf file.


USAGE
#####

Run the script with the following command:
    python temporal_plot.py -i path/to/tweets.pkl -a path/to/areas.gpkg -o path/to/figure.pdf


@author: Tuomas Väisänen
"""

import pandas as pd
import geopandas as gpd
import seaborn as sns
from collections import Counter
from itertools import dropwhile
import skbio.diversity.alpha as sk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from sklearn import preprocessing
import pickle5 as pickle
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define input file path
ap.add_argument("-i", "--input", required=True,
                help="Path to input pickle containing geotagged tweets.")

# Define path to area polygons
ap.add_argument("-a", "--areas", required=True,
                help="Path to input file containing area polygons.")

# Define path for output pdf
ap.add_argument("-o", "--output", required=True,
                help="Path to the output pdf.")

# Parse arguments
args = vars(ap.parse_args())

# function to count languages occurring more than once
def langcount(langlist):
    '''
    
    Parameters
    ----------
    langlist : List of languages.

    Returns
    -------
    A dictionary with languages and counts without singletons.

    '''
    cnt = Counter(langlist)
    for key, count in dropwhile(lambda key_count: key_count[1] >= 2, cnt.most_common()):
        del cnt[key]
    return list(cnt.keys()), list(cnt.values())

# initialize min max scaler
min_max_scaler = preprocessing.MinMaxScaler()

# function to scale values
def scale_minmax(series):
    '''
    Parameters
    ----------
    series : pandas series containing lang counts in a list.

    Returns
    -------
    Min-Max scaled values

    '''
    # get values
    series = series.values
    # reshape for scikit
    series = series.reshape(-1,1)
    # scale
    scaled = min_max_scaler.fit_transform(series)
    return scaled

# read pickle in
print('[INFO] - Reading data in...')
with open(args['input'], 'rb') as handle:
    df = pickle.load(handle)

# read spatial units in
areas = gpd.read_file(args['areas'])

# spatial join attributes
print('[INFO] - Spatial joining...')
tweetdf = gpd.sjoin(df, areas)

# drop geometry and convert to pandas dataframe, otherwise .explode wont work
tweetdf = pd.DataFrame(tweetdf.drop(columns='geometry'))

# extract hours and months
print('[INFO] - Extracting temporal data...')
tweetdf['hour'] = tweetdf['created_at'].dt.hour
tweetdf['month'] = tweetdf['created_at'].dt.month
tweetdf['week'] = tweetdf['created_at'].dt.week

# drop week 53 which is one day and only present in 2018
tweetdf = tweetdf[tweetdf['week'] != 53]

# explode tweets
tweetdf = tweetdf.explode('langs')

# get diversity order
print('[INFO] - Preparing data for plotting...')
divord = tweetdf.groupby('nimi')['langs'].apply(list).rename('langs').reset_index()
divord['counts'] = divord['langs'].apply(lambda x: langcount(x)[1])
divord['shannon'] = divord['counts'].apply(lambda x: np.exp(sk.shannon(x, base=np.e)))
divord = divord.sort_values(by=['shannon'], ascending=False)
divord = divord['nimi'].tolist()

# get unique user counts per spatial unit
users = tweetdf.groupby('nimi')['user_id'].apply(list).rename('users').reset_index()
users['count'] = users['users'].apply(lambda x: len(Counter(x)))
users = pd.Series(users['count'].values,index=users.nimi).to_dict()

# group areas and calculate langs
tweetareas = tweetdf.groupby(['nimi','week'])['langs'].apply(list).rename('langs').reset_index()

# count langs
print('[INFO] - Calculating language counts and Shannon diversity...')
tweetareas['counts'] = tweetareas['langs'].apply(lambda x: langcount(x)[1])

# count shannon entropy and scale
tweetareas['shannon'] = tweetareas['counts'].apply(lambda x: np.exp(sk.shannon(x, base=np.e)))

# add scaled values
tweetareas['shannon_minmax'] = scale_minmax(tweetareas['shannon'])

# get spatial units
spatial_units = divord

# dictionary for english region names
newnames = {'Pohjanmaa':'Ostrobothnia', 'Ahvenanmaa ':'Åland islands',
         'Keski-Pohjanmaa':'Central Ostrobothnia'}

# rename regions
tweetareas = tweetareas.replace({'nimi':newnames})

# special cases
cases = ['Vaasa', 'Pohjanmaa', 'Helsinki', 'Porvoo', 'Ahvenanmaa ', 'Joensuu',
         'Seinäjoki', 'Keski-Pohjanmaa', 'Kuopio', 'Pori', 'Pirkanmaa', 'Oulu']

newcases = ['Vaasa', 'Ostrobothnia', 'Helsinki', 'Porvoo', 'Åland islands', 'Joensuu',
            'Seinäjoki', 'Central Ostrobothnia', 'Kuopio', 'Pori', 'Pirkanmaa', 'Oulu']

# set plotting context
sns.set(font_scale=1.3)
sns.set_style("white", {'axes.grid': False})

#initiate palette
palette = sns.cubehelix_palette(len(cases), start=0.1, rot=.98, reverse=True)

# initiate sub plots
print('[INFO] - Starting to plot figure...')
fig, axs = plt.subplots(4,3, figsize=(15, 17))

# loop over spatial units
for i, ax in zip(range(len(newcases)), axs.flat):
    # fetch spatial unit names
    su = newcases[i]
    fisu = cases[i]
    # get data for current spatial unit
    data = tweetareas[tweetareas['nimi'] == su]
    # plot shannon entropy over chosen grouping type
    sns.regplot(data=data, x='week', y='shannon_minmax', ax=ax,
                order=5, fit_reg=True, x_ci=99.9, n_boot=1000,
                scatter_kws={'s':10}, scatter=True, 
                color=palette[i]).set_title(su)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.text(x=0.98, y=0.91, s='Users: ' + str(users[fisu]), fontsize=11, alpha=0.75, ha='right', va='bottom', transform=ax.transAxes)
    ax.set(ylim=(0, 1))
    ax.set_ylabel('Shannon entropy\n(normalized)', fontsize=11)
    ax.set_xlabel('Week', fontsize=11)
    ax.tick_params(axis='x',color='black')
plt.subplots_adjust(hspace=0.42, wspace=0.36)
plt.savefig(args['output'], dpi=300, bbox_inches='tight')
print('[INFO] - Figure saved to disk!')

print('[INFO] - ... done!')