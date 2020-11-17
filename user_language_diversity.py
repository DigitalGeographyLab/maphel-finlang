#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

INFO
####

This script calculates diversity metrics for unique users based on their
language profile and plots Figure 1. The output is saved as a pickled dataframe
and the figure and its parts as pdf files. The pdf files go to current
working directory.


USAGE
#####

Run the script with the following command:
    python user_language_diversity.py -i path/to/tweets.pkl -ho path/to/home_detections.csv
    -o path/to/output.pkl


NOTE
####

This script filters out users whose home location has not been detected to be
in Finland. To change this behavior alter code at lines 134-139 accordingly.

@author: Tuomas Väisänen
"""

import pandas as pd
import skbio.diversity.alpha as sk
from collections import Counter
from itertools import dropwhile
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define input file path
ap.add_argument("-i", "--input", required=True,
                help="Path to input pickle user language profiles.")

# Define user home detection csv path
ap.add_argument("-ho", "--home", required=True,
                help="Path to csv containing home locations for unique user ids.")

# Define path for output pickle
ap.add_argument("-o", "--output", required=True,
                help="Path to the output pickle.")


# Parse arguments
args = vars(ap.parse_args())

# function to count languages occurring more than once
def langcount(langlist):
    # count languages
    cnt = Counter(langlist)
    # remove languages that occur only once
    for key, count in dropwhile(lambda key_count: key_count[1] >= 2, cnt.most_common()):
        del cnt[key]
        
    return list(cnt.keys()), list(cnt.values())

# function to calculate diversity metric
def lang_entropy(dictionary):
    '''
    Info
    ----
    Calculates entropy, diversity class and dominant language following
    Holloway (2018)
    
    Parameters
    ----------
    dictionary : Dictionary containing languages and counts

    Returns
    -------
    A list containing the following diversity metrics:
        [Personal scaled entropy,
         share of dominant language,
         share of 2 most popular languages,
         diversity classification,
         dominant language]
    '''
    # get sum of language uses
    langsum = sum(dictionary.values())
    # get unique language count
    uniques = len(dictionary)
    # set dictionaries up
    origdict = dictionary
    entropydict = {}
    # get dominant language
    domlang = Counter(origdict).most_common()[0][0]
    # create language specific entropies
    if uniques > 1:
        for key, value in origdict.items():
            entropydict[key] = (value / langsum) * (np.log(1 / (value / langsum)))
        # get personal entropy
        entropysum = sum(entropydict.values()) / np.log(uniques)
        # count share of largest language
        share = Counter(origdict).most_common()[0][1] / langsum
        # count share of second largest language
        share2 = share + (Counter(origdict).most_common()[1][1] / langsum)
        # first diversity classification
        if entropysum <= 0.3707:
            entclass = 'Low'
        elif entropysum > 0.3707:
            entclass = 'Moderate'
        # classify high diversity
        if (entclass == 'Moderate') & (entropysum > 0.7414) & (share < 0.45) & (share2 < 0.8):
            entclass = 'High'
    else:
        entropysum = 0.0
        share = 1.0
        share2 = 1.0
        entclass = 'Low'
    return [entropysum, share, share2, entclass, domlang]

# read user language profile data in
print('[INFO] - Reading data in..')
data = pd.read_pickle(args['input'])

# read user home locations
users = pd.read_csv(args['home'])

# join origin data
data = pd.merge(data, users, on='user_id')

# rename home detection column
data = data.rename(columns={'home_unique_weeks':'home_country'})

# filter only users who most likely live in Finland
data = data[data['home_country'].str.contains('Finland')]

# count language use without singletons
print('[INFO] - Calculating language diversities...')
data['ulangs'] = data['langs'].apply(lambda x: langcount(x)[0])
data['counts'] = data['langs'].apply(lambda x: langcount(x)[1])
data = data[data['counts'].map(lambda d: len(d)) > 0] # drop empties if any exist

# calculate diversity metrics
data['dominance'] = data['counts'].apply(lambda x: sk.dominance(x))
data['berger'] = data['counts'].apply(lambda x: sk.berger_parker_d(x))
data['menhinick'] = data['counts'].apply(sk.menhinick)
data['simpson'] = data['counts'].apply(sk.simpson)
data['singles'] = data['counts'].apply(sk.singles)
data['shannon'] = data['counts'].apply(lambda x: np.exp(sk.shannon(x, base=np.e)))
data['unique'] = data['counts'].apply(sk.observed_otus)

# language counts to dictionary
data['langdict'] = data.apply(lambda x: dict(zip(x['ulangs'], x['counts'])), axis=1)

# calculate ellis et al diversity metrics
data['divs'] = data['langdict'].apply(lang_entropy)

# calculate number of sentences per user
data['sents'] = data['langs'].apply(len)

# parse results
for i, row in data.iterrows():
    data.at[i, 'ellis_entropy'] = row['divs'][0]
    data.at[i, 'domshare'] = row['divs'][1]
    data.at[i, 'dom2share'] = row['divs'][2]
    data.at[i, 'div_class'] = row['divs'][3]
    data.at[i, 'domlang'] = row['divs'][4]

# save to pickle
print('[INFO] - Saving to pickle...')
data.to_pickle(args['output'])

# Ellis diversity classes
hcnt = data['div_class'].value_counts()['High']
mcnt = data['div_class'].value_counts()['Moderate']
lcnt = data['div_class'].value_counts()['Low']

# Plot figure 1 as unified and separate files
print('[INFO] - Plotting Figure 1...')

plt.figure(figsize=(20,16))
plt.rcParams.update({'font.size': 17})
ax1 = plt.subplot(3, 3, 1)
ax1 = data['dominance'].plot(kind='hist', bins=40, title='Dominance index', grid='both', weights=np.ones(len(data)) / len(data))
ax2 = plt.subplot(3, 3, 2)
ax2 = data['berger'].plot(kind='hist', bins=40, title='Berger-Parker dominance', grid='both', weights=np.ones(len(data)) / len(data))
ax3 = plt.subplot(3, 3, 3)
ax3 = data['menhinick'].plot(kind='hist', bins=40, title='Menhinick richness index', grid='both', weights=np.ones(len(data)) / len(data))
ax4 = plt.subplot(3, 3, 4)
ax4 = data['unique'].plot(kind='hist', bins=26, title='Unique languages', grid='both', weights=np.ones(len(data)) / len(data))
ax5 = plt.subplot(3, 3, 5)
ax5 = data['shannon'].plot(kind='hist', bins=40, title='Shannon entropy', grid='both', weights=np.ones(len(data)) / len(data))
ax6 = plt.subplot(3, 3, 6)
ax6 = sns.violinplot(x='div_class', y='ellis_entropy', hue='div_class', data=data, ax=ax6, dodge=False)
handles, labels = ax6.get_legend_handles_labels()
labels = ['Low (n=20405)', 'Moderate (n=12841)', 'High (n=76)']
ax6.legend(handles,labels, loc=4, fontsize='medium').set_title('')
ax6.set_title('User diversity classes')
ax6.set_xlabel('')
ax6.set_ylabel('Scaled entropy')
plt.tight_layout()
plt.savefig('Figure1.pdf', bbox_inches='tight')

# actual plot used in NeuphMitt but separate files for each part
plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 30})
ax1 = data['dominance'].plot(kind='hist', bins=40,
                             grid='both', weights=np.ones(len(data)) / len(data))
plt.tight_layout()
plt.savefig('Figure_1a.pdf', dpi=200, bbox_inches='tight')


plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 30})
ax2 = data['berger'].plot(kind='hist', bins=40,
                          grid='both', weights=np.ones(len(data)) / len(data))
plt.tight_layout()
plt.savefig('Figure_1b.pdf', dpi=200, bbox_inches='tight')

plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 30})
ax3 = data['menhinick'].plot(kind='hist', bins=40,
                             grid='both', weights=np.ones(len(data)) / len(data))
plt.tight_layout()
plt.savefig('Figure_1c.pdf', dpi=200, bbox_inches='tight')

plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 30})
ax4 = data['unique'].plot(kind='hist', bins=26,
                          grid='both', weights=np.ones(len(data)) / len(data))
plt.tight_layout()
plt.savefig('Figure_1d.pdf', dpi=200, bbox_inches='tight')

plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 30})
ax5 = data['shannon'].plot(kind='hist', bins=40,
                           grid='both', weights=np.ones(len(data)) / len(data))
plt.tight_layout()
plt.savefig('Figure_1e.pdf', dpi=200, bbox_inches='tight')

plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 30})
ax6 = sns.violinplot(x='div_class', y='ellis_entropy', hue='div_class', data=data, dodge=False)
handles, labels = ax6.get_legend_handles_labels()
labels = ['Low (n=' + str(hcnt) + ')',
          'Moderate (n=' + str(mcnt) + ')',
          'High (n=' + str(lcnt) + ')']
ax6.legend(reversed(handles), reversed(labels), loc=4, fontsize='medium').set_title('')
ax6.set_xlabel('')
ax6.set_ylabel('Scaled entropy')
plt.tight_layout()
plt.savefig('Figure_1f.pdf', dpi=200, bbox_inches='tight')

print('[INFO] - ... done!')