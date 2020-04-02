#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:14:31 2020

@author: waeiski
"""

import pandas as pd
import psycopg2
import fasttext
from urllib.parse import urlparse
import emoji
import re
from sqlalchemy.engine.url import URL
from sqlalchemy import create_engine
from sqlalchemy import MetaData
from sqlalchemy.orm import sessionmaker
import argparse

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define the path to input database
ap.add_argument("-ho", "--host", required=True,
                help="Address of your database server")

# Name of the database to connect to on host server
ap.add_argument("-db", "--database", required=True,
                help="Database name")

# User name for the database on the host server
ap.add_argument("-u", "--user", required=True,
                help="Your username in the database")

# Password for user
ap.add_argument("-pw", "--password", required=True,
                help="Your password in the database")

# Define the preprocessing strategy
ap.add_argument("-p", "--preprocessing", required=True,
                help="Selected preprocessing strategy: valid values include "
                     "'no_preprocessing', 'rm_all' and 'rm_trail'.")

# Define input column manually
ap.add_argument("-c", "--column", required=False,
                help="The name of the column containing the texts to process.")

# Define input table name
ap.add_argument("-i", "--input", required=True,
                help="The name of the input table in database.")

# Define output table name
ap.add_argument("-o", "--output", required=True,
                help="The name of the output table in database.")

# Define output table name
ap.add_argument("-mp", "--modelpath", required=True,
                help="The path to fastText model.")

# Parse arguments
args = vars(ap.parse_args())

# Assign arguments to variables
prep = args['preprocessing']
database = args['database']
host = args['host']
user = args['user']
pw = args['password']
tablename = args['input']
outtable = args['output']
modelpath = args['modelpath']

# Check if DataFrame input column has been set manually
if args['column'] is not None:
    inputcol = args['column']
else:
    inputcol = 'text'
    
ft_model = fasttext.load_model(modelpath)

# Define the preprocessing function
def preprocess_caption(row, mode):
    """Applies the selected preprocessing steps to the text.

     Args:
         row: A UTF-8 string.
         mode: A string indicating the selected preprocessing strategy.
               Valid values include: 'no_preprocessing' (no preprocessing),
               'rm_all' (remove all hashtags) and 'rm_trail' (remove trailing
               hashtags).

     Returns:
         A string containing the preprocessed text.
    """
    # Check if preprocessing has been requested.
    if mode != 'no_preprocessing':

        # Convert unicode emoji to shortcode emoji
        row = emoji.demojize(row)

        # Remove single emojis and their groups
        row = re.sub(r':(?<=:)([a-zA-Z0-9_\-&\'’]*)(?=:):', '', row)

        # Apply the selected preprocessing strategy defined in the variable
        # 'mode'. This defines how each row is processed. The selected mode
        # defines the preprocessing steps applied to the data below by
        # introducing different conditions.

        # Remove all mentions (@) in the caption
        row = re.sub(r'@\S+ *', '', row)

        # If mode is 'rm_all', remove all hashtags (#) in the caption
        #if mode == 'rm_all':
        #    row = re.sub(r'#\S+ *', '', row)
        
        # remove hashes from hastags
        row = row.replace('#', '')
        
        # remove old school heart emojis <3
        row = row.replace('&lt;3', '')
        
        # remove greater than symbols >
        row = row.replace('&gt;', '')

        # Split the string into a list
        row = row.split()

        # Remove all non-words such as smileys etc. :-)
        row = [word for word in row if re.sub('\W', '', word)]

        # Check the list of items for URLs and remove them
        try:
            row = [word for word in row if not urlparse(word).scheme]
        except ValueError:
            pass

        # Attempt to strip extra linebreaks following any list item
        row = [word.rstrip() for word in row]

        # If mode is 'rm_trail', remove hashtags trailing the text, e.g.
        # "This is the caption and here are #my #hashtags"
        if mode == 'rm_trail':
            while len(row) != 0 and row[-1].startswith('#'):
                row.pop()

        # Reconstruct the row
        row = ' '.join(row)

        # If mode is 'rm_trail', drop hashes from any remaining hashtags
        if mode == 'rm_trail':
            row = re.sub(r'g*#', '', row)

    # Simplify punctuation, removing sequences of exclamation and question
    # marks, commas and full stops, saving only the final character
    row = re.sub(r'[?.!,_]+(?=[?.!,_])', '', row)

    # Return the preprocessed row
    return row

def detect_ft(caption, preprocessing):
    """Identifies the language of a text using fastText.

    Args:
        caption: A string containing UTF-8 encoded text.
        preprocessing: A string indicating the selected preprocessing strategy.
                       Valid values include: 'no_preprocessing'
                       (no preprocessing),  'rm_all' (remove all hashtags) and
                       'rm_trail' (remove trailing hashtags).

    Returns:
        Saves the prediction into a column named 'langid' in the pandas
        DataFrame as a list of three tuples. The three tuple consists of an
        ISO-639 code, its associated probability and character length of the
        string input to fastText, e.g. ('en', 0.99999, 21).
    """
    # If the caption is None, return None
    if caption == 'None' or caption is None:
        return

    # Preprocess the caption
    caption = preprocess_caption(caption, preprocessing)

    # Perform sentence splitting for any remaining text
    if len(caption) == 0:
        return None

    else:
        # Calculate the character length of each sentence
        char_len = len(caption)

        # Make predictions
        prediction = ft_model.predict(caption)

        # Get the predicted languages and their probabilities
        language = prediction[0][0][-2:]
        probability = prediction[1][0]

        # Return languages and probabilities
        return [language, probability, char_len]
    
# Database info
print("[INFO] - Setting up database URL...")
db_url = URL(drivername='postgresql+psycopg2', host=host, database=database,
                   username=user, port=5432, password=pw)

# Create engine
print("[INFO] - Creating database engine...")
engine = create_engine(db_url, use_batch_mode=True)

# set up database connection
con = psycopg2.connect(database=database, user=user, password=pw,
                       host=host)
print('[INFO] - Connected to ' + str(database) + ' at ' + str(host))

# Init Metadata
meta = MetaData()

# Create session
print("[INFO] - Launching database session...")
Session = sessionmaker(bind=engine)
session = Session()

# sql to count rows
sql_rc = "SELECT COUNT(*) FROM " + tablename

# read user list in
userdf = pd.read_excel('UserList_FinEst.xlsx')

# define userlist
userlist = userdf['userID'].values.tolist()

# List of common bots
botlist = [61043461,61043172,126049550,618294231,1148156568,1447948944,
           756175396031397888,3048544857,429803867,248121914,1216653738,
           1959974761,3207467247,4130482420,2831214083,3291286474,210777901,
           844600442189332480,214474283,6456652,2436084463,3306180851,
           3291649719,782268722006425600,4829228050,912447782,2202066812,
           910848812431937536,1075013455,3109126514,416661318,808745240089808896,
           322598698]

# Define chunksize
chunksize = 500000

# Get row count
print('[INFO] - Getting row count of ' + str(tablename) + '...')
row_count = int(pd.read_sql(sql_rc, con).values)
print('[INFO] - Row count: ' + str(row_count))

# Get approximation of detection runs
run_count = int(row_count / chunksize) + 1
print('[INFO] - Table chunk count: ' + str(run_count))

# Read chunks in, detect languages and update database
for i in range(int(row_count / chunksize) + 1):
    print('[INFO] - Querying chunk {}/{}'.format(i, run_count))

    query = "SELECT * FROM {tablename} ORDER BY row_id LIMIT {chunksize} OFFSET {offset};".format(
        tablename=tablename, chunksize=chunksize, offset=i * chunksize)

    # read queary into dataframe
    df = pd.read_sql(query,con=con)

    # drop known bots
    print('[INFO] - Removing bots from chunk {}'.format(i))
    df = df[~df['user_id'].isin(botlist)]
    
    # dropping unwwanted users
    #print('[INFO] - Selecting identified users from chunk {}'.format(i))
    df = df[df['user_id'].isin(userlist)]
    
    # detect languages
    print('[INFO] - Detecting languages in chunk {}'.format(i))
    df['langid'] = df[inputcol].apply(lambda x: detect_ft(x, prep))
    
    # drop rows without language detections
    print('[INFO] - Dropping rows without language detections')
    df = df[df['langid'].notnull()]
    
    # parse results
    print('[INFO] - Parsing results to improve readability')
    df['language'] = df['langid'].apply(lambda x: x[0])
    df['prob'] = df['langid'].apply(lambda x: x[1])
    df['charlen'] = df['langid'].apply(lambda x: x[2])
    
    # drop list, postgreSQL doesn't support it
    df = df.drop(columns=['langid'])
    
    # push language detection results to database
    print('[INFO] - Pushing results from chunk {} results to table'.format(i))
    df.to_sql(outtable, schema='public', con=engine, if_exists='append',
              index=False, method='multi')

print('[INFO] - Table ' + str(tablename) + ' processed and pushed to ' + str(outtable))
print('[INFO] - ...done!')
# Close session
session.close()