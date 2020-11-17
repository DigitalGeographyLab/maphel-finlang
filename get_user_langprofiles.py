#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 15:19:55 2020

INFO
####

This script reads language detected tweets from a PostgreSQL database table to
a pandas dataframe and saves it locally to disk as a pickled dataframe. 


USAGE
#####

Run the script with the following command:
    python get_user_langprofiles.py -ho your.host.com -db databasename
    -u username -pw password -tb table -o path/to/tweets.pkl

NOTE
####

The output pickle contains language detections per user in a python list.

@author: Tuomas Väisänen
"""

import pandas as pd
import psycopg2
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

# Table in database
ap.add_argument("-tb", "--table", required=True,
                help="Table your data is in")

# Output file
ap.add_argument("-o", "--output", required=True,
                help="Path to output pickle file")
# Parse arguments
args = vars(ap.parse_args())

# Assign arguments to variables
database = args['datanase']
host = args['host']
user = args['user']
pw = args['password']
tablename = args['table']

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

# sql to get user language uses
sql = "SELECT user_id, string_agg(language::character varying,';') as langs"\
    " FROM twitter_histories_full_fin_est_fulltext_multilang"\
        " WHERE prob >= 0.7 "\
            " GROUP BY user_id"

# retrieve data
print("[INFO] - Querying to dataframe...")
df = pd.read_sql(sql, con=con)

# convert languages to list
df['langs'] = df['langs'].apply(lambda x: x.split(';'))

# save to csv
print("[INFO] - Saving to disk...")
df.to_pickle(args['output'])

print("[INFO] - ... done!")
