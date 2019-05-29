#!/usr/bin/env python
# coding: utf-8

# In[658]:


import pandas as pd
import numpy as np
import os
from datetime import timedelta, date
import sys

'''This script produces a file which shows availability. Each row is a booking
date, and each col is a (room type, stay date) tuple'''

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)

# load raw hotel data
df_data = pd.read_csv('data/2018_clean.csv', parse_dates=['ARRIVAL', 'DEPART', 'LOOK_DATE'])
df_data = (df_data.sort_values(['LOOK_DATE', 'RESNO'])
           .drop(['RESNO'], axis=1))

df_caps = pd.read_csv('data/capacities.csv', index_col='UNIT')


# In[659]:


# helper vars
look_start = df_data['LOOK_DATE'].min()
look_end = df_data['LOOK_DATE'].max()
look_range = pd.date_range(look_start, look_end)

ssn_start = df_data['ARRIVAL'].min()
ssn_end = df_data['ARRIVAL'].max()
ssn_range = pd.date_range(ssn_start, ssn_end)

# get max periods per day from transactions script
periods_per_day = (df_data[df_data['CANCEL_INDICATOR'] == 0] # look at bookings only
                   .drop(['CANCEL_INDICATOR'], axis=1)
                   .groupby('group_id').first() # collapse groups
                   .sort_values('LOOK_DATE')
                   .groupby(by=['LOOK_DATE', 'ARRIVAL']).count()
                   .max()[0])
intraday_range = range(0, periods_per_day)


# In[660]:


# initializing blank df with same height as transactions

# generate room type list
df_grouped_types = df_data.groupby(by='UNIT').count()
unit_list = df_grouped_types.index.tolist()

# create blank df
df_avail = pd.DataFrame(index=pd.MultiIndex.from_product([look_range, ssn_range, intraday_range], 
                                                           names=['LOOK_DATE', 'ARRIVAL', 'INTRADAY']),
                        columns=unit_list)
df_avail = df_avail.fillna(0.0)

# add default capacities
for index, row in df_capa.iterrows():
    df_avail[index] = df_avail[index] + row['CAPACITY']


# In[662]:


row = df_data.iloc[1,:]
df_test = df_avail.copy()
cur_idx = (row['LOOK_DATE'], row['ARRIVAL'], 1)
cur_idx2 = (row['LOOK_DATE'], row['ARRIVAL'])
# df_test.loc[cur_idx2, row['UNIT']]


# In[663]:


len(df_data)


# In[664]:


# fill in capacities by iterating over transaction data

# helper load bar function
def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()
    
total = len(df_data)
    
for index, row in df_data.iterrows():
    progress(index, total, status='Filling in availability')  
    cap_change = 1/row['grp_size']
    
    # book arrive delta, must subtract capacity for all these dates
    ba_delta = pd.date_range(row['LOOK_DATE'], row['ARRIVAL'])

    # subtract or add capacity for each day b/w look and arrive
    if row['CANCEL_INDICATOR'] == 0:
        # loop subtract (i know this is bad practice lol)
        for date in ba_delta:
            for i in intraday_range:
                cur_idx = (date, row['ARRIVAL'], i)
                df_avail.loc[cur_idx, row['UNIT']] -= cap_change
    else:
        for date in ba_delta:
            for i in intraday_range:
                cur_idx = (date, row['ARRIVAL'], i)
                df_avail.loc[cur_idx, row['UNIT']] += cap_change


# In[677]:


df_avail.to_csv('data/availability_general.csv')


# In[714]:


# collapse index, map to 1's and 0's and prod nums
df_avail_t = pd.read_csv('data/availability_general.csv')

df_avail_t = df_avail_t.drop(['LOOK_DATE', 'ARRIVAL', 'INTRADAY'], axis=1)
df_avail_t.index = df_avail_t.index.rename('T')
df_avail_t = df_avail_t.applymap(lambda x: 0 if x <= 0 else 1)

df_avail_t = df_avail_t.rename(index=str, columns={"CD": "prod_1",
                                                   "DD": "prod_2",
                                                   "CK": "prod_3",
                                                   "DK": "prod_4",
                                                   "DKB": "prod_5",
                                                   "2BV": "prod_6",
                                                   "4BV": "prod_7"})

df_avail_t = df_avail_t.reindex(sorted(df_avail_t.columns), axis=1)


# In[719]:


df_avail_t.to_csv('data/availability_sprint1.csv')


# In[720]:


df_avail


# In[ ]:




