#!/usr/bin/env python
# coding: utf-8

# In[155]:


import pandas as pd
import numpy as np
import os
from datetime import timedelta, date

'''This script produces a file which shows transactions. Each row is a transaction
with arrival, depart, and room type.'''

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)

# load raw hotel data
purchases_df = pd.read_csv('data/2018_clean.csv', parse_dates=['ARRIVAL', 'DEPART', 'LOOK_DATE'])
purchases_df = (purchases_df[purchases_df['CANCEL_INDICATOR'] == 0] # look at bookings only
                .drop(['RESNO', 'CANCEL_INDICATOR'], axis=1)
                .groupby('group_id').first() # collapse groups
                .sort_values('LOOK_DATE'))

# initialize some helper vars
df_grouped = purchases_df.groupby(by=['LOOK_DATE', 'ARRIVAL']).count()
periods_per_day = df_grouped['DEPART'].max()

look_start = purchases_df['LOOK_DATE'].min()
look_end = purchases_df['LOOK_DATE'].max()

ssn_start = purchases_df['ARRIVAL'].min()
ssn_end =http://localhost:8888/notebooks/EM_extension/preprocessing/create_transactions.ipynb# purchases_df['ARRIVAL'].max()


# In[157]:


# build empty trans df
intraday_range = range(0, periods_per_day)
look_range = pd.date_range(look_start, look_end)
ssn_range = pd.date_range(ssn_start, ssn_end)

trans_df = pd.DataFrame(index = pd.MultiIndex.from_product([look_range, ssn_range, intraday_range], 
                                                           names=['LOOK_DATE', 'ARRIVAL', 'INTRADAY']),
                        columns=['UNIT'])


# In[ ]:


# fill in purchases
for index, row in purchases_df.iterrows():
    intra_day_counter = 0 # counter var for intra day
    cur_cell = (row['LOOK_DATE'], row['ARRIVAL'], intra_day_counter)
    while(trans_df.loc[cur_cell].notnull()[0]): # find empty cell
        intra_day_counter += 1
        cur_cell = (row['LOOK_DATE'], row['ARRIVAL'], intra_day_counter)
    trans_df.loc[cur_cell] = row['UNIT']
    print(row['LOOK_DATE'])


# In[ ]:


trans_df.to_csv('data/transactions_pre_dates.csv')


# In[210]:


# create alternate version with t seasons and numbered products
trans_df_t = pd.read_csv('data/transactions_pre_dates.csv')
trans_df_t = (trans_df_t.reset_index()
              .drop(['LOOK_DATE', 'ARRIVAL', 'INTRADAY'], axis=1))
trans_df_t.index = trans_df_t.index.rename('T')
trans_df_t['UNIT'] = trans_df_t['UNIT'].map({'CD': 1, 'DD': 2, 'CK': 3, 'DK': 4,
                                             'DKB': 5, '2BV': 6, '4BV': 7})

trans_df_t = trans_df_t.fillna(int(0))
trans_df_t['UNIT'] = trans_df_t['UNIT'].astype(int)


# In[ ]:


trans_df_t.to_csv('data/transactions_pre_t.csv')

