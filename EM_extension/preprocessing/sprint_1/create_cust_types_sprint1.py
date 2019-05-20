#!/usr/bin/env python
# coding: utf-8

# In[90]:


import pandas as pd


# In[91]:


# blank df
cust_types = pd.DataFrame(columns=['rank_' + str(i) for i in range(1,9)])
cust_types['cust_type'] = ''
cust_types = cust_types.set_index('cust_type')

# add customer types

# flexible customer, prefers cheap
cust_types = cust_types.append({'rank_1': 'CD', 'rank_2': 'CK', 'rank_3': 'DD',
                   'rank_4': 'DK', 'rank_5': 'DKB'}, ignore_index=True)

# flexible customer, prefers expensive
cust_types = cust_types.append({'rank_1': 'DKB', 'rank_2': 'DK', 'rank_3': 'DD',
                   'rank_4': 'CK', 'rank_5': 'CD'}, ignore_index=True)

# double only customer, prefers cheap
cust_types = cust_types.append({'rank_1': 'CD', 'rank_2': 'DD'}, ignore_index=True)

# double only customer, prefers expensive
cust_types = cust_types.append({'rank_1': 'DD', 'rank_2': 'CD'}, ignore_index=True)

# king only customer, prefers cheap
cust_types = cust_types.append({'rank_1': 'CK', 'rank_2': 'DK', 'rank_3': 'DKB'}, ignore_index=True)

# king only customer, prefers expensive
cust_types = cust_types.append({'rank_1': 'DKB', 'rank_2': 'DK', 'rank_3': 'CK'}, ignore_index=True)

# 2BV customer
cust_types = cust_types.append({'rank_1': '2BV'}, ignore_index=True)

# 4BV customer
cust_types = cust_types.append({'rank_1': '4BV'}, ignore_index=True)

cust_types = cust_types.fillna(0)
cust_types.index += 1
cust_types.index.name = 'cust_type'


# map to numbers
cust_types = cust_types.apply(lambda x: x.map({'CD': 1, 'DD': 2, 'CK': 3, 'DK': 4,
                                  'DKB': 5, '2BV': 6, '4BV': 7, 0: 0}),
                 axis=1)

# output
cust_types.to_csv('data/types_sprint1.csv')


# In[ ]:




