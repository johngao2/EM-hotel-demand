import pandas as pd
import numpy as np
import os

print(pd.__path__)
print(os.getcwd())

'''This script produces a file which shows room availability. Each column is
a date/room type tuple, and each row shows all the room type/date pairs
that have available spots'''

# load raw hotel data
df = pd.read_csv('data/2017_booking_data.csv')
df = df.drop(['RESNO'], axis=1)

print(df)