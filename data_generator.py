# simulate later, hand generate for now
# covariates only include price for the time being, where price = 2 * option ID

import pandas as pd

if __name__ == '__main__':
    df = pd.DataFrame({'time': [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4],
          'option_ID': [1, 2, 3, 5, 6, 2, 3, 4, 6, 1, 2, 1, 3, 5],
          'price': [2, 4, 6, 10, 12, 4, 6, 8, 12, 2, 4, 2, 6, 10],
          'count': [10, 2, 4, 11, 23, 8, 3, 2, 1, 6, 32, 66, 21, 19]})
    df = df.set_index(['time', 'option_ID'])
    df.to_csv("handmade_data.csv")