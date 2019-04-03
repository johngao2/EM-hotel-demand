import pandas as pd
import itertools
import numpy as np
import random
import os

def shuffle_list(item_list, types_list):
    '''Helper function to generate list permutations'''
    
    randomized_list = item_list[:]
    while True:
        random.shuffle(randomized_list)
        for prev_type in types_list:
            if randomized_list == prev_type:
                break
        else:
            return randomized_list

def generate_types(N, items):
    '''Generates N random customer types'''
    
    types_list = []
    while len(types_list) < N:
        new_type = shuffle_list(items, types_list)
        if new_type[0] != 0: # making sure customer wants to at least buy something
            types_list.append(shuffle_list(items, types_list))

    types = np.array(types_list)

    for i in range(0, N):
        irrelevant = 0
        for n in range (0, 16):
            if irrelevant == 1:
                types[i][n] = 0
            elif types[i][n] == 0:
                irrelevant = 1
    
    types = pd.DataFrame({'rank1':types[:,0], 'rank2':types[:,1], 'rank3':types[:,2],
              'rank4':types[:,3], 'rank5':types[:,4], 'rank6':types[:,5],
              'rank7':types[:,6], 'rank8':types[:,7], 'rank9':types[:,8],
              'rank10':types[:,9], 'rank11':types[:,10], 'rank12':types[:,11],
              'rank13':types[:,12], 'rank14':types[:,13], 'rank15':types[:,14],
              'rank16':types[:,15]})
    
    types['cust_type'] = types.index + 1
    types = types.set_index('cust_type')

    types = types.astype(int)
        
    return types

def generate_true_proportions(N):
    '''Creates list of tuples, where first element of tuple is type prob, and second element
    represents the type.'''
    
    probs = np.random.dirichlet(np.ones(N),size=1).tolist()[0]
    
    probs = pd.DataFrame({'prob':probs})
    probs['type'] = probs.index + 1
    probs = probs.set_index('type')
    return probs

def generate_availability(T, n, items):
    '''Generates availability data for products, choosing between 
    2 and 10 items (inclusive) per time period'''
    
    avail_arr = np.zeros((T, n))
    for t in range(0, T):
        num_avail = random.randint(2,10)
        avail_list = random.sample(items, num_avail)
        for item in avail_list:
            avail_arr[t][item-1] = 1
    # convert to pd
    avail = pd.DataFrame({'prod1':avail_arr[:,0], 'prod2':avail_arr[:,1], 'prod3':avail_arr[:,2],
                         'prod4':avail_arr[:,3], 'prod5':avail_arr[:,4], 'prod6':avail_arr[:,5],
                         'prod7':avail_arr[:,6], 'prod8':avail_arr[:,7], 'prod9':avail_arr[:,8],
                         'prod10':avail_arr[:,9], 'prod11':avail_arr[:,10], 'prod12':avail_arr[:,11],
                         'prod13':avail_arr[:,12], 'prod14':avail_arr[:,13], 'prod15':avail_arr[:,14]})
    
    # add time index
    avail['T'] = avail.index + 1
    avail = avail.set_index('T')

    avail = avail.astype(int)

    return avail

def generate_transactions(T, types, avail, lam, probs):
    '''Generates transaction data by generating types based on true probs, and assigning
    their highest choice in each time period as the true transaction'''
    
    true_types = np.random.choice(10, T, p=probs['prob']) # generate true types from probabilities
    true_types = true_types
    true_types

    trans_data = np.zeros(T)

    for t in range(0, T):
        lambda_sample = random.random()
        if lambda_sample > lam:
            trans_data[t] = 0
        else:
            cur_type = types.iloc[true_types[t], :]
            cur_avail = avail.iloc[t, :]
            for n in range(0, 15): # iterate from top choice to bottom
                choice = cur_type[n]
                if choice == 0:
                    trans_data[t] = 0
                    break
                if cur_avail[choice-1] == 1: # this choice is available
                    trans_data[t] = choice
                    break

    trans = pd.DataFrame({'prod_num':trans_data[:]})

    # add time index
    trans['T'] = trans.index + 1
    trans = trans.set_index('T')

    trans = trans.astype(int)
    
    return trans

def simulate_dataset(n, lam, T, N, sample_num):
    '''Generates one dataset given parameters'''
    
    items = [i for i in range(0, n+1)]
    probs = generate_true_proportions(N)
    types = generate_types(N, items)
    avail = generate_availability(T, n, items)
    trans = generate_transactions(T, types, avail, lam, probs)
    
    path = 'data/simulated_data/l' + str(lam) + '/' + str(T) + '/' + str(sample_num) + '/'
    
    os.makedirs(path)
    
    probs_path = path + 'probs.csv'
    types_path = path + 'types.csv'
    avail_path = path + 'avail.csv'
    trans_path = path + 'trans.csv'
    
    probs.to_csv(probs_path)
    types.to_csv(types_path)
    avail.to_csv(avail_path)
    trans.to_csv(trans_path)
    
    print('Finished generating sample', sample_num, 'for T=' + str(T) + ', lambda=' + str(lam))

if __name__ == '__main__':
    # simulation parameters
    n = 15 # num items
    lam = [0.8,0.2] # lambda
    T = [10000, 50000, 100000] # selling season lengths
    N = 10 # number of customer types

    for t in T:
        for lambd in lam:
            for sample_num in range (1, 31):
                simulate_dataset(n, lambd, t, N, sample_num)
                
    print('DONE!')



