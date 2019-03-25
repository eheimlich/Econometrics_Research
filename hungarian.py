# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:22:51 2019

@author: Ethan H
"""

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import scipy.optimize

#def hungarian_matching(n, mean1, mean2, sd1, sd2, percentile):

n = 20

# Generate each column of data:
a = np.random.normal(100, 10, n)
b = np.random.normal(100, 10, n)
    
# Create pandas Dataframe:
data = pd.DataFrame({'X1': a,'X2': b})

# Create dissimilarity matrix:
dist = pdist(data, 'euclidean')
df_dist = pd.DataFrame(squareform(dist))

# Create Hungarian Matrix 
hungarian_matrix =  df_dist.iloc[:int(n/2), int((n/2)):]

optimum_pairs = scipy.optimize.linear_sum_assignment(hungarian_matrix)

pairs_dict = {}

for i in range(int(n/2)):
    pairs_dict[i] = optimum_pairs[1][i]

distance_dict = {}

for i in range(int(n/2)):
    distance_dict[i] = hungarian_matrix.iloc[i, optimum_pairs[1][i]]
    
#print(optimum_pairs[1][0])
#print(optimum_pairs)
distance_dict = sorted(distance_dict.items(),key = lambda x : x[1])

    
print (distance_dict)
print(distance_dict[0][0])
print(pairs_dict)
print(pairs_dict[distance_dict[0][0]])
    
    
    