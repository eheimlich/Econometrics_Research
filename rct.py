# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:19:44 2019

@author: Ethan H
"""

import hungarian

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import scipy.optimize

from sklearn.linear_model import LogisticRegression

from sklearn.cross_validation import train_test_split


class RCT:
    
        
    def __init__(self, n): 
        self.n = n
        a = np.random.normal(100, 10, self.n)
        b = np.random.normal(100, 10, self.n)
        assignment = self.treatment_assignment()
        
        
        # Create pandas Dataframe:
        self.data =  pd.DataFrame({'Assignment': assignment, 'X1': a ,'X2': b})
        
    
    def treatment_assignment(self):
        """
        TODO: May want to change this so that the treatment and control ratio
        is balanced.
        
        Returns a treatment assignment vector.
        @rtype: ndarray
        """
        return(np.random.binomial(1, .5, self.n))
        
    def generate_prop(self, x_features):
        """
        Generates propensity scores.
        @rtype: ndarray
        """
        x = self.data[x_features] # Features
        y = self.data['Assignment'] # Target variable
        
        #Split data to have all in the training set. 
        x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=1)
   
        # fit the model with data
        logreg = LogisticRegression()

        logreg.fit(x_train, y_train)

        y_pred=logreg.predict_proba(self.data[x_features])
   
        return(y_pred[:,1])
        



