# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 02:27:40 2024

@author: GBAGUIDI Mannondé D.
"""

## Necessary libraries are imported:The pandas, numpy, and scipy libraries are imported.
## pandas for data manipulation, numpy for mathematical operations,
## rankdata for assigning ranks to data, and chi2 for conducting the chi2 test.
##Impoertation of lidrarie math to help to clean the data.
import pandas as pd
import numpy as np
from scipy.stats import rankdata, chi2
from collections import Counter
import math


def rank_correction(rm):
  """Calculates the rank correction factor for a list of ranks.

      rm: A list of ranks.

      The rank correction factor, a float between 0 and 1. A value closer to 1
          indicates less correction is needed.
 
    # Example usage
          ranks = [1, 2, 2, 3, 1]
          correction_factor = rank_correction(ranks)
          print(correction_factor)  # Output: 0.9
  """
 
  
  sv = sorted(rm)  # Sort the input list
  ct = Counter(sv)  # Create a dictionary of rank counts
  n = len(sv)  # Get the number of elements
  corr = 0.0  # Initialize correction factor

  for value, count in ct.items():
    if count > 1:
      corr += (count**3 - count)  # Update correction for tied ranks

  if n < 2:
    return 1.0  # If less than 2 elements, correction factor is 1

  else:
    return 1.0 - corr / (n**3 - n)  # Calculate rank correction factor


## Definition of a function KW_test(data) which takes a pandas DataFrame as parameter.

def KW_test(data, interpretation = True, alpha = 0.05):
    
    """Compute the Kruskal-Wallis H-test for independent samples.

    The Kruskal-Wallis H-test tests the null hypothesis that the population
    median of all of the groups are equal.  It is a non-parametric version of
    ANOVA.  The test works on 2 or more independent samples, which may have
    different sizes.  Note that rejecting the null hypothesis does not
    indicate which of the groups differs.  Post hoc comparisons between
    groups are required to determine which groups are different. This function,
    performs the Kruskal-Wallis test on a DataFrame to compare medians
     across groups.
 Arguments:
       data: A pandas DataFrame containing the data for the test.
       This data will have 2 or more columns.
       
       alpha:(Optional) To give a interpretation, this parameter is use like
       the value of the error. The value of error depend of the domaine
       of application.The following options are available (default is : 0.05)
       
       interpretation:(Optional)
        Give a intrepretion of the result.
        The following options are available (default is 'True'):
            * True: print in the terminal a interpretation of the result
            * False: this option dont give the interpretation of the result
return:
     This fonction return The p-value of the Kruskal-Wallis test,
     a float between 0 and 1.
    """
    ## Check validity of interpretation and alpha   
    if interpretation not in (True, False):
        raise ValueError("interpretation must be True or False")
    if alpha>0 and alpha<1:
        pass
    else:
        raise ValueError("alpha is a float number betwen 0 and 1, like: 0.05")


    ## The number of groups or variables in the DataFrame is retrieved and stored in `p`.
    ## A empty list `n` is initialized to store the number of observations per group.
    ## The total sum of observations `N` is initialized to 0.
    p = len(data.columns.tolist())   
    n = [] 
    N = 0
    
    ## Check for at least two group
    if p < 2:
        raise ValueError("Need at least two groups in stats.kruskal()")
    
    ## Iterate through columns (groups) to count non-missing values
    for i in range(1, p+1) :
        ## For each column in the DataFrame, the number of non-NaN (Not a Number) observations
        ## is counted and stored in `n`, and the total number of non-NaN observations is stored in `N`.
        n = n + [len(list(filter(lambda x: not isinstance(x,float) or not math.isnan(x),data[data.columns.tolist()[i-1]])))]
        N += len(list(filter(lambda x: not isinstance(x,float) or not math.isnan(x),data[data.columns.tolist()[i-1]]))) 
    
    # Check for at least one data point in each group
    for i in range(1, p+1):
        if n[i-1]==0:
            raise ValueError("Need at least one data in each groups")
   
    # Collect all data points from non-missing values
    Alldata = []
    for i in range(1, p+1) :
    ## All observations are extracted into a list called `Alldata`.
        Alldata  = Alldata  + list(filter(lambda x: not isinstance(x,float) or not math.isnan(x),data[data.columns.tolist()[i-1]]))
    
    ## The mean ranks (`rm`) of the observations are calculated using the `rankdata` function from scipy.
    rm = rankdata(Alldata, method = 'average')
    R = []
    start_index = 0
    dic = {}
    for taille in n :
        R.append(np.sum(rm[start_index:start_index+taille]))
        start_index += taille
    ## The sum of ranks for each group is calculated and stored in a list called `R`.
    
    ## A dictionary `dic` is created to store the sum of ranks and the number of observations for each group.
    for i in range(1 , p+1) :
        dic[f"V{i}"] = [R[i-1], n[i-1]]

    ## The Kruskal-Wallis correction constant (`factor`) is calculated.
    factor = 12/(N*(N+1))
    
    ## The weighted sum of squares of ranks is calculated.
    S = 0
    for j in range(1, p+1) :
        V = dic[f'V{j}']
        S = S + (V[0]*V[0])/V[1]
    
    ## The Kruskal-Wallis statistic (`H`) is calculated.
    H = factor*S-3*(N+1)
    H /= rank_correction(rm)
  
    ## The degrees of freedom (`df`) are calculated.
    df = p-1
    
    ## The p-value is calculated from the statistic `H` and the degrees of freedom using the chi-squared distribution function.
    p_value = 1 - chi2.cdf(H, df)
    
    ## Interpretation of results
    if interpretation:
        
        print('               Kruskal-Wallis   ')
        print('             ___________________');
        print('             ___________________');
        print('\n \n')
        print(f'H = {H}, df = {df}, p-value = {p_value} \n')
    
        if (p_value <= alpha):
            print(' Alternative hypothesis: True \n')
            print(' So the medians are not all equal. \n')
        else:  
            print(' Null hypothesis: True \n')
            print(' So the medians are all equal. \n')
        
    return p_value

