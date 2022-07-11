# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 18:37:02 2021

@author: jcholbro
"""
import os
import pandas as pd 
import isbnlib 
import numpy as np
import time
import requests
from bs4 import BeautifulSoup
import datetime
import re
from sympy import symbols, Eq, solve


data = pd.read_excel(r"C:\Users\jcholbro\OneDrive - University Of Houston\School_UH\Macro II\Macro_PS8_Data.xlsx", ignore_index=True)

data = data.drop(columns='Unnamed: 2')

mu = np.mean(data.Consumption)
delta = np.std(data.Consumption)
ave_stocks = np.mean(data.Equity)
ave_bonds = np.mean(data.Bonds)
true_premium = ave_stocks-ave_bonds

lambda1 = 1+mu+delta
lambda2 = 1+mu-delta

data['lag_consumption'] = data['Consumption'].shift(1)

cov = data[['Consumption','lag_consumption']].cov()
cov = cov.iloc[0,1]
var = data['Consumption'].var()

rho=(cov+var)/(2*var)
phi=(1+rho)/2

#phi = something
phi1 = 1-phi
phi2 = phi


beta = 0.99
sigma = 10

#x, y = symbols('x y')
w1, w2 = symbols('w1 w2')


# =============================================================================
# w1 = 20
# w2=19
# =============================================================================





eq1 = beta*(phi1*w1*lambda1**(1-sigma)+phi2*w2*lambda2**(1-sigma))+beta*(phi1*lambda1**(1-sigma)+phi2*lambda2**(1-sigma))-w1
eq2 = beta*(phi2*w1*lambda1**(1-sigma)+phi1*w2*lambda2**(1-sigma))+beta*(phi2*lambda1**(1-sigma)+phi1*lambda2**(1-sigma))-w2
solve((eq1,eq2), (w1, w2))
solution = solve((eq1,eq2), (w1, w2))

s1=solution[w1]
s2=solution[w2]
returns_stocks = ((lambda1*(s1+1))/s1 + (lambda2*(s2+1))/s2)*(.5)

# =============================================================================
# Part A
# =============================================================================
returns_stocks


returns_bonds = 1/(beta*(phi1*lambda1**(-sigma)+(phi2*lambda2**(-sigma))))

premium  = (returns_stocks-1)-(returns_bonds-1)
print(premium*100)
how_close = true_premium-premium
print(true_premium)
print(how_close)
# =============================================================================
# Part B
# =============================================================================
returns_bonds

sigma_list = list(range(1,11))
beta_list = list(np.linspace(0.90,1.0,11))
delta_list = list(np.linspace(0.25,1.25,11))

results_rows = []

w1, w2 = symbols('w1 w2')

for beta in beta_list:
    for sigma in sigma_list:
        for delta in delta_list:
            results_col= {}

            #print(delta)
            #print(sigma)
            #print(beta)
            
            M1 = ((lambda1-delta)/(1-(delta/lambda1)))**(-sigma)
            M2 = ((lambda2-delta)/(1-(delta/lambda2)))**(-sigma)
            
            eq3 = beta*(phi1*w1*lambda1*M1+phi2*w2*lambda2*M2)+beta*(phi1*lambda1*M1+phi2*lambda2*M2)-w1      
            eq4 = beta*(phi2*w1*lambda1*M1+phi1*w2*lambda2*M2)+beta*(phi2*lambda1*M1+phi1*lambda2*M2)-w2
            solve((eq3,eq4), (w1, w2))
            solution2 = solve((eq3,eq4), (w1, w2))
            s1=solution2[w1]
            s2=solution2[w2]
            
            returns_stocks = ((lambda1*(s1+1))/s1 + (lambda2*(s2+1))/s2)*(.5)
            
            returns_bonds = 1/(beta*(phi1*lambda1**(-sigma)+(phi2*lambda2**(-sigma))))
            
            premium  = returns_stocks-returns_bonds
            #print(premium)
            how_close = true_premium-premium
            
            results_col['model premium'] = premium*100
            results_col['data_premium'] = true_premium*100
            results_col['difference'] = how_close
            results_col['parameters'] = [beta,sigma, delta]
            results_rows.append(results_col)
            
            
table = pd.DataFrame(results_rows)         
                
#table.dtypes    
table['difference'] = pd.to_numeric(table['difference'])     
table['difference'] = table['difference'].astype(float)
table['difference'].idxmin()
table.loc[table['difference'].idxmin()]


