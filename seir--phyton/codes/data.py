# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 21:55:11 2020
@author: Phurichai Rungcharoenkitkul

-- Modified --
on: Apr 2021
by: Ricardo Aguirre Leal
"""

import pandas as pd
# import numpy  as np
import pickle
import os
from datetime import date
# ------------------
dirpath = os.path.abspath('..')

# =================================================
# Basic settings
# =================================================

chosen_date = '2021-04-21' 
arq_cov     = 'dados_brasil_spline_2021-04-21'
arq_cov_reg = 'dados_brasil_reg_spline_2021-04-21'

# =================================================
# Read datas .csv 
# =================================================

filepathcovid    = os.path.join(dirpath, f"data\\{arq_cov}.csv"    )
filepathcovidreg = os.path.join(dirpath, f"data\\{arq_cov_reg}.csv")
df  = pd.read_csv(filepathcovid   )
dfr = pd.read_csv(filepathcovidreg)

# =================================================
# Adjusts in dfs
# =================================================

newcols  = {'total_vaccinations_per_hundred':      'vac_total'  ,
            'people_vaccinated_per_hundred' :      'vac_partial',
            'people_fully_vaccinated_per_hundred': 'vac_fully'  ,
            'google_retail_and_recreation':        'mobility_L' ,
            'google_grocery_and_pharmacy' :        'mobility_F' ,
            'google_workplaces':                   'mobility_W' ,
            'google_parks':                        'mobility_P' ,
            'google_transit':                      'mobility_T' ,
            'google_residential':                  'mobility_R' }
pivcols  = ['total_cases', 'total_deaths',
            'new_cases'  , 'new_deaths'  , 'population',  
            'mobility_L' , 'mobility_W'  , 'mobility_F',
            'mobility_P' , 'mobility_T'  , 'mobility_R',
            'vac_total'  , 'vac_partial' , 'vac_fully' ]
# ------------------
df  = df.rename(columns=newcols)
df1 = df.pivot_table(values=pivcols, index='date', columns=['iso2'], dropna=False)
df1.index = pd.to_datetime(df1.index)
df1['mobility_L'] = df1['mobility_L']/100
df1['mobility_W'] = df1['mobility_W']/100
df1['mobility_F'] = df1['mobility_F']/100 
df1['mobility_P'] = df1['mobility_P']/100
df1['mobility_T'] = df1['mobility_T']/100
df1['mobility_R'] = df1['mobility_R']/100
# ------------------
dfr  = dfr.rename(columns=newcols)
dfr1 = dfr.pivot_table(values=pivcols, index='date', columns=['iso2'], dropna=False)
dfr1.index = pd.to_datetime(dfr1.index)
dfr1['mobility_L'] = dfr1['mobility_L']/100
dfr1['mobility_W'] = dfr1['mobility_W']/100
dfr1['mobility_F'] = dfr1['mobility_F']/100 
dfr1['mobility_P'] = dfr1['mobility_P']/100 
dfr1['mobility_T'] = dfr1['mobility_T']/100 
dfr1['mobility_R'] = dfr1['mobility_R']/100 

# =================================================
# Save datas into files
# =================================================

namef  = dirpath + "\\data\\data_daily_"     + chosen_date +'.pkl'
namefr = dirpath + "\\data\\data_daily_reg_" + chosen_date +'.pkl'
pickle.dump(df1,  open(namef, 'wb'))
pickle.dump(dfr1, open(namefr,'wb'))

