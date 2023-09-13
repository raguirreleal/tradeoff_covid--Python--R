# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 10:39:18 2021
@author: Phurichai Rungcharoenkitkul

-- Modified --
on: Apr 2021
by: Ricardo Aguirre Leal
"""

from datetime import datetime
# from datetime import date
from pathlib  import Path
# import numpy  as np
import pandas as pd
import os
# ------------------
dirpath = os.path.abspath('..')

# =================================================
# Basic settings
# =================================================

chosen_date  = '2021-04-21' 
default_maxT = datetime(2021, 4, 21)  
virus_thres  = 50     # Epidemic starting point
pd_reg       = False   # The parameter p_d is by state (False) or by mean region (True)?

# =================================================
# Load data and set folders
# =================================================

df1  = pd.read_pickle(dirpath+f'\\data\\data_daily_{chosen_date}.pkl')
dfr1 = pd.read_pickle(dirpath+f'\\data\\data_daily_reg_{chosen_date}.pkl')
# Percentage of infection cases detected p_d:
df_ifr = pd.read_csv(dirpath+f"\\data\\ifr_brasil_epicovid.csv")
df_ifr.index = df_ifr['iso2']  
# Parameter p_d, conform choose by state or by mean region, and set folders
if pd_reg:
    pd_set = pd.DataFrame(df_ifr['p_notif_region_R23'])
    pd_set = pd_set.rename(columns={'p_notif_region_R23': 'region_R23'})
    # ------------------
    optm_param_folder = dirpath+'\\results_optm\\pd_reg'
    param_save_folder = dirpath+f'\\params\\param_pdreg_{chosen_date}'
    out_save_folder = dirpath+f'\\output\\output_pdreg_{chosen_date}'
else:
    pd_set = pd.DataFrame(df_ifr['p_notif_R23'])
    pd_set = pd_set.rename(columns={'p_notif_R23': 'R23'})
    # ------------------
    optm_param_folder = dirpath+'\\results_optm\\pd_states'
    param_save_folder = dirpath+f'\\params\\param_pdstates_{chosen_date}'
    out_save_folder = dirpath+f'\\output\\output_pdstates_{chosen_date}'
Path(f'{param_save_folder}').mkdir(exist_ok=True)
Path(f'{out_save_folder}').mkdir(exist_ok=True)
Path(f'{optm_param_folder}').mkdir(exist_ok=True)

# =================================================
# Epidemiological assumptions
# =================================================

# Specific days:
IncubeD   = 6     
RecoverID = 10 
RecoverHD = 8    
VentilatedD = 10  # Recovery Time when Ventilated
days_lost_immunity = 365 # 
contagion_to_death = 23
infection_to_death = contagion_to_death - IncubeD  # = 17
# ------------------
r_i  = 2**(1/IncubeD)-1       # Incubation-to-infection transition rate
r_ri = 2**(1/RecoverID)-1     # Rate of recovery not under infection
r_rh = 2**(1/RecoverHD)-1     # Rate of recovery under hospitalization
r_rv = 2**(1/VentilatedD)-1   # Rate of recovery under ventilation
# ------------------
# Percentage of ventilated:
p_v  = 0.035  
# ------------------
# Probabilidade de ser hospitalizado, condicional a ter sido detectado:
# 22252 receiving ongoing care; 87515 in-hospital deaths; 321893 hospital admissions for COVID-19
p_hd = 0.34  # = (22252 + 87515)/321893 
# ------------------
# Percentage of detected cases hospitalized:
if pd_reg:     
    p_h  = pd_set['region_R23'] * p_hd  
else:
    p_h  = pd_set['R23'] * p_hd  
# ------------------
effi_one = 0.5  # vaccine efficacy after one dose
effi_two = 0.95 # vaccine efficacy after two doses

# =================================================
#   Alert current time
# =================================================

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Loading as of", current_time)
print(' ')
