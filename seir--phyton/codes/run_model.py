# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 16:15:14 2021
@author: Ricardo
"""

from joblib import Parallel, delayed
import winsound
import pandas as pd 
import numpy  as np
import time
import os
# ------------------
import data                                                                                   
from seir import *
# ------------------
tic = time.perf_counter()
np.warnings.filterwarnings('ignore') 
plt.style.use('ggplot')

# =================================================
# Choosing and setting run
# =================================================

# Fit the mobiblity (mob)? Plot graphics (graph)? 
# What the number of cores for parallel processment (ncores)?  
mob     = True 
graph   = True
ncores  = 4
# Run for State and/or Region?
states  = True
regions = True
# Save and rewrite new results?
save    = True
# ------------------
csetreg = ['Norte', 'Nordeste', 'Centro-Oeste', 'Sudeste', 'Sul', 'Brasil']
cset    = ['AC', 'AM', 'AP', 'CE', 'DF', 'GO', 'MG', 'MS', 'MT', 'PA', 'PE', 'PR', 'RN', 
           'RO', 'RR', 'SE', 'TO','AL', 'BA', 'ES', 'MA', 'PB', 'PI', 'RJ', 'RS', 'SC', 'SP']  

# =================================================
# Deffining functions 
# =================================================

def for_base(c, r_dth, r_d, p_d, mob, lreg):
    tmp = solveCovid(c, r_dth, r_d, p_d)
    tmp.prelim(lreg)
    tmp.compute()
    if mob:
        tmp.fitmob()
        p_dict = np.append(tmp.best_betas, [1e9*tmp.phiL, 1e9*tmp.phiW, 1e9*tmp.phiF])              
        ret = [tmp, p_dict]
    else:           
        ret = [tmp]
    return ret

def run_set(cr_set, mob, cset_name):
    if (cset_name != 'csetreg'):
        set_temp = cr_set[cr_set.iso2.isin(cset)]
        lreg = False
    else:
        set_temp = cr_set[cr_set.iso2.isin(csetreg)]
        lreg = True
    cset_run = set_temp.iso2
    if not cr_set.iso2.isin(cset_run).all():
        cset_name = 'sparse'    
    pd_cj  = set_temp.p_d
    dth_cj = set_temp.r_dth
    d_cj   = set_temp.r_d
    dth_cj = 2**(1/dth_cj)-1
    d_cj   = 2**(1/d_cj)-1
    Is = cset_run.index.values
    # ------------------ 
    results = Parallel(n_jobs=ncores)(delayed(for_base)(cset_run[i], dth_cj[i], d_cj[i], pd_cj, mob, lreg) for i in Is)
    # ------------------ 
    if save:
        for j in np.arange(len(results)):
            results[j][0].df3.to_csv(f'{out_save_folder}\\df3_{results[j][0].iso2}.csv')   
        if mob:
            p_dict = {'Parameters': ['beta0','beta1L','beta2L','beta1W','beta2W','beta1F','beta2F','phiL','phiW','phiF']}
            for j in np.arange(len(results)):  
                p_dict[results[j][0].iso2] = results[j][1]
            pd.DataFrame(p_dict).to_csv(f'{param_save_folder}\\param_est_for_{cset_name}.csv', float_format='%.4f',index=False)
    return results

# =================================================
# Run model
# =================================================

if os.path.exists(f'{optm_param_folder}\\df_r_csetreg.csv') and regions:
    cr_set = pd.read_csv(f'{optm_param_folder}\\df_r_csetreg.csv', header = 0)
    cr_set.index = cr_set['iso2']
    results_reg = run_set(cr_set, mob, 'csetreg')
    if graph:
        for i in np.arange(len(results_reg)):
            results_reg[i][0].graphics(mob)
        print("Graphics for csetreg")

if states:
    results_states = dict(cset0='', cset1='', cset2='')
    ln_cset = list(results_states)   
    for iset in ln_cset:
        if os.path.exists(f'{optm_param_folder}\\df_r_{iset}.csv'):
            cr_set = pd.read_csv(f'{optm_param_folder}\\df_r_{iset}.csv', header = 0)
            cr_set.index = cr_set['iso2']
            results_states[iset] = run_set(cr_set, mob, iset)
            if graph:
                for i in np.arange(len(results_states[iset])):
                    results_states[iset][i][0].graphics(mob)
                print(f'Graphics for {iset}')

plt.show()

# =================================================
#   Alert finalization and total time
# =================================================

winsound.Beep(2500, 100)
toc = time.perf_counter()
print(' ')
print(f'Calculation done || Total Time: {(toc - tic)/60:0.2f} minutes')
print(' ')

# =================================================
#   Mannualy run model
# =================================================

if False:  
    # --- Setting
    c          = 'AL'
    days_r_dth = 1
    days_r_d   = 1
    p_d        = 0.16314
    # --- Execution
    p_dx = pd_set[pd_set.columns[0]]
    p_dx.at[c] = p_d
    ret = for_base(c, (2**(1/days_r_dth)-1), (2**(1/days_r_d)-1), p_dx, True, False)
    ret[0].graphics(mob)
    
# --- OPTIMAL PARAMETERS ESTIMATES:
#                 p_d_name    p_d    r_dth  r_d
# ---  csetreg  (same for R23)       
# Norte         region_R23  0.161859   1.0  1.0
# Nordeste      region_R23  0.163974   1.0  1.0
# Centro-Oeste  region_R23  0.431995   1.0  3.0
# Sudeste       region_R23  0.204892   1.0  1.0
# Sul           region_R23  0.215327   1.0  1.0
# Brasil        region_R23  0.184704   1.0  1.0
# -----  cset0  (same for R23)
# AC            region_R23  0.161859   1.0  5.0
# AL            region_R23  0.163974   1.0  1.0
# AM            region_R23  0.161859   1.0  1.0
# AP            region_R23  0.161859   1.0  1.0
# CE            region_R23  0.163974   1.0  1.0
# DF            region_R23  0.431995   1.0  3.0
# GO            region_R23  0.431995   1.0  3.0
# MG            region_R23  0.204892   1.0  2.0
# MS            region_R23  0.431995   1.0  7.0
# MT            region_R23  0.431995   1.0  4.0
# PB            region_R23  0.163974   1.0  1.0
# -----  cset1  (for R23, change only in RR, whith r_d=3)
# PA            region_R23  0.161859   1.0  1.0
# PE            region_R23  0.163974   1.0  1.0
# PR            region_R23  0.215327   1.0  2.0
# RN            region_R23  0.163974   1.0  2.0
# RO            region_R23  0.161859   1.0  2.0
# RR            region_R23  0.161859   1.0  2.0
# SE            region_R23  0.163974   1.0  1.0
# TO            region_R23  0.161859   1.0  1.0
# -----  cset2  (same for R23)
# BA            region_R23  0.163974   1.0  2.0
# ES            region_R23  0.204892   1.0  2.0
# MA            region_R23  0.163974   1.0  1.0
# PI            region_R23  0.163974   1.0  1.0
# RJ            region_R23  0.204892   1.0  1.0
# RS            region_R23  0.215327   1.0  2.0
# SC            region_R23  0.215327   1.0  2.0
# SP            region_R23  0.204892   1.0  1.0

