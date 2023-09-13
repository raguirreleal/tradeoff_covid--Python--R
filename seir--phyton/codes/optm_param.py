# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 16:15:14 2021
@author: Ricardo
"""

import time
import winsound
import pandas as pd 
import numpy  as np
import os
# ------------------
import data                                                                                   
from seir import *
# ------------------
tic = time.perf_counter()
np.warnings.filterwarnings('ignore')

# =================================================
#   Choosing and setting run
# =================================================

# Choose the sets for run and the number of cores for parallel processment:
run_region = True
run_cset0  = True
run_cset1  = True
run_cset2  = True
ncores     = 4
# Save and rewrite new parameters?
save       = False
# ------------------
csetreg = ['Norte', 'Nordeste', 'Centro-Oeste', 'Sudeste', 'Sul', 'Brasil' ] # 6
cset0   = ['AC', 'AL', 'AM', 'AP', 'CE', 'DF', 'GO', 'MG', 'MS', 'MT', 'PB'] # 11
cset1   = ['PA', 'PE', 'PR', 'RN', 'RO', 'RR', 'SE', 'TO'] # 8
cset2   = ['BA', 'ES', 'MA', 'PI', 'RJ', 'RS', 'SC', 'SP'] # 8

# =================================================
#    Deffining functions 
# =================================================

def for_base(c, r_dth, r_d, p_d, lreg):
    tmp = solveCovid(c, r_dth, r_d, p_d)
    tmp.prelim(lreg)
    tmp.compute()        
    return tmp.df3

def run_baseline(cset, r_dth, r_d, p_d, lreg):
    results_base = Parallel(n_jobs=ncores)(delayed(for_base)(c, r_dth, r_d, p_d, lreg) for c in cset)
    MAE_DD = np.repeat(0, len(cset)).astype(float)
    MAE_DT = np.repeat(0, len(cset)).astype(float)
    for i in np.arange(len(cset)):  
        df_temp = results_base[i].dropna()
        MAE_DD[i] = ((df_temp['DD'] - df_temp['total_deaths'])/df_temp['total_deaths']).abs().mean()*100
        MAE_DT[i] = ((df_temp['DT'] - df_temp['total_cases' ])/df_temp['total_cases' ]).abs().mean()*100
    return [MAE_DD, MAE_DT]

def optm_par(cset):
    lreg = False  # logical for region/state run
    if (cset == csetreg):
        lreg = True
    days_r_dth = np.array([1])      # 46
    days_r_d   = np.arange(1, 17)   # 16
    inter =[pd_set.columns.values.astype(str), days_r_dth, days_r_d]
    idx = pd.MultiIndex.from_product(inter, names=["p_d" ,"r_dth", "r_d"])
    df_DD = pd.DataFrame(columns=idx, index=cset).astype(float)
    df_DT = pd.DataFrame(columns=idx, index=cset).astype(float)
    N_c = len(days_r_dth)*len(days_r_d)*len(pd_set.columns)
    n_c = 0
    for z in pd_set.columns:
        p_d = pd_set[z]
        for i in days_r_dth:
            r_dth = 2**(1/i)-1  
            for j in days_r_d:
                r_d = 2**(1/j)-1    
                ret = run_baseline(cset, r_dth, r_d, p_d, lreg)
                df_DD[z, i, j] = ret[0]
                df_DT[z, i, j] = ret[1]
                n_c = n_c + 1
                tocm = time.perf_counter()
                print(f'{n_c}of{N_c} || p_d: {z},  dr_dth: {i},  dr_d: {j} || Total Time: {(tocm - tic)/60:0.2f} min')
    df_DD2 = df_DD.T
    df_DT2 = df_DT.T
    df_DDT2 = (df_DD2 + df_DT2).round(3)   # 3 
    # ------------------    
    df_r = pd.DataFrame(columns=cset, index=["p_d_name" ,"p_d" ,"r_dth", "r_d"])
    for c in cset:
        idx_rs = df_DDT2[c].idxmin()
        df_r.at['p_d_name', c] = idx_rs[0]        
        df_r.at['p_d', c] = pd_set[idx_rs[0]][c].astype(float)
        df_r.at['r_dth', c] = idx_rs[1].astype(float)
        df_r.at['r_d', c] = idx_rs[2].astype(float)     
    df_r = df_r.T
    df_r.index= df_r.index.rename('iso2')
    return [df_r, df_DD2, df_DT2, df_DDT2]

# =================================================
#   Run optimization
# =================================================

l_cset = dict(csetreg=csetreg, cset0=cset0, cset1=cset1, cset2=cset2)
whatcset = np.where([run_region, run_cset0, run_cset1, run_cset2])[0]
ln_cset = np.array(list(l_cset))[whatcset]

for iset in ln_cset:
    ticm = time.perf_counter()
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(f'Initiation of optimization for {iset} in', current_time)
    print(' ')  
    lrs = optm_par(l_cset[iset]) # For 675 combinations, about 1. 4 hours in a Intel i5-6600K @4.6 GHz 4 cores
    # ------------------           
    winsound.Beep(2500, 400)
    toc = time.perf_counter()
    print(' ')
    print(f'Optmization for {iset} done || Total Time: {(toc - ticm)/60:0.2f} minutes')
    print(lrs[0])
    print(' ')
    # ------------------
    if save:
        lrs[0].to_csv(f'{optm_param_folder}\\df_r_{iset}.csv')
        lrs[1].to_csv(f'{optm_param_folder}\\df_DD2_{iset}.csv')
        lrs[2].to_csv(f'{optm_param_folder}\\df_DT2_{iset}.csv')  
        lrs[3].to_csv(f'{optm_param_folder}\\df_DDT2_{iset}.csv') 
        print(f'Parameters saved in {optm_param_folder}')
        print(' ')
    
# =================================================
#   Alert finalization and total time
# =================================================

winsound.Beep(2500, 300)
winsound.Beep(2500, 300)
winsound.Beep(2500, 300)
toc = time.perf_counter()
print(' ')
print(f'All optmizations done || Total Time: {(toc - tic)/60:0.2f} minutes')
print(' ')

