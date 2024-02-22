# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 09:27:49 2020
@author: Phurichai Rungcharoenkitkul
github.com/phurichai/covid19macro
-- (very) Modified --
on: Apr 2021
by: Ricardo Aguirre Leal
"""

from scipy.optimize import dual_annealing
# from datetime import date, datetime
from datetime import timedelta
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from scipy.linalg import lstsq
# from pathlib import Path
import pandas as pd
import numpy as np
# import pickle  
import time
import math
# ------------------
from param import *
# ------------------
np.warnings.filterwarnings('ignore') 

# ============================
# Calc functions  
# ============================

def calc_step(x, pars):
    """
    SEIR model building on DELPHI v.3
    Features 16 distinct states, taking into account undetected, 
    deaths, hospitalized and recovered
    [0 S, 1 E, 2 I, 3 UR, 4 DHR, 5 DQR, 6 UD, 7 DHD, 8 DQD, 9 R, 10 D, 
     11 TH, 12 DVR,13 DVD, 14 DD, 15 DT, 16 V]
    """    
    S, E, I, AR, DHR, DQR, AD, DHD, DQD, R, D, TH, DVR, DVD, DD, DT, V = x
    gamma_t, p_dth, r_d, r_dth, r_v, r_re, N, p_d_i, p_h_i = pars   
    # ------------------
    S1 = S - gamma_t * S * I / N + r_re * R + r_re * V - r_v
    if S1 < 0: # Vaccination reaches saturating point
        S1 = 0
        r_v = S - gamma_t * S * I / N + r_re * R + r_re * V
    E1 = E + gamma_t * S * I / N - r_i * E
    I1 = I + r_i * E - r_d * I
    AR1 = AR + r_d * (1 - p_dth) * (1 - p_d_i) * I - r_ri * AR
    DHR1 = DHR + r_d * (1 - p_dth) * p_d_i * p_h_i * I - r_rh * DHR
    DQR1 = DQR + r_d * (1 - p_dth) * p_d_i * (1 - p_h_i) * I - r_ri * DQR
    AD1 = AD + r_d * p_dth * (1 - p_d_i) * I - r_dth * AD
    DHD1 = DHD + r_d * p_dth * p_d_i * p_h_i * I - r_dth * DHD
    DQD1 = DQD + r_d * p_dth * p_d_i * (1 - p_h_i) * I - r_dth * DQD
    R1 = R + r_ri * (AR + DQR) + r_rh * DHR - r_re * R
    D1 = D + r_dth * (AD + DQD + DHD)
    # Helper states 
    TH1 = TH + r_d * p_d_i * p_h_i * I
    DVR1 = DVR + r_d * (1 - p_dth) * p_d_i * p_h_i * p_v * I - r_rv * DVR
    DVD1 = DVD + r_d * p_dth * p_d_i * p_h_i * p_v * I - r_dth * DVD
    DD1 = DD + r_dth * (DHD + DQD)
    DT1 = DT + r_d * p_d_i * I
    V1 = V + r_v - r_re * V
    x1 = [S1, E1, I1, AR1, DHR1, DQR1, AD1, DHD1, DQD1,
          R1, D1, TH1, DVR1, DVD1, DD1, DT1, V1]
    return x1

def calc_initial_states_func(params1, params2):
    N, PopulationCI, PopulationR, PopulationD, PopulationI, p_d_i, p_h_i, p_v = params1
    k, r_dth, p_d_i, p_h_i, p_dth0 = params2
    # ------------------
    E_0 = PopulationCI / p_d_i * k
    I_0 = PopulationCI / p_d_i * k
    UR_0 = (PopulationCI / p_d_i - PopulationCI) * (1 - p_dth0)
    DHR_0 = (PopulationCI * p_h_i) * (1 - p_dth0)
    DQR_0 = PopulationCI * (1 - p_h_i) * (1 - p_dth0)
    UD_0 = (PopulationCI / p_d_i - PopulationCI) * p_dth0
    DHD_0 = PopulationCI * p_h_i * p_dth0
    DQD_0 = PopulationCI * (1 - p_h_i) * p_dth0
    R_0 = PopulationR / p_d_i
    D_0 = PopulationD / p_d_i
    S_0 = N - (E_0 +I_0 +UR_0 +DHR_0 +DQR_0 +UD_0 +DHD_0 +DQD_0 +R_0 +D_0)
    TH_0 = PopulationCI * p_h_i
    DVR_0 = (PopulationCI * p_h_i * p_v) * (1 - p_dth0)
    DVD_0 = (PopulationCI * p_h_i * p_v) * p_dth0
    DD_0 = PopulationD
    DT_0 = PopulationI
    V_0 = 0
    x_init = [S_0, E_0, I_0, UR_0, DHR_0, DQR_0, UD_0, DHD_0, DQD_0, 
              R_0, D_0, TH_0, DVR_0, DVD_0, DD_0, DT_0, V_0] 
    return x_init

# ============================
# Main class 'solveCovid'   
# ============================

class solveCovid:  
    def __init__(self,iso2: str, r_dth, r_d, p_d):  
        self.iso2 = iso2
        self.r_dth = r_dth
        self.r_d = r_d
        self.p_d = p_d
        self.phi_min = 1e-13  # Lowerbound for phi - authorities care about output
        self.pdth_min = 0.005 # Lowerbound on death probability - society still think there is death probability
        self.gamma_max = 15 
        self.effi_one = effi_one # Efficacy after one dose
        self.effi_two = effi_two # Efficacy after two doses
        self.r_re = 2**(1/days_lost_immunity)-1 # Baseline: lost immunity after x days

    def prelim(self, lreg):
        iso2 = self.iso2
        if lreg:
            self.N = dfr1.fillna(method='ffill')['population'][iso2].iloc[-1]
            df2 = dfr1.iloc[:,dfr1.columns.get_level_values(1)==iso2][[
                    'total_cases','total_deaths','new_cases','new_deaths',
                    'mobility_L', 'mobility_W', 'mobility_F',
                    'mobility_P', 'mobility_T', 'mobility_R', 'vac_total',
                    'vac_partial','vac_fully']][dfr1['total_cases'][iso2] > virus_thres] 
        else:
            self.N = df1.fillna(method='ffill')['population'][iso2].iloc[-1]
            df2 = df1.iloc[:,df1.columns.get_level_values(1)==iso2][[
                    'total_cases','total_deaths','new_cases','new_deaths',
                    'mobility_L', 'mobility_W', 'mobility_F',
                    'mobility_P', 'mobility_T', 'mobility_R', 'vac_total',
                    'vac_partial','vac_fully']][df1['total_cases'][iso2] > virus_thres] 
        df2['vac_total']   = df2['vac_total'].interpolate()
        df2['vac_partial'] = df2['vac_partial'].interpolate()
        df2['vac_fully']   = df2['vac_fully'].interpolate()
        if np.isnan(df2['vac_partial'].iloc[-1].values[0]) and ~np.isnan(df2['vac_total'].iloc[-1].values[0]):
            df2['vac_partial'] = 0.8 * df2['vac_total'] # If no data on breakdowns exist, do manual approximation
            df2['vac_fully'] = 0.2 * df2['vac_total']
        df2 = df2.fillna(0).droplevel('iso2',axis=1) 
        PopulationI = df2['total_cases'][0]
        PopulationD = df2['total_deaths'][0]
        if PopulationD==0:
            PopulationD = 0
            PopulationR = 5
        else:
            PopulationR = PopulationD * 5
        PopulationCI = PopulationI - PopulationD - PopulationR # Undetected and infectious cases
        self.cases_data_fit = df2['total_cases'].tolist()
        self.deaths_data_fit = df2['total_deaths'].tolist()        
        self.newcases_data_fit = df2['new_cases'].tolist()
        self.newdeaths_data_fit = df2['new_deaths'].tolist() 
        self.maxT = len(df2['total_cases'])
        self.T = len(df2)
        self.t_cases = np.arange(0,self.T)
        self.GLOBAL_PARAMS = [self.N, PopulationCI, PopulationR, PopulationD, PopulationI, self.p_d[iso2], p_h[iso2], p_v]
        self.gamma_0_days = 10  # average of gamma_t during first n days becomes the target 
        self.vac_partial = df2['vac_partial'].values
        self.vac_fully = df2['vac_fully'].values     
        df2['V_'] = self.N * (self.effi_one*df2['vac_partial']+self.effi_two*df2['vac_fully'])/100 # V = expected number of effectively vaccinated persons
        default_maxT1 = default_maxT + timedelta(days=1)
        ix = pd.date_range(start=df2.index[0], end=default_maxT1, freq='D') # Expand time-sample, to include forecast later
        df_v = df2.reindex(ix)
        df_v['V_'] = df_v['V_'].interpolate().clip(0,self.N)
        self.df2 = df2
        self.df_v = df_v
        self.p_dth0 = self.newdeaths_data_fit[0] / (self.r_dth * self.GLOBAL_PARAMS[1]) # Set p_dth0 to match D1-D0 to newdeaths_data_fit        
        
    def step_seir(self, t, x, gamma_t, p_dth) -> list:    
        r_v   = self.df_v['V_'].iloc[t+1] - self.df_v['V_'].iloc[t]
        pars = [gamma_t, p_dth, self.r_d, self.r_dth, r_v, self.r_re, self.N, self.p_d[self.iso2], p_h[self.iso2]]
        x1 = calc_step(np.asarray(x), np.asarray(pars))
        return x1
    
    def loss_gamma0(self,k):
        r_d   = self.r_d
        r_dth = self.r_dth
        p_d_i   = self.p_d[self.iso2]
        newcases  = np.asarray(self.newcases_data_fit)
        newdeaths = np.asarray(self.newdeaths_data_fit)
        gamma_t_vec = np.repeat(float('inf'), self.gamma_0_days)  
        x_init = calc_initial_states_func(np.asarray(self.GLOBAL_PARAMS), np.asarray([k, r_dth, p_d_i, p_h[self.iso2], self.p_dth0]))
        S_0, E_0, I_0, UR_0, DHR_0, DQR_0, UD_0, DHD_0, DQD_0, R_0, D_0, TH_0, DVR_0, DVD_0, DD_0, DT_0, V_0 = x_init
        newcases_sm2 = np.append(newcases, newcases[-2:]) # Extend the list for forward projection below
        newdeaths_sm2 = np.append(newdeaths, newdeaths[-1])
        x_0 = x_init.copy()
        for t in range(self.gamma_0_days): # Target first n days
            gamma_t = np.clip((newcases_sm2[t+2]/(r_d*p_d_i) - (1-r_d)**2 *I_0 - r_i*(2-r_d-r_i)*E_0 )*self.N/(r_i*S_0*I_0), 0.01, self.gamma_max)
            p_dth = np.clip((newdeaths_sm2[t+1] - r_dth*(1-r_dth)*(DHD_0 + DQD_0))/(r_dth*r_d*p_d_i*I_0), 0, 1)
            x_1 = self.step_seir(t, x_0, gamma_t, p_dth)
            x_0 = x_1
            gamma_t_vec[t] = gamma_t      
        loss = (np.mean(gamma_t_vec) - (r_d*6) )**2  # gamma_0 equivalent to R0=6 is 2.08
        return loss  
    
    def get_initial_conditions(self):
        output = dual_annealing(self.loss_gamma0, x0=[5], bounds=[(1,50)], seed=1234)   
        self.kstar = output.x[0]  # kstar that matches gamma_0 to target
        x_init = calc_initial_states_func(np.asarray(self.GLOBAL_PARAMS), np.asarray([self.kstar, self.r_dth, self.p_d[self.iso2], p_h[self.iso2], self.p_dth0]))
        return x_init

    def compute(self):
        r_d   = self.r_d
        r_dth = self.r_dth
        p_d_i = self.p_d[self.iso2]
        newcases  = np.asarray(self.newcases_data_fit)
        newdeaths = np.asarray(self.newdeaths_data_fit)
        lenew = len(newcases)
        # ------------------
        default_maxT1 = default_maxT + timedelta(days=1)
        ix = pd.date_range(start=self.df2.index[0], end=default_maxT1, freq='D') 
        df3 = self.df2.reindex(ix)
        col_temp = ['S', 'E', 'I', 'AR', 'DHR', 'DQR', 'AD', 'DHD', 'DQD', 'R', 'D', 'TH', 'DVR', 'DVD', 'DD', 'DT', 'V']
        df4 = pd.DataFrame(columns=col_temp, index=df3.index)
        # ------------------
        self.x_init = self.get_initial_conditions() 
        S_0, E_0, I_0, AR_0, DHR_0, DQR_0, AD_0, DHD_0, DQD_0, R_0, D_0, TH_0, DVR_0, DVD_0, DD_0, DT_0, V_0 = self.x_init
        gamma_t_vec = np.repeat(float('inf'), lenew)    
        p_dth_vec = np.repeat(float('inf'), lenew)
        newcases_sm2 = np.append(newcases,  newcases[-2:]) # Extend the list for forward projection below
        newdeaths_sm2 = np.append(newdeaths, newdeaths[-1])
        x_0 = self.x_init.copy()
        x_data = self.x_init.copy()
        for t in range(lenew):
            gamma_t = np.clip((newcases_sm2[t+2]/(r_d*p_d_i) - (1-r_d)**2 *I_0 - r_i*(2-r_d-r_i)*E_0 )*self.N/(r_i*S_0*I_0), 0.01, self.gamma_max)
            p_dth = np.clip((newdeaths_sm2[t+1] - r_dth*(1-r_dth)*(DHD_0 + DQD_0))/(r_dth*r_d*p_d_i*I_0), 0, 1)
            x_1 = self.step_seir(t, x_0, gamma_t, p_dth)
            S_0, E_0, I_0, AR_0, DHR_0, DQR_0, AD_0, DHD_0, DQD_0, R_0, D_0, TH_0, DVR_0, DVD_0, DD_0, DT_0, V_0 = x_1
            x_0 = x_1
            gamma_t_vec[t] = gamma_t
            p_dth_vec[t]   = p_dth
            df4.iloc[t]    = x_1
        df3 = df3.merge(df4, how='left', left_index=True, right_index=True)
        df3 = df3.iloc[:(len(df3)-1), :]
        df3['gamma_t'] = gamma_t_vec
        df3['pdth_t']  = p_dth_vec
        self.df3 = df3
        return df3

    def fitmob(self):
        r_d = self.r_d  
        mL_t = self.df3['mobility_L'].values
        mW_t = self.df3['mobility_W'].values
        mF_t = self.df3['mobility_F'].values
        independent = np.column_stack([mL_t**0, mL_t, mL_t**2, mW_t, mW_t**2, mF_t, mF_t**2])
        dependent = self.df3['gamma_t'].values[:len(mL_t)]
        best_betas, resid_reg, rnk, s = lstsq(independent, dependent)   
        beta0, beta1L, beta2L, beta1W, beta2W, beta1F, beta2F = best_betas
        self.best_betas = best_betas
        # ------------------
        self.df3['gammaL_mob'] = beta1L*mL_t + beta2L*(mL_t**2) 
        self.df3['gammaW_mob'] = beta1W*mW_t + beta2W*(mW_t**2) 
        self.df3['gammaF_mob'] = beta1F*mF_t + beta2F*(mF_t**2) 
        self.df3['gamma_mob'] = beta0 + self.df3['gammaL_mob'] + self.df3['gammaW_mob'] + self.df3['gammaF_mob']
        self.df3['gamma_tilde'] = self.df3['gamma_t'] - self.df3['gamma_mob']
        # ------------------
        mL = self.df3['mobility_L'][-1]
        mW = self.df3['mobility_W'][-1]
        mF = self.df3['mobility_F'][-1]
        s  = self.df3['S'][-1]/self.N
        i  = self.df3['I'][-1]/self.N    
        gamma_tilde = self.df3['gamma_tilde'][-1]
        pdth = self.df3['pdth_t'][-1]
        pdth = max(pdth, self.pdth_min) # Get around cases where pdth=0 for countries with very few cases
        # ------------------
        gamma_m = beta0 + beta1L*mL + beta2L*(mL**2) + beta1W*mW + beta2W*(mW**2) + beta1F*mF + beta2F*(mF**2)
        LHS2 = pdth*r_d*i*(1-r_d + (gamma_m + gamma_tilde)*s) # Eq. 2.23
        LHS1L = pdth*r_d*i*s*(beta1L + 2*beta2L*mL) # derivada na eq. 2.24           
        phiL = -(LHS1L * LHS2)/mL  # Eq. 2.24
        LHS1W = pdth*r_d*i*s*(beta1W + 2*beta2W*mW) # derivada na eq. 2.24           
        phiW = -(LHS1W * LHS2)/mW  # Eq. 2.24
        LHS1F = pdth*r_d*i*s*(beta1F + 2*beta2F*mF) # derivada na eq. 2.24            
        phiF = -(LHS1F * LHS2)/mF  # Eq. 2.24
        self.phiL = max(phiL, self.phi_min)
        self.phiW = max(phiW, self.phi_min)
        self.phiF = max(phiF, self.phi_min)
        return best_betas        
    
    def graphics(self, mob):
        df = self.df3
        transpa = 0.0
        if mob:
            fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(15,8), constrained_layout=True)
        else:
            fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,8), constrained_layout=True)
        # ------------------
        ax[0,0].plot(df.index, 100*df['total_cases']/self.N, linewidth = 1, label='Obs Total Cases', color='blue')
        ax[0,0].plot(df.index, 100*df['DT']/self.N, linewidth = 1, label='$DT_t$', color='red')
        ax[0,0].set_title('Cases',fontsize='x-large')
        ax[0,0].set(ylabel = '% of population')
        ax2 = ax[0,0].twinx()
        ax2.plot(df.index, 100*df['I']/self.N, linewidth = 1, label='$I_t$ (rhs)',color='green',linestyle='--')
        lines, labels = ax[0,0].get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='best', framealpha=transpa,fontsize='x-large')
        # ------------------
        ax[1,0].plot(df.index, 100*df['total_deaths']/self.N, linewidth = 1, label='Obs Total Deaths', color='blue')
        ax[1,0].plot(df.index, 100*df['DD']/self.N, linewidth = 1, label='$DD_t$', color='red')
        ax[1,0].set_title('Deaths',fontsize='x-large')
        ax[1,0].set(ylabel='% of population')
        ax4 = ax[1,0].twinx()
        ax4.plot(df.index, df['new_deaths'], linewidth = 1, label='Obs New Deaths',color='orange',linestyle='--')
        lines3, labels3 = ax[1,0].get_legend_handles_labels()
        lines4, labels4 = ax4.get_legend_handles_labels()
        ax4.legend(lines4 + lines3, labels4 + labels3, loc='best', framealpha=transpa,fontsize='x-large')
        # ------------------
        ax[0,1].plot(df.index, df['gamma_t'], linewidth = 1, label=r'$\gamma_t$',color='red')
        if mob:
            ax[0,1].plot(df.index, df['gamma_mob'], linewidth = 1, label=r'$\gamma^{m}_t$', color ='black') 
            ax[0,1].plot(df.index, df['gamma_tilde'], linewidth = 1, label=r'$\gamma^{d}$', color='orange')
        ax[0,1].set_title('Infection rate',fontsize='x-large')
        ax[0,1].legend(loc='best',framealpha=transpa ,fontsize='x-large')
        ax[0,1].axhline(1, linewidth = 2, color='gray', linestyle=':')
        # ------------------        
        ax[1,1].plot(df.index, 100*df['pdth_t'], linewidth = 1, label='Death probability', color='black')
        ax[1,1].legend(loc=0,framealpha=transpa ,fontsize='x-large')
        ax[1,1].set_title('Death probability',fontsize='x-large')
        ax[1,1].set(ylabel='%')
        # ------------------
        if mob:
            ax[1,2].plot(df.index, 100*df['mobility_L'], linewidth = 1, label='Leisure', color='blue')
            ax[1,2].plot(df.index, 100*df['mobility_W'], linewidth = 1, label='Work', color='green')
            ax[1,2].plot(df.index, 100*df['mobility_F'], linewidth = 1, label='Pharmacy', color='orange')
            ax[1,2].legend(loc=0,framealpha=transpa ,fontsize='x-large')
            ax[1,2].set_title('Mobility',fontsize='x-large')
            ax[1,2].set(ylabel='% deviations from norm')
            ax[1,2].axhline(0, linewidth = 2, color='gray', linestyle=':')
            # ------------------
            ax[0,2].plot(df.index, df['gammaL_mob'], linewidth = 1, label=r'$\gamma^{L}_t$', color ='blue')
            ax[0,2].plot(df.index, df['gammaW_mob'], linewidth = 1, label=r'$\gamma^{W}_t$', color ='green')
            ax[0,2].plot(df.index, df['gammaF_mob'], linewidth = 1, label=r'$\gamma^{F}_t$', color ='orange')
            ax[0,2].set_title('Infection rate by mobilities.',fontsize='x-large')
            ax[0,2].legend(loc='best',framealpha=transpa ,fontsize='x-large')
            ax[0,2].axhline(0, linewidth = 2, color='gray', linestyle=':')
        # ------------------
        plt.setp(ax[0,0].get_xticklabels(), rotation=30, horizontalalignment='right')
        plt.setp(ax[0,1].get_xticklabels(), rotation=30, horizontalalignment='right')        
        plt.setp(ax[1,0].get_xticklabels(), rotation=30, horizontalalignment='right')
        plt.setp(ax[1,1].get_xticklabels(), rotation=30, horizontalalignment='right')
        if mob:
            plt.setp(ax[0,2].get_xticklabels(), rotation=30, horizontalalignment='right')
            plt.setp(ax[1,2].get_xticklabels(), rotation=30, horizontalalignment='right')
        # ------------------
        fig.suptitle(f'Graphic for {self.iso2}',fontsize='xx-large')


        




