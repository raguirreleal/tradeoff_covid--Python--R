# mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
# Script:  ----
# For: relate mobility and contagious rate
# By: Ricardo Aguirre Leal
# Date: May 2021
# mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries and functions ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(readxl);    library(foreach)
library(tsDyn);     library(urca)
library(lmtest);    library(strucchange)
library(vars);      library(forecast)
library(tidyverse); library(magrittr)
library(doParallel)
# +++++++++++
source('21_funcs_mob_gamma.R')


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Options ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cj_pdreg      = 'output_pdreg_2021-04-21'    # set (directory) of 'df3' data files whith 'pd_reg'
cj_pdstates   = 'output_pdstates_2021-04-21' # set (directory) of 'df3' data files whith 'pd_states'
dir_data      = '03_df_from_seir/'           # directory of 'cj_data'
dir_save_data = '04_output_mob_gamma/'       # directory for save files
beg_time      = '2020-04-01'  # NULL
end_time      = NULL  # '2020-12-31'  # NULL
# +++++++++++
ini.dir  = getwd()
# +++++++++++
options(stringsAsFactors = F) 


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get, check and adjust datas ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Get list of dfs for 'cj_pdreg' (date, gamma_t, mobility_L, mobility_W and mobility_F):
ldf_pdreg = getlist_df3s(cj_pdreg)
# Get list of dfs for 'cj_pdstates':
ldf_pdstt = getlist_df3s(cj_pdstates)

# Initial dates for each state/region (in seir there is a 'virus_thres' param for epidemic starting point)
# Change to TRUE for print
if(F) lapply(ldf_pdreg, function(L) L$date[1]) %>% as.data.frame() %>% t() # is the same for 'cj_pdstates'

# Eliminate first and last 'el_x' observations (days) in each df 
# (gamma_t in seir can suffer distortions in theses)
# And ensure that observations begin not before 'beg_time' (if beg_time is NULL, then there is no cutoff)
# Also ensure the same for 'end_time'
el_x = 5 # 7
ldf_pdreg %<>% elim(el_x) 
ldf_pdstt %<>% elim(el_x) 

# Print new initial dates - Change to TRUE for print
if(F) lapply(ldf_pdreg, function(L) L$date[1]) %>% as.data.frame() %>% t() # is the same for 'cj_pdstates'

# Vectors of units
regioes = c("Brasil", "Centro-Oeste", "Nordeste", "Norte", "Sudeste", "Sul")
ufs     = names(ldf_pdreg)[nchar(names(ldf_pdreg))==2] %>% sort
# Quantity of units
nunits  = length(ldf_pdreg)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# List whith time series first difference or nivel ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# List whith time series for 'pdreg'
lts_pdreg  = makelist_ts(ldf_pdreg, dif = FALSE) 
# List whith time series for 'pdstt'
lts_pdstt  = makelist_ts(ldf_pdstt, dif = FALSE) 
# ...in first difference for 'pdreg'
ltsd_pdreg = makelist_ts(ldf_pdreg, dif = TRUE) 
# ...in first difference for 'pdstt'
ltsd_pdstt = makelist_ts(ldf_pdstt, dif = TRUE) 

# Plot mobility_W and gamma_t for regions
if(FALSE){   # run inside or change to TRUE
  plot(lts_pdreg$Brasil$mobility_W); for (i in regioes[-1]) lines(lts_pdreg[[i]]$mobility_W)
  plot(lts_pdreg$Brasil$gamma_t, ylim = c(0, 3)); for (i in regioes[-1]) lines(lts_pdreg[[i]]$gamma_t)
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Check unit root in time series ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++
# +++ In nivel:
# ADF or ERS unit root test (augmented Dickey-Fuller) (Elliott, Rothenberg & Stock)
# Matrix whith indicator of unit root (T or F) for each serie
ur0_pdreg = testlist_ur(lts_pdreg, 'ers')  # select 'adf' or 'ers'
ur0_pdstt = testlist_ur(lts_pdstt, 'ers')  # select 'adf' or 'ers'

# +++++++++++
# +++ In first difference:
# ADF or ERS test
ur1_pdreg = testlist_ur(ltsd_pdreg, 'ers')  # select 'adf' or 'ers'
ur1_pdstt = testlist_ur(ltsd_pdstt, 'ers')  # select 'adf' or 'ers'
# Zivot-Andrews test (considering structural change)
urza1_pdreg  = testlist_urza(ltsd_pdreg)
urza1_pdstt  = testlist_urza(ltsd_pdstt)
# +++ 
# Obs: Structural breaks in many series. All units have breaks in some serie.
#      Alowing for a unique structural break in series, all not have unit root in first difference
# +++

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VAR models in first difference and diagnostics ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Select information criterium
criterio_var = "SC"  # "AIC", "HQ", "SC", "FPE"
# Select mobility variables for VAR: 
# 'mobility_L', 'mobility_W', 'mobility_F', 'mobility_P' , 'mobility_T', 'mobility_R'
vars_mobil = c('mobility_L', 'mobility_W', 'mobility_F', 'mobility_P', 'mobility_R')

lvar_pdreg = makelist_var(ltsd_pdreg, vars_mobil)
lvar_pdstt = makelist_var(ltsd_pdstt, vars_mobil)

# Test for structural change -- type 'F' (Fstats) or 'efp' (see 'strucchange' package)
# List of p-value (H0: no structural change) for each resid in VARs
ltscVAR_pdreg = makelist_sc(lvar_pdreg, test = "efp", print = T) 
ltscVAR_pdstt = makelist_sc(lvar_pdstt, test = "efp", print = T) 

# Plot stability test -- method from the generalised fluctuation test framework
if(FALSE) {  # run inside or change to TRUE
  for(i in 1:nunits) stability(lvar_pdreg[[i]], type = "Rec-CUSUM") %>%  plot()
  for(i in 1:nunits) stability(lvar_pdreg[[i]], type = "OLS-CUSUM") %>%  plot()
}

# +++ 
# Obs: Possible whith structural breaks in models
# +++

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Test for cointegration of series ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Select mobility variables for VAR: 
vars_mobil = c('mobility_L', 'mobility_W', 'mobility_F', 'mobility_P', 'mobility_R')
# Get nlags in VAR nivel for cointegration test
lnlags_pdreg = get_nlags(lts_pdreg, vars_mobil) 
lnlags_pdstt = get_nlags(lts_pdstt, vars_mobil)

# Johansen Procedure for VAR, with constant:
# H0:r<=r0; H1:r>r0
# Os testes são realizados em sequência, de forma crescente, até que a hipótese nula não seja rejeitada
# If r=0, no cointegration
# If r=k (where k is #variables), then full rank and stationary variables (VAR!, nor VEC)
lcajo_pdreg = makelist_cajo(lts_pdreg, vars_mobil, lnlags_pdreg, print = F)
lcajo_pdstt = makelist_cajo(lts_pdstt, vars_mobil, lnlags_pdstt, print = F)

# Results of rank selection by ca.jo
lrcajo_pdreg = makelist_rcajo(lcajo_pdreg, pct=5, print = T)
lrcajo_pdstt = makelist_rcajo(lcajo_pdstt, pct=5, print = T)

# Test if is need include trend in VECs
llttest_pdreg = makelist_lttest(lcajo_pdreg, lrcajo_pdreg, print = T)
llttest_pdstt = makelist_lttest(lcajo_pdstt, lrcajo_pdstt, print = T)

# +++ 
# Obs: The models are not VAR, but VEC (whitout include trend)
# +++ 

# +++++

# Test for structural change -- type 'F' (Fstats) or 'efp' (see 'strucchange' package)
# List of p-value (H0: no structural change) for each resid in VECs from 'ca.jo'
ltscVEC_pdreg = makelist_sc(lcajo_pdreg, lrcajo_pdreg, test = "efp", print = T) 
ltscVEC_pdstt = makelist_sc(lcajo_pdstt, lrcajo_pdstt, test = "efp", print = T) 

# Plot 'strucchange' test for series whith p-value below 0.05
# The 'test' must be the same of that create ltscVEC
plot_ltscVEC(ltscVEC_pdreg)

# +++ 
# Obs: Most of breaks are in around early november 
# +++ 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VEC models and diagnostics ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The models VEC which structural change (for those have)
scVEC_pdreg = makelist_scVEC(ltscVEC_pdreg)
scVEC_pdstt = makelist_scVEC(ltscVEC_pdstt)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# IFR and VD for VAR models ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

scVEC_pdreg$lVEC2VAR$Brasil %>% irf(response = 'gamma_t', n.ahead = 90) %>% plot() 
scVEC_pdreg$lVEC2VAR$Brasil %>% irf(response = 'gamma_t', n.ahead = 60, cumulative = T) %>% plot() 
scVEC_pdreg$lVEC2VAR$Brasil %>% fevd(n.ahead=10) %>% plot()



scVEC_pdreg$lVEC2VAR$Brasil %>% logLik() %>% BIC() # com T: -36035.75; sem T: 




