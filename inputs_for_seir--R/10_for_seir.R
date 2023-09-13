# mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
# Script:  ----
# For: create data files for use in SEIR model whith Python
# By: Ricardo Aguirre Leal
# Date: May 2021
# mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries and functions ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(readxl);    library(beepr)
library(foreach);   library(forecast)
library(tidyverse); library(magrittr)
library(lubridate)
# +++++++++++
source('11_funcs_for_seir.R')


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Options ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cj_data  = '26abr2021'        # '12abr2021'  '22mar2021'  
data_ini = '2020-03-01'       # week 10 - begin sunday - total 53 weeks per year
dir_data = '01_data/'         # directory for data files, except for vaccine file
dir_vac  = 'E:/covid_19_vac/' # directory for vaccine file
dir_save_data = '02_output_for_seir/' # directory for save files
# +++++++++++
ini.dir  = getwd()
data_ini = as.Date(data_ini)
# +++++++++++
options(stringsAsFactors = F) 


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Get, check and join datas ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Get df for cases and deaths:
df_covid = getdf_covid()              
# Check regularity conditions in df_covid. 'view_neg' for view lines whith negative values
check_df_covid(view_neg = FALSE)       

ufs = unique(df_covid$uf) %>% sort()  # brasilian states by alphabetic order

# Get df for vaccines:
df_vac = getdf_vac()  # require a lot of PC memmory and some minutes
# Check regularity conditions in df_vac. 'view_neg' for view lines whith negative values
check_df_vac(view_neg = TRUE)   
# Remove vaccinated before 2021 (that also will remove negative entries)
df_vac %<>% filter(year(data) > 2020)

# Same dates for cases, deaths and vaccines:
if(max(df_vac$data) <= max(df_covid$data)){
  df_covid %<>% filter(data <= max(df_vac$data))
} else {
  df_vac %<>% filter(data <= max(df_covid$data))
}
# Join cases, deaths and vaccines  
dados = left_join(df_covid, df_vac) %>% arrange(uf, data)

# +++++++++++ 

# Get df for Google mobility:
df_mobil = getdf_mobil()
# Check regularity conditions in df_mobil. 'view_neg' for view lines whith negative values
check_df_mobil(view_neg = TRUE)

# Same dates for dataframes:
if(max(df_mobil$data) <= max(dados$data)){
  dados %<>% filter(data <= max(df_mobil$data))
} else {
  df_mobil %<>% filter(data <= max(dados$data))
}
# Join dataframes
dados %<>% left_join(df_mobil) %>% arrange(uf, data)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Adjusts and new variables ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dados$population %<>% as.integer()
names(dados)[1:3] = c('region', "iso2", 'date')

dados %<>% mutate(total_vaccinations_per_hundred = total_vaccinations / (population/100) )
dados %<>% mutate(people_vaccinated_per_hundred  = people_vaccinated / (population/100) )
dados %<>% mutate(people_fully_vaccinated_per_hundred = people_fully_vaccinated / (population/100) )

# Turn NA in vaccinations variables equal to zero (in previous period vaccination)
dados$total_vaccinations %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
dados$people_vaccinated  %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
dados$people_fully_vaccinated  %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Adjust data seasonally and by splines ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Seasonally (dates are shortened for only entire weeks)
dados_seas = adj_df_seas()

# Smooth by splines. 's_par' is the used parameter measure of smooth
dados_spline = adj_df_splines(s_par = 0.4)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Datas by region ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# (geographic groups of states and the country)

regioes  = unique(dados$region)     # vector of regions
dias     = unique(dados$date)       # vector of days in 'dados'
dias_saz = unique(dados_seas$date)  # vector of days in 'dados_seas'

# Construct region dataframe
dados_reg = group_reg(dados, dias)

# Construct region dataframe whith seasonally adjust
dados_reg_seas = group_reg(dados_seas, dias_saz)

# Construct region dataframe whith smooth splines
dados_reg_spline = group_reg(dados_spline, dias)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Dataframe for proportion of infected people ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Table of 'p_a' and 'IFR' by states. Copy from PDF - https://doi.org/10.1101/2020.08.18.20177626 (supplementary)
# The 'p_a' "[...] is the percentage of Brazilians that have been infected by SARS-CoV-2 [...]" (p. 4)
# Values by rounds (1, 2 and 3) of the EPICOVID19-BR survey (see the cited article for get the reference)
{
  x = list()
  x[1] = "AC 3.9 (1.7—6.9) 5.7 (3.4—8.6) 4.9 (2.7—7.7) 0.35 (0.32—0.37) 0.53 (0.51—0.54) 0.58 (0.58—0.58) 0.72 (0.4—1.6) 0.81 (0.53—1.4) 1. (0.58—1.9) 0.85 (O.63—l.2)"
  x[2] = "AL 0.53 (0.—2.4) 9.6 (6.2—14.) 13. (9.3—18.) 0.26 (0.24—0.28) 0.43 (0.42—0.45) 0.49 (0.48—0.49) 1.6 (0.48—25.) 0.41 (0.29—0.64) 0.35 (0.23—0.51) 0.47 (0.37—0.61)"
  x[3] = "AM 11. (7.6—16.) 13. (9.4—18.) 7.5 (4.4—11.) 0.64 (0.62—0.66) 0.69 (0.69—0.69) 0.55 (0.53—0.58) 0.51 (0.36—0.77) 0.48 (0.29—0.71) 0.47 (0.28—1.3) 0.5 (0.38—0.65)"
  x[4] = "AP 8.6 (5.1—13.) 14. (9.3—19.) 13. (9.—18.) 0.3 (0.27—0.34) 0.43 (0.41—0.45) 0.49 (0.47—0.5) 0.32 (0.21—0.53) 0.29 (0.2—0.43) 0.34 (0.21—0.52) 0.31 (0.24—0.4)"
  x[5] = "BA 0.88 (0.16—1.8) 3.9 (1.9—6.7) 2.3 (0.91—4.2) 0.12 (0.11—0.13) 0.22 (0.21—0.23) 0.3 (0.3—0.31) 0.99 (0.47—3.9) 0.47 (0.27—0.95) 1. (0.54—2.5) 0.81 (0.57—1.2)"
  x[6] = "CE 6.7 (3.8—10.) 13. (9.2—17.) 18. (13.—22.) 0.74 (0.71—0.78) 0.96 (0.95—0.96) 0.89 (0.87—0.91) 0.98 (0.63—1.7) 0.68 (0.47—1.) 0.52 (0.29—0.7) 0.67 (0.55—0.82)"
  x[7] = "DF 0. (0.—1.5) 0. (0.—2.2) 0. (0.—2.2) 0.068 (0.062—0.075) 0.16 (0.14—0.17) 0.28 (0.26—0.29) 0.69 (0.13—27.) 1.1 (0.25—31.) 1.9 (0.44—54.) 1.3 (0.61—3.5)"
  x[8] = "ES 1.1 (O.—2.8) 1.5 (0.43—2.8) 2.6 (1.3—4.3) 0.18 (0.17—0.2) 0.37 (0.35—0.39) 0.5 (0.49—0.51) 0.92 (0.34—7.) 1.9 (0.95—5.6) 1.6 (0.94—3.2) 1.6 (1.1—2.5)"
  x[9] = "GO 0.3 (0.—1.2) 0.38 (0.—1.2) 0.95 (0.02—2.3) 0.034 (0.032—0.036) 0.068 (0.061—0.076) 0.15 (0.14—0.16) 0.39 (0.13—5.9) 0.82 (0.29—9.9) 0.96 (0.39—5.9) 0.78 (0.47—1.6)"
  x[10] = "MA 1.1 (0.—3.6) 6.7 (4.4—9.4) 7.1 (4.9—9.6) 0.4 (0.39—0.42) 0.53 (0.52—0.53) 0.54 (0.54—0.54) 1.6 (0.54—16.) 0.71 (0.45—1.1) 0.7 (0.43—1.1) 0.82 (0.61—1.1)"
  x[11] = "MG 0.95 (0.23—1.9) 0.46 (0.08—0.92) 0.63 (0.1—1.3) 0.024 (0.023—0.026) 0.049 (0.045—0.053) 0.097 (0.09—0.1) 0.18 (0.09—0.61) 0.77 (0.37—3.2) 1.1 (0.51—4.6) 0.78 (0.52—1.3)"
  x[12] = "MS 0. (0.—2.4) 0.2 (0.—1.4) 0.49 (0.—1.5) 0.0074 (0.0069—0.0079) 0.021 (0.018—0.025) 0.068 (0.052—0.083) 0.04 (0.01—1.2) 0.22 (0.06—4.3) 0.64 (0.23—7.2) 0.4 (0.22—0.89)"
  x[13] = "MT 4.3 (0.12—13.) 0.86 (0.04—2.) 0.61 (0.—1.5) 0.035 (0.029—0.041) 0.15 (0.13—0.17) 0.38 (0.35—0.41) 0.04 (0.01—0.26) 1.1 (0.48—6.3) 3.7 (1.5—26.) 2.3 (1.4—4.4)"
  x[14] = "PA 12. (9.1—15.) 13. (10.—16.) 6.3 (4.5—8.4) 0.7 (0.67—0.74) 0.9 (0.9—0.91) 0.83 (0.81—0.86) 0.56 (0.44—0.75) 0.66 (0.5—0.87) 1.3 (O.63—1.9) 0.66 (0.56—0.78)"
  x[15] = "PB 1.5 (0.05—3.8) 4.8 (2.9—7.3) 1.7 (0.45—3.3) 0.15 (0.14—0.16) 0.28 (0.26—0.29) 0.42 (0.4—0.43) 0.57 (0.22—3.4) 0.51 (0.34—0.86) 1.8 (0.88—5.8) 0.89 (0.64—1.3)"
  x[16] = "PE 2.4 (0.76—4.8) 1.8 (0.34—4.) 0.64 (0.—2.) 0.48 (0.45—0.51) 0.67 (0.66—0.68) 0.71 (0.7—0.71) 1.5 (0.74—4.2) 2.5 (1.1—9.7) 5. (1.5—54.) 2.7 (1.7—4.8)"
  x[17] = "PI 3. (0.4—7.5) 1.3 (0.28—2.7) 8.1 (5.4—11.) 0.12 (0.1—0.13) 0.3 (0.28—0.32) 0.49 (0.47—0.51) 0.24 (0.1—1.1) 1.7 (0.79—6.) 0.56 (0.39—0.84) 0.78 (0.58—1.1)"
  x[18] = "PR 0.53 (0.02—1.2) 0.64 (O.—l.8) 0.79 (0.21—1.5) 0.025 (0.023—0.027) 0.052 (0.048—0.056) 0.1 (0.092—0.11) 0.29 (0.12—1.7) 0.42 (0.15—3.9) 0.97 (0.49—3.1) 0.69 (0.45—1.2)"
  x[19] = "RJ 1.6 (0.11—3.8) 5.5 (2.9—9.) 7.9 (4.8—12.) 0.51 (0.48—0.54) 0.72 (0.71—0.74) 0.79 (0.78—0.79) 2. (0.8—11.) 1.1 (0.67—2.1) 0.87 (0.5—1.5) 1.2 (0.85—1.6)"
  x[20] = "RN 2. (0.36—4.4) 3.2 (1.3—5.7) 4.9 (2.7—7.7) 0.15 (0.14—0.17) 0.33 (0.3—0.36) 0.55 (0.53—0.57) 0.52 (0.23—2.) 0.84 (0.46—2.) 0.98 (0.61—1.8) 0.85 (0.61—1.3)"
  x[21] = "RO 0.05 (0.—2.2) 1.9 (0.32—4.3) 5.3 (2.8—8.6) 0.25 (0.22—0.27) 0.48 (0.47—0.5) 0.6 (0.59—0.61) 1.7 (0.43—40.) 1.7 (0.73—6.8) 0.95 (0.56—1.8) 1.3 (0.85—2.)"
  x[22] = "RR 3.4 (1.2—6.7) 24. (18.—29.) 20. (15.—26.) 0.24 (0.23—0.26) 0.46 (0.43—0.48) 0.67 (0.65—0.68) 0.55 (0.28—1.5) 0.19 (0.15—0.24) 0.31 (0.22—0.42) 0.26 (0.22—0.32)"
  x[23] = "RS 0.43 (0.01—0.93) 0.45 (0.01—0.99) 0.46 (0.04—0.98) 0.025 (0.024—0.026) 0.038 (0.036—0.04) 0.066 (0.062—0.07l) 0.37 (0.16—2.3) 0.54 (0.22—3.3) 0.97 (0.42—5.1) 0.67 (0.43—1.2)"
  x[24] = "SC 0.5 (0.09—0.99) 0.51 (0.1—0.99) 0.59 (0.13—1.1) 0.022 (0.021—0.022) 0.031 (0.029—0.033) 0.054 (0.049—0.058) 0.31 (0.15—1.2) 0.44 (0.21—1.7) 0.68 (0.32—2.4) 0.49 (0.33—0.82)"
  x[25] = "SE 0. (0.—l.6) 0.33 (0.—2.) 2.7 (0.75—5.5) 0.11 (0.1—0.13) 0.29 (0.27—0.32) 0.52 (0.5—0.55) 1.1 (0.27—29.) 2.1 (0.62—39.) 1.4 (0.68—4.4) 1.5 (0.92—3.)"
  x[26] = "SP 2. (0.41—4.3) 1.1 (0.—2.8) 0.69 (0.—2.1) 0.27 (0.26—0.29) 0.36 (0.35—0.36) 0.4 (0.39—0.4) 0.93 (0.42—3.4) 1.8 (0.65—12.) 2.7 (0.88—28.) 1.7 (1.1—3.3)"
  x[27] = "TO 0.35 (0.—1.1) 0.58 (0.—1.5) 1.2 (0.35—2.3) 0.058 (0.053—0.063) 0.12 (0.11—0.12) 0.16 (0.16—0.17) 0.78 (0.29—9.1) 1.1 (0.43—9.2) 0.99 (0.48—3.) 0.97 (O.63—l.7)"
}

pa_ifr = foreach(i=1:27, .combine = rbind) %do% {
  c = x[[i]] %>% strsplit(" ") %>% unlist()
  int = c(substr(c, 1, 1) == "(") %>% which()
  d = c[-c(1, int)] %>% as.numeric()
  c[int[1:3]] = chartr("lO", "10", c[int[1:3]])
  e = chartr("'('—')'-", "        ", c[int[1:3]]) %>% 
        str_trim(side="both")%>% strsplit(" ") %>% unlist() %>% as.numeric()
  f = c( d[1], e[1:2], d[2], e[3:4], d[3], e[5:6], d[4:10] )
  v = matrix(f, nrow=1, dimnames = list(c[1], NULL))
}
rm(x, c, int, d, e, f, v, i)
pa_ifr  = cbind(ufs, as.data.frame(pa_ifr))
colnames(pa_ifr) = c("iso2", 
                     "pa_R1", "pa_inf_R1", "pa_sup_R1", 
                     "pa_R2", "pa_inf_R2", "pa_sup_R2", 
                     "pa_R3", "pa_inf_R3", "pa_sup_R3", 
                     "pd_R1", "pd_R2", "pd_R3", 
                     "IFR_R1", "IFR_R2", "IFR_R3", "IFR_all")
pa_ifr %<>% tibble()  # the dataframe of above table


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Dataframe for proportion of notified cases ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Dataframe for proportion of confirmed cases relative to population, in the rounds of EPICOVID19-BR survey
df_epi = make_df_epi(dados_spline)  # choose 'dados', 'dados_seas' or 'dados_spline'

# Bind dataframes
df_exp = bind_cols(pa_ifr, df_epi[,-1]) %>% arrange(iso2)

# Create variables for proportion of notified cases relative to all that have been infected
# One proportion to the 3th round, and another for mean of 2th and 3th rounds values
df_exp %<>% mutate(
  p_notif_R23 = ( (p_casos_R2/(pa_R2+p_casos_R2)) + (p_casos_R3/(pa_R3+p_casos_R3)) ) / 2,
  p_notif_R3  = p_casos_R3 / (pa_R3+p_casos_R3)
)

#++++++++++++  

# Include values for regions; and include region mean cols for states
df_exp = include_reg_df_exp()


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save datas ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ultimo_dia = max(dados$date)

write_csv(dados, paste0(dir_save_data, 'dados_brasil_', ultimo_dia, '.csv'))
write_csv(dados_seas, paste0(dir_save_data, 'dados_brasil_seas_', ultimo_dia, '.csv'))
write_csv(dados_spline, paste0(dir_save_data, 'dados_brasil_spline_', ultimo_dia, '.csv'))
#++++++++++++  
write_csv(dados_reg, paste0(dir_save_data, 'dados_brasil_reg_', ultimo_dia, '.csv'))
write_csv(dados_reg_seas, paste0(dir_save_data, 'dados_brasil_reg_seas_', ultimo_dia, '.csv'))
write_csv(dados_reg_spline, paste0(dir_save_data, 'dados_brasil_reg_spline_', ultimo_dia, '.csv'))
#++++++++++++  
write.csv(df_exp, paste0(dir_save_data, 'ifr_brasil_epicovid.csv'), row.names = F)

# csv whith a comma for the decimal point and a semicolon for the separator:
write_csv2(dados, paste0(dir_save_data, 'dados_brasil_', ultimo_dia, '(excel).csv'))
write_csv2(dados_seas, paste0(dir_save_data, 'dados_brasil_seas_', ultimo_dia, '(excel).csv'))
write_csv2(dados_spline, paste0(dir_save_data, 'dados_brasil_spline_', ultimo_dia, '(excel).csv'))
#++++++++++++
write_csv2(dados_reg, paste0(dir_save_data, 'dados_brasil_reg_', ultimo_dia, '(excel).csv'))
write_csv2(dados_reg_seas, paste0(dir_save_data, 'dados_brasil_reg_seas_', ultimo_dia, '(excel).csv'))
write_csv2(dados_reg_spline, paste0(dir_save_data, 'dados_brasil_reg_spline_', ultimo_dia, '(excel).csv'))
#++++++++++++  
write.csv2(df_exp, paste0(dir_save_data, 'ifr_brasil_epicovid(excel).csv'), row.names = F)



