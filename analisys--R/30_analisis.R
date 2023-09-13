

# +++++++++++
library(readxl)
library(foreach)
library(stringr)
library(tidyverse); library(magrittr); library(lubridate)
library(forecast)
# +++++++++++

ini.dir = getwd()
  dir_outs = "D:/OneDrive/3 - PESQUISA/01 - Producao/Covid/Artigo Brasil/01 Covid_Brasil - Atual/output/output_2021-03-13"
  setwd(dir_outs)
  arqs = list.files(dir_outs, full.names=T)
setwd(ini.dir)

inic = nchar(arqs[1]) - 5
endc = inic + 1
ufs = substr(arqs, inic, endc)
dfs = list()
for(i in 1:(length(arqs))) {
  dfs[[i]] = read_csv(arqs[i], col_names = TRUE)
}


for(i in ufs) {
  df_ = dfs[[which(ufs==i)]]
  plot(df_$total_cases, type = 'l', main = paste(i, ' - DT e Cases'))
  lines(df_$DT, col = 'red')
}

for(i in ufs) {
  df_ = dfs[[which(ufs==i)]]
  plot(df_$total_deaths, type = 'l', main = paste(i, ' - DD e Deaths'))
  lines(df_$DD, col = 'red')
}


