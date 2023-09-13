# mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
# Script:  ----
# For: Functions for script to create data files for use in SEIR model
# By: Ricardo Aguirre Leal
# Date: May 2021
# mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for get data in files ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# For covid-19 cases and deaths 
getdf_covid = function() {
  df_covid = read_csv2(paste0(dir_data, "HIST_PAINEL_COVIDBR_", cj_data, '.csv'), 
                       col_names = TRUE, col_types = "cc--i--Diiiiii---")
  df_covid %<>% filter(regiao!="Brasil", is.na(codmun), data>=data_ini) %>% select(-codmun)
  names(df_covid)[2] = "uf"
  names(df_covid)[4:9] = c("semana", "population", "total_cases", "new_cases", "total_deaths", "new_deaths")
  # +++++++++++
  cat("\n Dataframe for covid-19 cases and deaths created \n")
  return(df_covid)
}

# For covid-19 vaccines
getdf_vac = function() {
  # escolher "estabelecimento_uf" porque "paciente_endereco_uf" apresenta "XX" e NA
  locvars_vac = c(20, 28, 29)   # "estabelecimento_uf", "vacina_dataAplicacao", "vacina_descricao_dose"
  cty = rep("-", 34)
  cty[locvars_vac] = rep("c", 3)
  # +++++++++++
  df_vac = read_csv2(paste0(dir_vac, "completo_", cj_data,'.csv'), 
                     col_types = paste0(cty, collapse=""), col_names = TRUE)
  names(df_vac) = c("uf", "data", "dose")
  df_vac %<>% mutate(dose = str_trim(dose))
  xDs = unique(df_vac$dose) %>% sort()  # "1ª Dose" "2ª Dose" "Dose"    "Única"
  df_vac %<>% mutate(data = as_date(str_sub(data, 1, 10)))
  # +++++++++++
  df_vac_one = df_vac[df_vac$dose %in% xDs[1], ] %>% group_by(uf, data) %>% 
    summarise(people_vaccinated = n())  %>% arrange(uf, data) %>% ungroup()
  df_vac_full = df_vac[df_vac$dose %in% xDs[c(2, 4)], ] %>% group_by(uf, data) %>% 
    summarise(people_fully_vaccinated = n()) %>% arrange(uf, data) %>% ungroup()  
  df_vac2 = left_join(df_vac_one, df_vac_full) %>% arrange(uf, data)
  df_vac2$people_fully_vaccinated %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  df_vac2$people_vaccinated %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  df_vac2 %<>% group_by(uf) %>% mutate(people_vaccinated = cumsum(people_vaccinated)) %>%
    arrange(uf, data) %>% ungroup()
  df_vac2 %<>% group_by(uf) %>% mutate(people_fully_vaccinated = cumsum(people_fully_vaccinated)) %>%
    arrange(uf, data) %>% ungroup()
  df_vac2$total_vaccinations = df_vac2$people_vaccinated + df_vac2$people_fully_vaccinated
  # +++++++++++
  cat("\n Dataframe for covid-19 vaccines created \n")
  beep()
  return(df_vac2)
}

getdf_mobil = function() {
  cty = c(rep("-", 5), "c", "-", "-", "D", rep("n", 6)) 
  df_mobil20 = read_csv(paste0(dir_data, "2020_BR_Region_Mobility_Report_", cj_data, '.csv'), col_names = TRUE, 
                        col_types = paste0(cty, collapse = ""))
  df_mobil21 = read_csv(paste0(dir_data, "2021_BR_Region_Mobility_Report_", cj_data, '.csv'), col_names = TRUE, 
                        col_types = paste0(cty, collapse = ""))
  df_mobil = bind_rows(df_mobil20, df_mobil21)
  names(df_mobil) = c("uf", "data", "google_retail_and_recreation", "google_grocery_and_pharmacy", 
                      "google_parks", "google_transit", "google_workplaces", "google_residential")
  df_mobil %<>% filter(!is.na(uf)) %>% arrange(uf, data) 
  df_mobil %<>% filter(data >= data_ini)
  df_mobil %<>% mutate(uf = str_sub(uf, 4))
  # +++++++++++
  df_mobil$google_retail_and_recreation %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  df_mobil$google_grocery_and_pharmacy %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  df_mobil$google_parks %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  df_mobil$google_transit %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  df_mobil$google_workplaces %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  df_mobil$google_residential %<>% str_replace_na() %>% str_replace('NA', "0") %>% as.integer()
  # +++++++++++
  cat("\n Dataframe for Google mobility created \n")
  return(df_mobil)
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions check and adjust datas ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# For covid-19 cases and deaths 
check_df_covid = function(view_neg) {
  xx = df_covid %>% mutate( ldata=lag(.$data)) %>% group_by(uf) %>% mutate(ddate = data - ldata)
  res = (xx$data == xx$ldata) %>% sum(na.rm = T)
  if(res!=0) { cat("\n \n ERROR: The states have not the same dates \n")
  } else {
    cat("\n OK: All the states have the same dates \n")  }
  cat("\n Number of NA in data frame: ", sum(is.na(df_covid)), " \n")
  cat("\n Number of negative values in data frame: ", sum(df_covid<0, na.rm = TRUE), " \n \n")
  # +++++++++++
  if(view_neg) {
    df_covid[(which(df_covid<0, TRUE)[,1]), ]
    df_covid[(which(df_covid<0, TRUE)[,1]), ] %>% View()
  }
}

# For covid-19 vaccines 
check_df_vac = function(view_neg) {
  cat("\n \n Number of vaccinated before 2021: ", sum(year(df_vac$data) < 2021), "\n")
  cat("\n Number of vaccinated before 2021-01-18 (official date start vaccines): ", 
      sum(df_vac$data < "2021-01-18"), "\n")
  cat("\n Number of NA in data frame: ", sum(is.na(df_vac)), " \n")
  cat("\n Number of negative values in data frame: ", sum(df_vac<0, na.rm = TRUE), " \n \n")
  # +++++++++++
  if(view_neg) {
    df_vac[(which(df_vac<0, TRUE)[,1]), ]
    df_vac[(which(df_vac<0, TRUE)[,1]), ] %>% View()
  }
}

# For Google mobility 
check_df_mobil = function(view_neg) {
  xx = df_mobil %>% mutate(ldata=lag(.$data)) %>% group_by(uf) %>% mutate(ddate = data - ldata)
  res = (xx$data == xx$ldata) %>% sum(na.rm = T)
  if(res!=0) { cat("\n \n ERROR: The states have not the same dates \n")
  } else {
    cat("\n OK: All the states have the same dates \n")  }
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function for seasonally adjust data ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

adj_df_seas = function() {
  sfobs = epiweek(max(dados$date)) - 1
  seman_fim = epiweek('2020-12-31') + sfobs  
  semanas = rep((10:seman_fim), each=7, times=27)
  selfilter = which( 
    (year(dados$date)==2020) |                             
      (year(dados$date)==2021 & dados$semana==53) |        
      (year(dados$date)==2021 & dados$semana<=sfobs)
  ) 
  dados2 = dados[selfilter,] # dates are shortened for only entire weeks
  dados2 %<>% mutate(semana = semanas)
  semanas = unique(semanas) %>% sort()
  # +++++++++++
  tscov = function(sts) {
    sts = as.numeric(unlist(sts))
    msts = msts(sts, seasonal.periods = c(7),start = decimal_date(as.Date(data_ini)))
    msts = decompose(msts)$trend
    tira = c(1:7, sort(1:(length(sts)), T)[1:7] )
    msts = msts[-tira] 
    msts
  }
  # +++++++++++
  dados_seas = dados2 %>% select(-total_vaccinations_per_hundred, -people_vaccinated_per_hundred,
                                 -people_fully_vaccinated_per_hundred)
  dados_seas %<>% filter(semana != min(semanas) & semana != max(semanas))
  colfor = dados_seas %>% select(-region, -iso2, -date, -semana, -population) %>% colnames()
  dados_seas %<>% mutate_at(all_of(colfor), as.double)
  colfor = dados_seas %>% select(-region, -iso2, -date, -semana, -population, 
                                 -total_vaccinations, -total_cases, -total_deaths) %>% colnames()
  for(c in colfor) {
    for(uf in ufs) {
      linhas = dados2$iso2 == uf
      newdata = tscov(dados2[linhas,c])
      linhas = dados_seas$iso2 == uf
      dados_seas[linhas,c] = newdata
    }
  }
  for(uf in ufs) {
    linhas = dados_seas$iso2 == uf
    newdata_cases  = dados_seas$new_cases[linhas] %>% cumsum()
    newdata_deaths = dados_seas$new_deaths[linhas] %>% cumsum()
    newdata_vacc   = dados_seas$people_vaccinated[linhas] + dados_seas$people_fully_vaccinated[linhas]
    dados_seas[linhas,'total_cases'] = newdata_cases
    dados_seas[linhas,'total_deaths'] = newdata_deaths
    dados_seas[linhas,'total_vaccinations'] = newdata_vacc
  }
  dados_seas %<>% mutate(total_vaccinations_per_hundred = total_vaccinations / (population/100) )
  dados_seas %<>% mutate(people_vaccinated_per_hundred  = people_vaccinated / (population/100) )
  dados_seas %<>% mutate(people_fully_vaccinated_per_hundred = people_fully_vaccinated / (population/100) )
  return(dados_seas)
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function for smooth data whith splines ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

adj_df_splines = function(s_par) {
  dados_spline = dados %>% select(-total_vaccinations_per_hundred, -people_vaccinated_per_hundred,
                                  -people_fully_vaccinated_per_hundred)
  colfor = dados_spline %>% select(-region, -iso2, -date, -semana, -population) %>% colnames()
  dados_spline %<>% mutate_at(all_of(colfor), as.double)
  colfor = dados_spline %>% select(-region, -iso2, -date, -semana, -population, 
                                   -total_vaccinations, -total_cases, -total_deaths) %>% colnames()
  for(colvar in colfor) {
    for(uf in ufs) {
      linhas = dados$iso2 == uf
      newdata = smooth.spline(dados[linhas, colvar], spar = s_par)$y %>% round(5)
      if(!str_detect(colvar, 'google')) newdata %<>% abs()
      linhas = dados_spline$iso2 == uf
      dados_spline[linhas,colvar] = newdata
    }  
  }
  for(uf in ufs) {
    linhas = dados_spline$iso2 == uf
    newdata_cases  = dados_spline$new_cases[linhas] %>% cumsum()
    newdata_deaths = dados_spline$new_deaths[linhas] %>% cumsum()
    newdata_vacc = dados_spline$people_vaccinated[linhas] + dados_spline$people_fully_vaccinated[linhas]
    dados_spline[linhas,'total_cases'] = newdata_cases
    dados_spline[linhas,'total_deaths'] = newdata_deaths
    dados_spline[linhas,'total_vaccinations'] = newdata_vacc
  }
  dados_spline %<>% mutate(total_vaccinations_per_hundred = total_vaccinations / (population/100) )
  dados_spline %<>% mutate(people_vaccinated_per_hundred  = people_vaccinated / (population/100) )
  dados_spline %<>% mutate(people_fully_vaccinated_per_hundred = people_fully_vaccinated / (population/100) )
  return(dados_spline)
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for regional datas ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

group_reg = function(df, days) {
  # Function to calculate regional lines 
  linhaf = function(dadosd) {
    linha = c(reg, 
              d, 
              dadosd$semana[1],
              sum(dadosd$population), 
              sum(dadosd$total_cases), 
              sum(dadosd$new_cases), 
              sum(dadosd$total_deaths), 
              sum(dadosd$new_deaths), 
              sum(dadosd$people_vaccinated), 
              sum(dadosd$people_fully_vaccinated),
              sum(dadosd$total_vaccinations), 
              mean(dadosd$google_retail_and_recreation), 
              mean(dadosd$google_grocery_and_pharmacy), 
              mean(dadosd$google_parks), 
              mean(dadosd$google_transit), 
              mean(dadosd$google_workplaces), 
              mean(dadosd$google_residential) )
    return(linha)
  }
  # +++++++++++
  dados_grouped = foreach(reg=regioes, .combine = rbind) %do% {
    # exclude federal district (DF), already include in Goias state:
    dadosR = df %>% filter(iso2 != 'DF') %>% filter(region==reg) 
    dadosR_d = foreach(d=days, .combine = rbind) %do% {
      dadosd = dadosR %>% filter(date==d) %>% as.data.frame()
      linhaf(dadosd)
    }
  }
  colnames(dados_grouped) = df %>% select(-iso2, -total_vaccinations_per_hundred, 
                                          -people_vaccinated_per_hundred, 
                                          -people_fully_vaccinated_per_hundred) %>% colnames()
  dados_grouped = as_tibble(dados_grouped)
  dados_grouped %<>% mutate(date = as_date(as.numeric(date)))
  dados_grouped %<>% mutate(across(3:ncol(dados_grouped), as.double))
  # +++++++++++
  reg = 'Brasil'
  dadosBR_d = foreach(d=days, .combine = rbind) %do% {
    dadosd = dados_grouped %>% filter(date==d) %>% as.data.frame()
    linhaf(dadosd)
  }
  colnames(dadosBR_d) = colnames(dados_grouped)
  dadosBR_d = as_tibble(dadosBR_d)
  dadosBR_d %<>% mutate(date = as_date(as.numeric(date)))
  dadosBR_d %<>% mutate(across(3:ncol(dadosBR_d), as.double))
  # +++++++++++
  dados_grouped = bind_rows(dados_grouped, dadosBR_d)
  dados_grouped %<>% rename(iso2 = region)
  # +++++++++++
  dados_grouped %<>% mutate(total_vaccinations_per_hundred = total_vaccinations / (population/100) )
  dados_grouped %<>% mutate(people_vaccinated_per_hundred  = people_vaccinated / (population/100) )
  dados_grouped %<>% mutate(people_fully_vaccinated_per_hundred = people_fully_vaccinated / (population/100) )  # +++++++++++
  cat("\n Dataframe for regional data created \n")
  return(dados_grouped)
}


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for proportion of notified cases ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

make_df_epi = function(df) {
  # start and end dates in round of EPICOVID survey
  ini = as.Date(c("2020/05/14", "2020/06/04", "2020/06/21")) 
  end = as.Date(c("2020/05/21", "2020/06/07", "2020/06/24"))
  # +++++++++++
  df_epi = df %>% filter(date >= ini[1] & date <= end[1])
  df_epi %<>% group_by(iso2) %>% summarise(max = max(total_cases), population = max(population))
  df_epi %<>% mutate('p_casos_R1' = 100*max/population) %>% select(-2, -3)
  # +++++++++++
  df_epi2 = df %>% filter(date >= ini[2] & date <= end[2])
  df_epi2 %<>% group_by(iso2) %>% summarise(max = max(total_cases), population = max(population))
  df_epi2 %<>% transmute('p_casos_R2' = 100*max/population)
  # +++++++++++
  df_epi3 = df %>% filter(date >= ini[3] & date <= end[3])
  df_epi3 %<>% group_by(iso2) %>% summarise(max = max(total_cases), population = max(population))
  df_epi3 %<>% transmute('p_casos_R3' = 100*max/population)
  # +++++++++++
  df_epi %<>% bind_cols(df_epi2, df_epi3)
  return(df_epi)
}

include_reg_df_exp = function()  {
  uf_regiao = dados %>% select(iso2, region) %>% group_by(iso2) %>% 
    summarise(region = first(region)) %>% arrange(iso2)
  df_exp_reg = bind_cols( iso2 = uf_regiao$region, df_exp[-1]) 
  df_exp_reg %<>% filter(iso2 != 'DF') %>% group_by(iso2) %>% summarise_all(mean)
  df_exp_BR = df_exp %>% filter(iso2 != 'DF') %>% summarise_at(2:(ncol(df_exp)), mean) %>% 
                         bind_cols(iso2 = 'Brasil', .)
  df_exp_reg %<>% bind_rows(df_exp_BR)
  df_exp_reg = bind_cols(region = df_exp_reg$iso2, df_exp_reg)
  # include region mean cols for states:
  df_exp = bind_cols(region = uf_regiao$region, df_exp)
  df_exp %<>% bind_rows(df_exp_reg)
  df_exp$p_notif_region_R23 = NA
  df_exp$p_notif_region_R3  = NA
  regions_plus = c(regioes, 'Brasil')
  for (r in regions_plus) {
    linhas = which(df_exp$region == r)
    valueR23 = df_exp_reg %>% filter(iso2 == r) %>% select(p_notif_R23) %>% as.double()
    valueR3  = df_exp_reg %>% filter(iso2 == r) %>% select(p_notif_R3)  %>% as.double()
    df_exp$p_notif_region_R23[linhas] = valueR23
    df_exp$p_notif_region_R3[linhas]  = valueR3
  }
  df_exp %<>% select(-region)
  return(df_exp)
}






