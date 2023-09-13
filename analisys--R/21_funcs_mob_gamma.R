 # mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
# Script:  ----
# For: Functions for script to relate mobulity and contagious rate
# By: Ricardo Aguirre Leal
# Date: May 2021
# mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Util functions  ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function for suppressing output from cat() by Hadley Wickham:
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for get data in files and organize ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Get list of dfs for 'cj_pdreg' (date, gamma_t and mobility_'x'):
getlist_df3s = function(cj) {
  path_cj  = paste0(dir_data, cj, "/")
  files_   = list.files(path_cj)
  files_cj = str_subset(files_, 'df3_') %>% sort()
  ndfs     = str_replace_all(files_cj, 'df3_', '') %>% str_replace_all('.csv', '')
  # +++++++++++
  list_cj = foreach(f=files_cj) %do% {
    df = read_csv(paste0(path_cj, f), col_names = TRUE, col_types = cols()) %>% suppressWarnings()
    df %<>% rename(date = X1)
    df %>% select(date, gamma_t, mobility_L, mobility_W, mobility_F, mobility_P , mobility_T, mobility_R)
  }
  names(list_cj) = ndfs
  # +++++++++++
  cat("List of dataframes for", cj,"created. \n")
  cat("From", length(files_cj), "'df3_' files of", length(files_cj),"total files contained in directory. \n")
  return(list_cj)
}

# Eliminate first and last 'el_x' observations (days) in each df 
# And ensure that observations begin not before 'beg_time' (if beg is NULL, then there is no cutoff)
# Also ensure the same for 'end_time'
elim = function(ldf, el_x) {
  for (i in 1:length(ldf)) {
    ldf[[i]] %<>% slice( (el_x+1):(n()-el_x) )
    if(!is.null(beg_time)) ldf[[i]] %<>% filter(date >= beg_time)
    if(!is.null(end_time)) ldf[[i]] %<>% filter(date <= end_time)
  }
  return(ldf)
}

# Transform variables in percentual, multipling per 100 
in_perc = function(ldf) {
  for (i in 1:length(ldf)) {
    ldf[[i]] %<>% mutate(across(-1, ~ .*100))
  }
  return(ldf)
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for time series preliminaries ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# List whith time series first difference or nivel
makelist_ts = function(ldf, dif) {
  if (dif) {
    for (i in 1:nunits) ldf[[i]] = zoo(ldf[[i]][,-1], ldf[[i]]$date) %>% diff()
  } else {
    for (i in 1:nunits) ldf[[i]] = zoo(ldf[[i]][,-1], ldf[[i]]$date)
  }
  return(ldf)
}

# Unit root test for list of time series
# Performs the augmented Dickey-Fuller unit root test (ADF)
# or Elliott, Rothenberg \& Stock unit root test (ERS)
testlist_ur = function(lts, testur) {
  if(!testur %in% c('adf', 'ers')) stop("Type of test no allowed. Only 'adf' or 'ers'")
  urdf = foreach (iuf = 1:nunits, .combine = rbind) %do% {
    df = lts[[iuf]]
    v_break = rep(NA, nunits)
    urline = foreach (ivaruf = 1:ncol(df), .combine = cbind) %do% {
      if(testur=='adf') {
        urt = ur.df(df[,ivaruf], type = "none", selectlags = "BIC") 
      } else {
        urt = ur.ers(df[,ivaruf], type="DF-GLS", model="const", lag.max=6)
      }
      urt_p  = urt@teststat
      urt_c5 = urt@cval[2] # at 5%
      # H0: unit root. Reject H0 if test-statistic is < critical
      urt = ifelse(urt_p < urt_c5, FALSE, TRUE) # serie is unit root? FALSE or TRUE
    }
  }
  colnames(urdf) = colnames(df)
  urdf = cbind(unit = names(lts), as.data.frame(urdf))
  cat("The quantity of series whith unit root is:", sum(urdf == "TRUE")," \n \n")
  return(urdf)
}

# Unit root test for list of time series
# Zivot-Andrews unit root test, which allows a break at an unknown point in the intercept
testlist_urza = function(lts) {
  cl = detectCores() %>% makePSOCKcluster()
  registerDoParallel(cl)
  nvars = ncol(lts[[1]])
  L = foreach (iuf = 1:nunits) %dopar% {
    require(urca)
    df = lts[[iuf]]
    urline = rep(NA, nvars)
    breakline = rep(as.integer(0), nvars)
    for(j in 1:nvars) {
      urt = ur.za(df[,j], lag=4) 
      urt_p  = urt@teststat
      urt_c5 = urt@cval[2] # at 5%
      breakline[j] = urt@bpoint
      # H0: unit root. Reject H0 if test-statistic is < critical
      urline[j] = ifelse(urt_p < urt_c5, FALSE, TRUE) # serie have unit root? FALSE or TRUE
    }
    list(t(urline), t(breakline))
  }
  urdf     = sapply(L, function(x) x[[1]]) %>% t() %>% as.data.frame()
  df_break = sapply(L, function(x) x[[2]]) %>% t() %>% as.data.frame()
  colnames(urdf) = colnames(df_break) = colnames(lts[[1]])
  urdf     = cbind(unit=names(lts), urdf) 
  df_break = cbind(unit=names(lts), df_break) 
  # +++++
  cat("The quantity of series whith unit root is:", sum(urdf == "TRUE")," \n \n")
  stopCluster(cl)
  return(list(urdf = urdf, df_break = df_break))
}

# makelist_breaks = function(dfur, lurza) {
#   dfurza = lurza[[2]] %>% as.data.frame()
#   vars = colnames(dfur)
#   L = foreach(i = 2:length(vars)) %do% {
#     vectrue = dfurza[dfur[, i], i]
#   }
#   names(L) = vars[-1]
#   indexes = which(!as.matrix(lur[,-1]))
# }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for VAR models ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

makelist_var = function(lts, vars_mobil) {
  vars_used = lts[[1]] %>% colnames() %in% c('gamma_t', vars_mobil)
  L = foreach (iuf = 1:nunits) %do% {
    dfts = lts[[iuf]][, vars_used]
    VAR(dfts, type="const", lag.max=6, ic=criterio_var)
  }
  names(L) = names(lts)
  return(L)
}

makelist_sc = function(lmodel, lrVEC=NULL, test = "efp", print=T) {
  if(!test %in% c('F', 'efp')) stop("Type of test no allowed. Only 'F' or 'efp'")
  L = foreach (iuf = 1:nunits) %do% {
    if(class(lmodel[[1]])=="varest") {
      xx = lmodel[[iuf]] %>% residuals()
    } else {
      xx = cajorls(lmodel[[iuf]], lrVEC[[iuf]])$rlm$residuals
    }
    rowvars = foreach (j = 1:ncol(xx), .combine = cbind) %do% {
      if(test == 'F') {
        res = Fstats(xx[, j] ~ 1) %>% sctest()
      } else {
        res =    efp(xx[, j] ~ 1, type = 'OLS-CUSUM') %>% sctest()
      }
      round(res$p.value, 3)    # H0: no structural change
    }
    rowvars = cbind(Reject = any(rowvars<0.05), rowvars)
    colnames(rowvars) = c("Reject_?", colnames(xx)); rownames(rowvars) = "p-value" 
    rowvars
  }
  names(L) = names(lmodel)
  if(print) {
    df_print = L %>% unlist() %>% matrix(ncol = ncol(rowvars), byrow = T) 
    colnames(df_print) = colnames(rowvars)
    rownames(df_print) = names(lmodel)
    cat("\n There is", sum(df_print<0.05), "p-values below 0.05 \n \n")
    print(df_print)
  }
  lret = list( input = list(test=test, lmodel=lmodel, lrVEC=lrVEC), result = L )
  return(lret)
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for VEC models ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

get_nlags = function(lts, vars_mobil){
  if(!exists('criterio_var')) criterio_var = "SC"
  L = makelist_var(lts, vars_mobil)
  for(i in 1:length(L)) L[[i]] = L[[i]]$p
  return(L)
}

makelist_cajo = function(lts, vars_mobil, llags, print) {
  vars_used = lts[[1]] %>% colnames() %in% c('gamma_t', vars_mobil)
  for(i in 1:nunits) {
    lts[[i]] = ca.jo(lts[[i]][, vars_used], type = "eigen", K = llags[[i]], ecdet = "const")
    if(print) {
      print(names(lts)[i])
      print(cbind(teststat=lts[[i]]@teststat, lts[[i]]@cval))
    }
  }
  return(lts)
}

makelist_rcajo = function(lcajo, pct, print) {
  if(pct==10) a=1 else if(pct==5) a=2 else if(pct==1) a=3 else stop("Chose a valid 'pct'")
  k = length(lcajo[[1]]@teststat)
  for(i in 1:nunits) {
    r = k - match(FALSE, round(lcajo[[i]]@teststat,2) <= round(lcajo[[i]]@cval[,a],2)) + 1 
    lcajo[[i]] = r
    if(is.na(r)) lcajo[[i]] = 'zero_(VAR_diff)' else if(r==k) lcajo[[i]] = 'full_(VAR_nivel)'
  }
  if(print) {
    cat("Quantity of vectors of cointegration: \n")
    lcajo %>% as.matrix() %>% print()
  }
  return(lcajo)
}

makelist_lttest = function(lcajo, lrcajo, print) {
  for(i in 1:nunits) {
    if(is.numeric(lrcajo[[i]])) {
      tt = lttest(lcajo[[i]], lrcajo[[i]]) %>% quiet() 
      lcajo[[i]] = tt[2]
    } else {
      lcajo[[i]] = ifelse(lrcajo[[i]]=='zero_(VAR_diff)', yes = 'VAR_diff', no = 'VAR_nivel')
    }  
  }
  if(print) {
    cat("p-value of H0 for not including a linear trend: \n")
    lcajo %>% as.matrix() %>% print()
  }
  return(lcajo)
}

plot_ltscVEC = function(ltscVEC) {
  test   = ltscVEC$input$test
  lmodel = ltscVEC$input$lmodel 
  lrVEC  = ltscVEC$input$lrVEC
  lrestest = ltscVEC$result
  # +++
  for(iuf in 1:nunits) {
    if(lrestest[[iuf]][1]==0) next
    xx = cajorls(lmodel[[iuf]], lrVEC[[iuf]])$rlm$residuals
    jlog = which(lrestest[[iuf]][-1] < 0.05)
    for(j in jlog) {
      if(test == 'F') {
        res = Fstats(xx[, j] ~ 1)
        texttest = "Fstats --"
      } else {
        res = efp(xx[, j] ~ 1, type = 'OLS-CUSUM')
        texttest = "OLS-CUSUM --"
      }
      plot(res, main = paste(texttest, names(lrestest)[iuf], "--" ,colnames(xx)[j]))
      plot(xx[, j], type = "l", ylab = "Residual",
           main = paste("Residual --", names(lrestest)[iuf], "--" ,colnames(xx)[j])); abline(h=0, col=2)
    }
  }
}

makelist_scVEC = function(ltscVEC) {
  test     = ltscVEC$input$test
  lmodel   = ltscVEC$input$lmodel 
  lrVEC    = ltscVEC$input$lrVEC
  lrestest = ltscVEC$result
  nseries  = ncol(lrestest[[1]]) - 1
  # ++++++
  breaks = foreach(iuf = 1:nunits) %do% {
    if(lrestest[[iuf]][1]==1) {
      xx   = cajorls(lmodel[[iuf]], lrVEC[[iuf]])$rlm$residuals
      jlog = which(lrestest[[iuf]][-1] < 0.05)
      # ++++++
      min_grid = Inf; max_grid = 0
      for(j in jlog) {
        if(test == 'F') {
          sctest = Fstats(xx[, j] ~ 1)
          bs    = strucchange::boundary(sctest)
          grid_obs = which(sctest$Fstats  > bs | sctest$Fstats  < (-bs))
        } else {
          sctest = efp(xx[, j] ~ 1, type='OLS-CUSUM')
          bs    = strucchange::boundary(sctest)
          grid_obs = which(sctest$process > bs | sctest$process < (-bs))
        }
        min_grid = min(min_grid, grid_obs)
        max_grid = max(max_grid, grid_obs)
      }
      grid_obs = seq(round(min_grid*0.95, 0), round(max_grid*0.95, 0))
      # ++++++
      nobs = lmodel[[iuf]]@x %>% nrow()
      Mres = foreach(j = grid_obs, .combine = rbind) %do% {
        dummie = rep(0, nobs)
        dummie[j:nobs] = 1
        xcajo = ca.jo(lmodel[[iuf]]@x, type="eigen", K=lmodel[[iuf]]@lag, ecdet="const", 
                      dumvar=matrix(dummie, dimnames = list(NULL,"dummie_break")))
        bic = vec2var(xcajo, lrVEC[[iuf]]) %>% logLik() %>% BIC()
        cbind(breakpoint = j, bic = bic)
      }
      L = list(breakpoint = Mres[which.min(Mres[,2]), 1], bic_grid = Mres)
    } else L = NULL
    L
  }
  names(breaks) = names(lmodel)
  # ++++++
  lVEC = foreach(iuf = 1:nunits) %do% {
    nobs = lmodel[[iuf]]@x %>% nrow()
    if(!is.null(breaks[[iuf]])) {
      dummie = rep(0, nobs); dummie[breaks[[iuf]]$breakpoint:nobs] = 1
      xcajo = ca.jo(lmodel[[iuf]]@x, type="eigen", K=lmodel[[iuf]]@lag, ecdet="const", 
                    dumvar=matrix(dummie, dimnames = list(NULL,"dummie_break")))

    } else xcajo = ca.jo(lmodel[[iuf]]@x, type="eigen", K=lmodel[[iuf]]@lag, ecdet="const")
    xcajo
  }
  lVEC2VAR = lapply(lVEC, vec2var, lrVEC[[iuf]])
  lVEC     = lapply(lVEC, cajorls, lrVEC[[iuf]])
  names(lVEC) = names(lVEC2VAR) = names(lmodel)
  return(list(lVEC=lVEC, lVEC2VAR=lVEC2VAR, breaks=breaks))
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Functions for TVEC models ----
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++









