library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(here)
library(EBASE)
library(doParallel)
library(ggplot2)
library(patchwork)
library(oce)
library(SWMPr)
library(WtRegDO)

source(file = here('R/funcs.R'))

fwdat <- read_csv(here("data/apafwoxy2.csv")) 

# fwoxy for comparison
# convert areal to volumetric, extract b at right time step
fwdatcmp <- fwdat %>% 
  mutate(
    DateTimeStamp = dmy_hms(datet, tz = 'America/Jamaica'),
    Date = as.Date(DateTimeStamp, tz = 'America/Jamaica'),
    DO_obs = `oxy,mmol/m3`, 
    a = `aparam,(mmolO2/m2/d)/(W/m2)`,
    R = `er,mmol/m2/d`,
    P = `gpp,mmol/m2/d`,
    D = -1 * `gasex,mmol/m2/d`
  ) %>% 
  select(Date, DateTimeStamp, DO_obs, a, b = `bparam, (cm/hr)/(m2/s2)`, P, R, D)

# fwoxy for input to ebase
fwdatinp <- fwdat %>% 
  mutate(
    datet = dmy_hms(datet, tz = 'America/Jamaica')
  ) %>% 
  select(
    DateTimeStamp = datet,
    DO_obs = `oxy,mmol/m3`, 
    Temp = `temp,degC`, 
    Sal = `salt,ppt`, 
    PAR = `par,W/m2`, 
    WSpd = `wspd2,m2/s2`, 
    H = `ht,m`
  ) %>% 
  mutate(
    DO_obs = DO_obs / 1000 * 32, # to mg/L
    WSpd = sqrt(WSpd)
  )

# recovery default priors 7, 30 day opt -------------------------------------------------------

# 7 days opt
cl <- makeCluster(6)
registerDoParallel(cl)
  
res7 <- ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
             ndays = 7)
  
stopCluster(cl)

# 30 days opt
cl <- makeCluster(6)
registerDoParallel(cl)

res30 <- ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
              ndays = 30)

stopCluster(cl)

apadef <- list(
  opt7 = res7, 
  opt30 = res30
)

save(apadef, file = 'data/apadef.RData', compress = 'xz')

# recovery default priors 7, 30 day opt, fixed b ----------------------------------------------

# 7 days opt
cl <- makeCluster(6)
registerDoParallel(cl)

res7 <- ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
              ndays = 7, bprior = c(0.251, 1e-6))

stopCluster(cl)

# 30 days opt
cl <- makeCluster(6)
registerDoParallel(cl)

res30 <- ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
               ndays = 30, bprior = c(0.251, 1e-6))

stopCluster(cl)

apadeffix <- list(
  opt7 = res7, 
  opt30 = res30
)

save(apadeffix, file = 'data/apadeffix.RData', compress = 'xz')

# gridded comparisons, mean and sd, 1 day -----------------------------------------------------

# this takes about 24 hours to run
grd <- crossing(
  amean = c(2, 12),
  asd = c(0.4, 40),
  rmean = c(150, 900), 
  rsd = c(30, 3000),
  ndays = c(1),
  out = NA
)

str <- Sys.time()

# takes about 24 hours
for(i in 1:nrow(grd)){
  
  # counter
  cat(i, 'of', nrow(grd), '\n')
  print(Sys.time() - str)
  
  # get inputs
  selrow <- grd[i, ]
  aprior <- c(selrow$amean, selrow$asd)
  rprior <- c(selrow$rmean, selrow$rsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  ncores <- detectCores() - 2
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- try(ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
               aprior = aprior, rprior = rprior, bprior = c(0.251, 1e-6), ndays = ndays), 
             silent = T)

  stopCluster(cl)

  ind <- 1
  while(inherits(res, 'try-error') & ind <= 5){
    
    cat('retrying...\t')
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    res <- try(ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
                aprior = aprior, rprior = rprior, bprior = c(0.251, 1e-6), ndays = ndays), 
               silent = T)
    
    stopCluster(cl)
    
    ind <- ind + 1
    
  }
  if(ind > 5){
    cat('failed...\t')
    next()
  }
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

apagrd1 <- grd
save(apagrd1, file = 'data/apagrd1.RData', compress = 'xz')

# gridded comparisons, mean and sd, 7 day -----------------------------------------------------

# this takes about 24 hours to run
grd <- crossing(
  amean = c(2, 12),
  asd = c(0.4, 40),
  rmean = c(150, 900), 
  rsd = c(30, 3000),
  ndays = c(7),
  out = NA
)

str <- Sys.time()

# takes about 20 hours
for(i in 1:nrow(grd)){
  
  # counter
  cat(i, 'of', nrow(grd), '\n')
  print(Sys.time() - str)
  
  # get inputs
  selrow <- grd[i, ]
  aprior <- c(selrow$amean, selrow$asd)
  rprior <- c(selrow$rmean, selrow$rsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  ncores <- detectCores() - 2
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- try(ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
                   aprior = aprior, rprior = rprior, bprior = c(0.251, 1e-6), ndays = ndays), 
             silent = T)
  
  stopCluster(cl)
  
  ind <- 1
  while(inherits(res, 'try-error') & ind <= 5){
    
    cat('retrying...\t')
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    res <- try(ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
                     aprior = aprior, rprior = rprior, bprior = c(0.251, 1e-6), ndays = ndays), 
               silent = T)
    
    stopCluster(cl)
    
    ind <- ind + 1
    
  }
  if(ind > 5){
    cat('failed...\t')
    next()
  }
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

# load(file = 'data/agrd7to12.RData')
# apagrd7to12$out[13:16] <- grd$out[13:16]
# apagrd7 <- apagrd7to12
apagrd7 <- grd
save(apagrd7, file = 'data/apagrd7.RData', compress = 'xz')

# gridded comparisons, mean and sd, 30 day ----------------------------------------------------

# this takes about 48 hours to run
grd <- crossing(
  amean = c(2, 12),
  asd = c(0.4, 40),
  rmean = c(150, 900), 
  rsd = c(30, 3000),
  ndays = c(30),
  out = NA
)

str <- Sys.time()

# takes about 48 hours
for(i in 1:nrow(grd)){
  
  # counter
  cat(i, 'of', nrow(grd), '\n')
  print(Sys.time() - str)
  
  # get inputs
  selrow <- grd[i, ]
  aprior <- c(selrow$amean, selrow$asd)
  rprior <- c(selrow$rmean, selrow$rsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  ncores <- detectCores() - 2
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- try(ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
                   aprior = aprior, rprior = rprior, bprior = c(0.251, 1e-6), ndays = ndays), 
             silent = T)
  
  stopCluster(cl)
  
  ind <- 1
  while(inherits(res, 'try-error') & ind <= 5){
    
    cat('retrying...\t')
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    res <- try(ebase(fwdatinp, interval = 900, Z = fwdatinp$H, n.chains = 4, 
                     aprior = aprior, rprior = rprior, bprior = c(0.251, 1e-6), ndays = ndays), 
               silent = T)
    
    stopCluster(cl)
    
    ind <- ind + 1
    
  }
  if(ind > 5){
    cat('failed...\t')
    next()
  }
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

# load(file = 'data/apagrd30to8.RData')
# apagrd30to8$out[9:16] <- grd$out[9:16]
# apagrd30 <- apagrd30to8
apagrd30 <- grd
save(apagrd30, file = 'data/apagrd30.RData', compress = 'xz')

# evaluate fit all priors -------------------------------------------------

load(file = here('data/apagrd1.RData'))
load(file = here('data/apagrd7.RData'))
load(file = here('data/apagrd30.RData'))

apagrd <- bind_rows(apagrd1, apagrd7, apagrd30)

# # view incomplete groups, the last one
# lapply(apagrd$out, function(x) table(x[[1]]$grp))
apasumdat <- apagrd %>% 
  mutate(
    ind = 1:nrow(.),
    ests = purrr::pmap(list(ind, out), function(ind, out){
      
      cat(ind, '\t')

      # filter incomplete groups for 7, 30 day opt, the last one
      if(!ind %in% 1:16)
        out <- out[[1]] %>% 
          filter(grp != max(grp))
      else 
        out <- out[[1]]
      
      cmp <- inner_join(fwdatcmp, out, by = c('Date', 'DateTimeStamp')) %>%
        select(-Z, -converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>%
        rename(
          DO_mod.x = DO_obs.x,
          DO_mod.y = DO_mod
        ) %>%
        pivot_longer(!all_of(c('DateTimeStamp', 'Date', 'grp')), names_to = 'var', values_to = 'val') %>%
        separate(var, c('var', 'mod'), sep = '\\.') %>%
        mutate(
          mod = case_when(
            mod == 'x' ~ 'Fwoxy',
            mod == 'y' ~ 'EBASE'
          )
        ) %>%
        pivot_wider(names_from = 'mod', values_from = 'val')
      
      # average output by optimization period
      sumgrp <- cmp %>% 
        group_by(Date, var) %>% 
        summarise(
          Fwoxy = mean(Fwoxy, na.rm = T), 
          EBASE = mean(EBASE, na.rm = T), 
          .groups = 'drop'
        ) 
      
      sumcmp <- sumgrp %>% 
        group_by(var) %>% 
        nest() %>% 
        summarise(
          est = purrr::map(data, sumfun)
        ) %>% 
        unnest('est')
      
      return(sumcmp)
      
    })
  ) %>% 
  select(-out)

save(apasumdat, file = here('data/apasumdat.RData'), compress = 'xz')

# adding noise to fwoxy time series -----------------------------------------------------------

# for PAR conversion to W/m2
Jpmolph <-  0.2175e6 # 1 mol-photons = 0.2175e6 J for ave PAR wavelength of 550nm

# used weighted regression on 2021 apacp
apacpwq <- import_local('data/apa2021.zip', station_code = 'apacpwq') %>%
  qaqc(qaqc_keep = as.character(seq(-5, 5)))
apaebmet <- import_local('data/apa2021.zip', station_code = 'apaebmet') %>%
  qaqc(qaqc_keep = as.character(seq(-5, 5)))
apacp <- comb(apacpwq, apaebmet) %>%
  mutate(
    PAR = totpar * Jpmolph * 1e-3 / 15/ 60, # mmol/m2/15min to W/m2
  ) %>% 
  select(
    DateTimeStamp = datetimestamp,
    Temp = temp,
    Sal = sal,
    DO_obs = do_mgl,
    ATemp = atemp,
    BP = bp,
    WSpd = wspd,
    PAR,
    Tide = depth
  ) %>%
  filter(!is.na(Tide) | !is.na(DO_obs))

apacpdtd <- wtreg(apacp, wins = list(12, 12, 0.4), tz = 'America/Jamaica', lat = 29.7021, long = -84.8802)

save(apacpdtd, file = here('data/apacpdtd.RData'))

# get height data - not in exdat
load(file = here('data/apacpdtd.RData'))

nosdat <- apacpdtd %>% 
  mutate(
    tidnoise = DO_prd - DO_nrm, 
    obsnoise = DO_obs - DO_prd
  ) %>% 
  select(DateTimeStamp, tidnoise, obsnoise)

# add tidal noise to fwdatinp - need to figure out values less than zero (<1% of obs)
fwdatinpnos <- fwdatinp %>% 
  left_join(nosdat, by = 'DateTimeStamp') %>% 
  mutate(
    DO_nos = pmax(0, DO_obs + tidnoise + obsnoise)
  ) %>% 
  na.omit()

tomod <- fwdatinpnos %>% 
  select(-tidnoise, -obsnoise, -DO_obs) %>% 
  rename(DO_obs = DO_nos)

# run model for inputs
cl <- makeCluster(10)
registerDoParallel(cl)

# use interp for missing values
resnos <- ebase(tomod, interval = 900, Z = tomod$H, bprior = c(0.251, 1e-6), n.chains = 4, ndays = 7)

# save ebase results for noisy ts
save(resnos, file = here('data/resnos.RData'))

tomod <- fwdatinpnos %>% 
  select(-obsnoise, -tidnoise, -DO_nos) 

# run model for inputs
cl <- makeCluster(10)
registerDoParallel(cl)

# use interp for missing values
resobs <- ebase(tomod, interval = 900, Z = tomod$H, bprior = c(0.251, 1e-6), n.chains = 4, ndays = 7)

# save ebase results for observed (synthetic) ts
save(resobs, file = here('data/resobs.RData'))

# missing data summary for apacp --------------------------------------------------------------

# for PAR conversion to W/m2
Jpmolph <-  0.2175e6 # 1 mol-photons = 0.2175e6 J for ave PAR wavelength of 550nm

# used weighted regression on 2021 apacp
apacpwq <- import_local('data/apa2021.zip', station_code = 'apacpwq') %>%
  qaqc(qaqc_keep = as.character(seq(-5, 5)))
apaebmet <- import_local('data/apa2021.zip', station_code = 'apaebmet') %>%
  qaqc(qaqc_keep = as.character(seq(-5, 5)))
cmbdat <- comb(apacpwq, apaebmet) %>%
  mutate(
    PAR = totpar * Jpmolph * 1e-3 / 15/ 60, # mmol/m2/15min to W/m2
  ) %>% 
  select(
    DateTimeStamp = datetimestamp,
    Temp = temp,
    Sal = sal,
    DO_obs = do_mgl,
    WSpd = wspd,
    PAR
  )
misdat <- cmbdat %>% 
  reframe(across(-DateTimeStamp, function(x) sum(is.na(x)) / length(x))) %>% 
  mutate(across(everything(), function(x) round(100 * x, 1)))

save(misdat, file = here('data/misdat.RData'))          