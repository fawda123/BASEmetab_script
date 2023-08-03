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
    # DateTimeStamp = floor_date(DateTimeStamp, unit = 'hours'),
    DO_obs = DO_obs / 1000 * 32, # to mg/L
    WSpd = sqrt(WSpd)
  ) #%>% 
# summarise(across(everything(), ~mean(.x, na.rm = T)), .by = 'DateTimeStamp')

# # simple comparison ---------------------------------------------------------------------------
#
# # subset four days in June
# dat <- fwdatinp %>%
#   filter(month(fwdatinp$DateTimeStamp) == 6 & day(fwdatinp$DateTimeStamp) %in% 1:4)
# 
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# res <- ebase(dat, interval = 900, H = dat$H, progress = TRUE, n.chains = 4, interp = T)
# 
# stopCluster(cl)
# 
# fwdatinpcmp <- fwdatcmp %>% 
#   filter(Date <= max(res$Date) & Date >= min(res$Date))
# 
# cmp <- inner_join(fwdatcmp, res, by = c('Date', 'DateTimeStamp')) %>% 
#   select(-converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>% 
#   rename(
#     DO_mod.x = DO_obs.x, 
#     DO_mod.y = DO_mod
#   ) %>% 
#   pivot_longer(!all_of(c('DateTimeStamp', 'Date', 'grp')), names_to = 'var', values_to = 'val') %>% 
#   separate(var, c('var', 'mod'), sep = '\\.') %>% 
#   mutate(
#     mod = case_when(
#       mod == 'x' ~ 'Fwoxy', 
#       mod == 'y' ~ 'EBASE'
#     )
#   ) %>% 
#   pivot_wider(names_from = 'mod', values_from = 'val')
# 
# ggplot(cmp, aes(x = Fwoxy, y = EBASE, color = var)) + 
#   geom_point() +
#   facet_wrap(~var, ncol = 2, scales = 'free') +
#   theme_bw() + 
#   geom_abline(intercept = 0, slope = 1) +
#   theme(
#     legend.title = element_blank(),
#     strip.background = element_blank()
#   )

# gridded comparisons, mean and sd, 1 day -----------------------------------------------------

# this takes about 24 hours to run
grd <- crossing(
  amean = c(0, 2),
  asd = c(0.01, 1),
  rmean = c(0, 200), 
  rsd = c(0.5, 50),
  bmean = c(0, 0.502), 
  bsd = c(0.001, 0.1),
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
  bprior <- c(selrow$bmean, selrow$bsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- ebase(fwdatinp, interval = 900, H = fwdatinp$H, progress = TRUE, n.chains = 4, 
               aprior = aprior, rprior = rprior, bprior = bprior, ndays = ndays)
  
  stopCluster(cl)
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

apagrd <- grd
apagrd1a <- apagrd[1:32,]
apagrd1b <- apagrd[33:64,]
save(apagrd1a, file = 'data/apagrd1a.RData', compress = 'xz')
save(apagrd1b, file = 'data/apagrd1b.RData', compress = 'xz')

# gridded comparisons, mean and sd, 7 day -----------------------------------------------------

# this takes about 24 hours to run
grd <- crossing(
  amean = c(0, 2),
  asd = c(0.01, 1),
  rmean = c(0, 200), 
  rsd = c(0.5, 50),
  bmean = c(0, 0.502), 
  bsd = c(0.001, 0.1),
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
  bprior <- c(selrow$bmean, selrow$bsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- ebase(fwdatinp, interval = 900, H = fwdatinp$H, progress = TRUE, n.chains = 4, 
               aprior = aprior, rprior = rprior, bprior = bprior, ndays = ndays)
  
  stopCluster(cl)
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

apagrd <- grd
apagrd7a <- apagrd[1:32,]
apagrd7b <- apagrd[33:64,]
save(apagrd7a, file = 'data/apagrd7a.RData', compress = 'xz')
save(apagrd7b, file = 'data/apagrd7b.RData', compress = 'xz')

# gridded comparisons, mean and sd, 30 day ----------------------------------------------------

# this takes about 48 hours to run
grd <- crossing(
  amean = c(0, 2),
  asd = c(0.01, 1),
  rmean = c(0, 200), 
  rsd = c(0.5, 50),
  bmean = c(0, 0.502), 
  bsd = c(0.001, 0.1),
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
  bprior <- c(selrow$bmean, selrow$bsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  cl <- makeCluster(6)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- ebase(fwdatinp, interval = 900, H = fwdatinp$H, progress = FALSE, n.chains = 4, 
               aprior = aprior, rprior = rprior, bprior = bprior, ndays = ndays)
  
  stopCluster(cl)
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

apagrd <- grd
apagrd30a <- apagrd[1:32,]
apagrd30b <- apagrd[33:64,]
save(apagrd30a, file = 'data/apagrd30a.RData', compress = 'xz')
save(apagrd30b, file = 'data/apagrd30b.RData', compress = 'xz')

# synthetic with EBASE defaults -------------------------------------------

# this takes about 2 hours to run
grd <- crossing(
  ndays = c(1, 7, 30),
  out = NA
)

str <- Sys.time()

for(i in 1:nrow(grd)){
  
  # counter
  cat(i, 'of', nrow(grd), '\n')
  print(Sys.time() - str)
  
  # get inputs
  selrow <- grd[i, ]
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- ebase(fwdatinp, interval = 900, H = fwdatinp$H, progress = FALSE, n.chains = 4, 
               ndays = ndays)
  
  stopCluster(cl)
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

ebasedefault <- grd
save(ebasedefault, file = here('data/ebasedefault.RData'), compress = 'xz')

# evaluate fit all priors -------------------------------------------------

load(file = here('data/apagrd1a.RData'))
load(file = here('data/apagrd1b.RData'))
load(file = here('data/apagrd7a.RData'))
load(file = here('data/apagrd7b.RData'))
load(file = here('data/apagrd30a.RData'))
load(file = here('data/apagrd30b.RData'))

apagrd <- bind_rows(apagrd1a, apagrd1b, apagrd7a, apagrd7b, apagrd30a, apagrd30b)

apasumdat <- apagrd %>% 
  mutate(
    ind = 1:nrow(.),
    ests = purrr::pmap(list(ind, out), function(ind, out){
      
      cat(ind, '\t')

      cmp <- inner_join(fwdatcmp, out[[1]], by = c('Date', 'DateTimeStamp')) %>%
        mutate(
          a.x = a.x / H # need to convert a to m-3, ebase output is m-3, fwoxy is m-2
        ) %>% 
        select(-H, -converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>%
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
        group_by(grp, var) %>% 
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

# gridded comparison, b mean changes ----------------------------------------------------------

# this takes about ten hours to run
grd <- crossing(
  amean = 0.2, #c(0, 0.2, 2),
  asd = 0.1,
  rmean = 20, #c(0, 20, 200), 
  rsd = 5, 
  bmean = c(0.1255, 0.251, 0.502), # 0.251 / 2, 0.251, 0.251 * 2
  bsd = c(0.001, 0.01, 0.1),
  ndays = c(1, 7),
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
  bprior <- c(selrow$bmean, selrow$bsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  cl <- makeCluster(6)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- ebase(fwdatinp, interval = 900, H = fwdatinp$H, progress = TRUE, n.chains = 4, 
               aprior = aprior, rprior = rprior, bprior = bprior, ndays = ndays)
  
  stopCluster(cl)
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

apagrdmean <- grd
save(apagrdmean, file = 'data/apagrdmean.RData', compress = 'xz')

# evaluate fit, b mean changes ----------------------------------------------------------------

load(file = here('data/apagrdmean.RData'))

apasumdatmean <- apagrdmean %>% 
  mutate(
    ind = 1:nrow(.),
    ests = purrr::pmap(list(ind, out), function(ind, out){
      
      cat(ind, '\t')
      
      cmp <- inner_join(fwdatcmp, out[[1]], by = c('Date', 'DateTimeStamp')) %>%
        select(-converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>%
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
      
      sumgrp <- cmp %>% 
        filter(var %in% c('a', 'Rt_vol', 'b')) %>% 
        group_by(grp, var) %>% 
        summarise(
          Fwoxy = mean(Fwoxy, na.rm = T), 
          EBASE = mean(EBASE, na.rm = T), 
          .groups = 'drop'
        ) 
      sumtms <- cmp %>% 
        filter(!var %in% c('a', 'Rt_vol', 'b'))
      
      sumcmp <- bind_rows(sumgrp, sumtms) %>% 
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

save(apasumdatmean, file = here('data/apasumdatmean.RData'), compress = 'xz')

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

apacpdtd <- wtreg(apacp, wins = list(12, 12, 0.4), tz = 'America/Jamaica', lat = 29.7021, long = -84.8802, progress = T)

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
resnos <- ebase(tomod, interval = 900, H = tomod$H, progress = TRUE, n.chains = 4, ndays = 7)

# save ebase results for noisy ts
save(resnos, file = here('data/resnos.RData'))

tomod <- fwdatinpnos %>% 
  select(-obsnoise, -tidnoise, -DO_nos) 

# run model for inputs
cl <- makeCluster(10)
registerDoParallel(cl)

# use interp for missing values
resobs <- ebase(tomod, interval = 900, H = tomod$H, progress = TRUE, n.chains = 4, ndays = 7)

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
