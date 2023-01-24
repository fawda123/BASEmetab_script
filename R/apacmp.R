library(R2jags)
library(foreach)
library(doParallel)
library(here)
library(WtRegDO)
library(EBASE)
library(tidyverse)

source(here('R/funcs.R'))

# depth at site
depth <- 1.852841

# Odum, obs ------------------------------------------------------------------------------------

load(file = here('data/DO_APNERR2012_6_12_0.8.RData'))
tmp <- DO_APNERR2012_6_12_0.8 %>%
  select(DateTimeStamp, DO_mgl = DO_obs, Temp, Sal, ATemp, BP, WSpd, totpar, Tide)

opmetab <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = 'daily', metab_units = 'mmol', 
                    depth_val = NULL, depth_vec = depth, instant = T)

# all units to mmol O2 m-2 d-1
opmetaball <- opmetab %>% 
  group_by(metab_date) %>% 
  summarize(
    D = mean(depth * D, na.rm = T), 
    P = mean(Pg, na.rm = T), 
    R = -1 * mean(Rt, na.rm = T), 
    NEM = mean(NEM, na.rm = T), 
    .groups = 'drop'
  ) %>% 
  rename(Date = metab_date) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum')

# BASEmetab, obs -------------------------------------------------------------------------------

# number of seconds between observations
interval <- 900

# number of MCMC iterations 
n.iter <- 10000

# run metab_update if T
update.chains <- T

# number of MCMC chains to delete
n.burnin <- n.iter*0.5

# initial values for inits
K.init <- 0.7993042 / 1.852841

# should k be estimated with uninformative priors?
K.est <- F

# mean for the informed normal prior distribution if K.est = F
# 0.80 is m/d mean wanninkhof for the year from odum above
# 1.85 is depth at the site, BASE model uses k as d-1
K.meas.mean <-   0.7993042 / 1.852841

# sd for the informed normal prior distribution if K.est = F, set as (0.01 / .2501) * (0.8040253 / 1.852841) as prop from ebase
K.meas.sd <- (0.01 / .2501) * (0.7993042 / 1.852841)

# should p be estimated?
p.est <- FALSE

# should theta be estimated?
theta.est <- FALSE 

# input dataset
load(file = here('data/APNERR2012.RData'))
assign('data', APNERR2012)

# add DO saturated
data$DO.sat <- dosat_fun(data$tempC, data$salinity, data$atmo.pressure)

# Select dates
data$Date <- factor(data$Date, levels = unique(data$Date))
dates <- unique(data$Date)

# evaluate dates with complete record
n.records <- tapply(data$Date, INDEX=data$Date, FUN=length)
dates <- dates[n.records == (86400/interval)] # select only dates with full days

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 3)
registerDoParallel(cl)

# setup log file
strt <- Sys.time()

# process
output <- foreach(d = dates, .packages = c('here', 'R2jags')) %dopar% { 
  
  sink(here('log.txt'))
  cat('Log entry time', as.character(Sys.time()), '\n')
  cat(which(d == dates), ' of ', length(dates), '\n')
  print(Sys.time() - strt)
  sink()
  
  data.sub <- data[data$Date == d,]
  
  # Define data vectors
  num.measurements <- nrow(data.sub)
  tempC <- data.sub$tempC
  salinity <-data.sub$salinity
  atmo.pressure <- data.sub$atmo.pressure
  DO.meas <- data.sub$DO.meas
  PAR <- data.sub$I
  DO.sat <- data.sub$DO.sat
  
  # Initial values
  inits <- function(){
    list(
      K = K.init / (86400 / interval)
    )
  }
  
  # Different random seeds
  kern=as.integer(runif(1000,min=1,max=10000))
  iters=sample(kern,1)
  
  # Set 
  n.chains <- 3
  n.thin <- 10
  p.est.n <- as.numeric(p.est)
  theta.est.n <- as.numeric(theta.est)
  K.est.n <- as.numeric(K.est)
  K.meas.mean.ts <- K.meas.mean / (86400/interval)
  K.meas.sd.ts <- K.meas.sd / (86400/interval)
  data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","DO.sat", "K.init",
                    "K.est.n", "K.meas.mean.ts", "K.meas.sd.ts", "p.est.n", "theta.est.n")  
  
  # Define monitoring variables
  params <- c("A","R","K","K.day","p","theta","tau","ER","GPP","NEP","PR","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled",
              "gppts", "erpts", "kpts", "gets")
  
  ## Call jags ##
  metabfit <- do.call(R2jags::jags.parallel, 
                      list(data = data.list, inits = inits, parameters.to.save = params, model.file = here("BASE_metab_model_v2.3.txt"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                           n.thin = n.thin, n.cluster = n.chains, DIC = TRUE,
                           jags.seed = 123, digits=5)
  )
  
  # update metab if no convergence
  metabfit <- metab_update(metabfit, update.chains, n.iter)
  
  # check final convergence
  srf <- metabfit$BUGSoutput$summary[,8]
  Rhat.test <- ifelse(any(srf > 1.1, na.rm = T) == TRUE, "Check convergence", "Fine")
  
  # insert results to table and write table
  result <- data.frame(Date=as.character(d), 
                       GPP = metabfit$BUGSoutput$mean$GPP, 
                       ER = metabfit$BUGSoutput$mean$ER, 
                       NEP = metabfit$BUGSoutput$mean$NEP,
                       D = mean(metabfit$BUGSoutput$mean$gets * (86400 / interval))
                       )
  
  return(result)
  
}

# convert form g O2 m-3 d-1 to mmol O2 m-2 d-1
bsmetaball <- do.call('rbind', output) %>% 
  mutate(
    Date = lubridate::mdy(Date), 
    P = depth * GPP / 0.032,
    R = depth * ER / 0.032,
    NEM = depth * NEP / 0.032, 
    D = depth * D / 0.032
  ) %>% 
  select(Date, P, R, NEM, D) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'BASEmetab')


# EBASE, obs ------------------------------------------------------------------------------------

# setup parallel backend
cl <- makeCluster(5)
registerDoParallel(cl)

res <- ebase(exdat, interval = 900, H = 1.852841, progress = TRUE, n.chains = 5)

stopCluster(cl)

# all units mmol O2 m-2 d-1
ebmetaball <- res %>% 
  mutate(
    Date = as.Date(DateTimeStamp)
  ) %>% 
  group_by(Date) %>% 
  summarise(
    P = mean(P, na.rm = T), 
    R = mean(R, na.rm = T), 
    D = mean(D, na.rm = T),
    .groups = 'drop'
  ) %>% 
  mutate(
    NEM = P - R
  ) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'EBASE')

# combineall, obs -------------------------------------------------------------------------------

apaobscmp <- bind_rows(opmetaball, bsmetaball, ebmetaball)

# save(apaobscmp, file = here('data/apaobscmp.RData'))

# Odum, dtd ------------------------------------------------------------------------------------

load(file = here('data/DO_APNERR2012_6_12_0.8.RData'))
tmp <- DO_APNERR2012_6_12_0.8 %>%
  select(DateTimeStamp, DO_mgl = DO_nrm, Temp, Sal, ATemp, BP, WSpd, totpar, Tide)

opmetab <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = 'daily', metab_units = 'mmol', 
                    depth_val = NULL, depth_vec = depth, instant = T)

# all units to mmol O2 m-2 d-1
opmetaball <- opmetab %>% 
  group_by(metab_date) %>% 
  summarize(
    D = mean(depth * D, na.rm = T), 
    P = mean(Pg, na.rm = T), 
    R = -1 * mean(Rt, na.rm = T), 
    NEM = mean(NEM, na.rm = T), 
    .groups = 'drop'
  ) %>% 
  rename(Date = metab_date) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum')


# BASEmetab, dtd -------------------------------------------------------------------------------

# number of seconds between observations
interval <- 900

# number of MCMC iterations 
n.iter <- 10000

# run metab_update if T
update.chains <- T

# number of MCMC chains to delete
n.burnin <- n.iter*0.5

# initial values for inits
K.init <- 0.7993042 / 1.852841

# should k be estimated with uninformative priors?
K.est <- F

# mean for the informed normal prior distribution if K.est = F
# 0.80 is m/d mean wanninkhof for the year from odum above
# 1.85 is depth at the site, BASE model uses k as d-1
K.meas.mean <-   0.7993042 / 1.852841

# sd for the informed normal prior distribution if K.est = F, set as (0.01 / .2501) * (0.8040253 / 1.852841) as prop from ebase
K.meas.sd <- (0.01 / .2501) * (0.7993042 / 1.852841)

# should p be estimated?
p.est <- FALSE

# should theta be estimated?
theta.est <- FALSE 

# input dataset
load(file = here('data/APNERR2012dtd.RData'))
assign('data', APNERR2012dtd)

# add DO saturated
data$DO.sat <- dosat_fun(data$tempC, data$salinity, data$atmo.pressure)

# Select dates
data$Date <- factor(data$Date, levels = unique(data$Date))
dates <- unique(data$Date)

# evaluate dates with complete record
n.records <- tapply(data$Date, INDEX=data$Date, FUN=length)
dates <- dates[n.records == (86400/interval)] # select only dates with full days

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 3)
registerDoParallel(cl)

# setup log file
strt <- Sys.time()

# process
output <- foreach(d = dates, .packages = c('here', 'R2jags')) %dopar% { 
  
  sink(here('log.txt'))
  cat('Log entry time', as.character(Sys.time()), '\n')
  cat(which(d == dates), ' of ', length(dates), '\n')
  print(Sys.time() - strt)
  sink()
  
  data.sub <- data[data$Date == d,]
  
  # Define data vectors
  num.measurements <- nrow(data.sub)
  tempC <- data.sub$tempC
  salinity <-data.sub$salinity
  atmo.pressure <- data.sub$atmo.pressure
  DO.meas <- data.sub$DO.meas
  PAR <- data.sub$I
  DO.sat <- data.sub$DO.sat
  
  # Initial values
  inits <- function(){
    list(
      K = K.init / (86400 / interval)
    )
  }
  
  # Different random seeds
  kern=as.integer(runif(1000,min=1,max=10000))
  iters=sample(kern,1)
  
  # Set 
  n.chains <- 3
  n.thin <- 10
  p.est.n <- as.numeric(p.est)
  theta.est.n <- as.numeric(theta.est)
  K.est.n <- as.numeric(K.est)
  K.meas.mean.ts <- K.meas.mean / (86400/interval)
  K.meas.sd.ts <- K.meas.sd / (86400/interval)
  data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","DO.sat", "K.init",
                    "K.est.n", "K.meas.mean.ts", "K.meas.sd.ts", "p.est.n", "theta.est.n")  
  
  # Define monitoring variables
  params <- c("A","R","K","K.day","p","theta","tau","ER","GPP","NEP","PR","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled",
              "gppts", "erpts", "kpts", "gets")
  
  ## Call jags ##
  metabfit <- do.call(R2jags::jags.parallel, 
                      list(data = data.list, inits = inits, parameters.to.save = params, model.file = here("BASE_metab_model_v2.3.txt"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                           n.thin = n.thin, n.cluster = n.chains, DIC = TRUE,
                           jags.seed = 123, digits=5)
  )
  
  # update metab if no convergence
  metabfit <- metab_update(metabfit, update.chains, n.iter)
  
  # check final convergence
  srf <- metabfit$BUGSoutput$summary[,8]
  Rhat.test <- ifelse(any(srf > 1.1, na.rm = T) == TRUE, "Check convergence", "Fine")
  
  # insert results to table and write table
  result <- data.frame(Date=as.character(d), 
                       GPP = metabfit$BUGSoutput$mean$GPP, 
                       ER = metabfit$BUGSoutput$mean$ER, 
                       NEP = metabfit$BUGSoutput$mean$NEP,
                       D = mean(metabfit$BUGSoutput$mean$gets * (86400 / interval))
  )
  
  return(result)
  
}

# convert form g O2 m-3 d-1 to mmol O2 m-2 d-1
bsmetaball <- do.call('rbind', output) %>% 
  mutate(
    Date = lubridate::ymd(Date), 
    P = depth * GPP / 0.032,
    R = depth * ER / 0.032,
    NEM = depth * NEP / 0.032, 
    D = depth * D / 0.032
  ) %>% 
  select(Date, P, R, NEM, D) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'BASEmetab')

# EBASE, dtd ------------------------------------------------------------------------------------

load(file = here('data/APNERR2012dtd.RData'))
assign('data', APNERR2012dtd)
exdatdtd <- APNERR2012dtd %>% 
  rename(
    PAR = I,
    Sal = salinity,
    Temp = tempC, 
    DO_obs = DO.meas
  ) %>% 
  unite('DateTimeStamp', Date, Time, sep = ' ') %>% 
  mutate(
    DateTimeStamp = lubridate::ymd_hms(DateTimeStamp, tz = 'America/Jamaica')
  ) %>% 
  select(DateTimeStamp, DO_obs, Temp, Sal, PAR, WSpd) 

# setup parallel backend
cl <- makeCluster(5)
registerDoParallel(cl)

res <- ebase(exdatdtd, interval = 900, H = 1.852841, progress = TRUE, n.chains = 5)

stopCluster(cl)

# all units mmol O2 m-2 d-1
ebmetaball <- res %>% 
  mutate(
    Date = as.Date(DateTimeStamp)
  ) %>% 
  group_by(Date) %>% 
  summarise(
    P = mean(P, na.rm = T), 
    R = mean(R, na.rm = T), 
    D = mean(D, na.rm = T),
    .groups = 'drop'
  ) %>% 
  mutate(
    NEM = P - R
  ) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'EBASE')

# combineall, dtd ------------------------------------------------------------------------------

apadtdcmp <- bind_rows(opmetaball, bsmetaball, ebmetaball)

# save(apadtdcmp, file = here('data/apadtdcmp.RData'))

# combine obs and dtd -------------------------------------------------------------------------

apacmp <- list(
    observed = apaobscmp, 
    detided = apadtdcmp
  ) %>% 
  enframe(name = 'dotyp') %>% 
  unnest('value')

save(apacmp, file = here('data/apacmp.RData'))
