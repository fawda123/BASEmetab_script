# setup -------------------------------------------------------------------

library(R2jags)
library(foreach)
library(doParallel)
library(WtRegDO)
library(dplyr)
library(tidyr)
library(lubridate)
library(fwoxy)
library(here)

source(here('R/funcs.R'))

# Set model parameters
oxy_ic <- 250           # (mmol/m^3), initial oxygen concentration
a_param <- 0.2          # ((mmol/m^3)/day)/(W/m^2), light efficiency
er_param <- 20          # (mmol/m^3/day), ecosystem respiration

# Constant Forcings
ht_const <- 3           # m, height of the water column
salt_const <- 25        # ppt, salinity
temp_const <- 25        # deg C, water temperature
wspd_const <- 3         # m/s, wind speed at 10 m

example <- fwoxy(oxy_ic = oxy_ic, a_param = a_param, er_param = er_param, 
                 ht_in = ht_const, salt_in = salt_const, temp_in = temp_const,
                 wspd_in = wspd_const)

# number of seconds between observations
interval <- 900

# number of MCMC iterations 
n.iter <- 10000

# run metab_update if T
update.chains <- T

# number of MCMC chains to delete
n.burnin <- n.iter*0.5

# should k be estimated with uninformative priors?
K.est <- F

# mean for the informed normal prior distribution if K.est = F
# 0.67 is m/d from Fwoxy, 3 is depth at the site, BASE model uses k as d-1
K.meas.mean <-   0.6702159 / ht_const

# sd for the informed normal prior distribution if K.est = F
# this sd is the same range as a proportion of the mean for ebase b parameter
K.meas.sd <- 0.00890061

# should p be estimated?
p.est <- FALSE

# should theta be estimated?
theta.est <- FALSE 

# site depth (m)
depth <- ht_const

# prep fwoxy output for BASEmetab
# DO.meas is mmol/m3, should be mg/l
# tempC is C, should be C
# WSpd is m/s, should be m/s
# salinity is ppt, should be ppt
# I (par) is W/m2, should be umol/m2/s (multiply by 4.57 https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf)
# atmo.pressure should be atm, 1 at sea level
data <- example %>% 
  mutate(
    DateTimeStamp = force_tz(as.POSIXct(`time, sec`, origin = Sys.Date(), tz = 'UTC'), tzone = 'America/Jamaica'),
    DateTimeStamp = as.character(DateTimeStamp),
    DO.meas = `oxy, mmol/m3` * 0.032, # to mg/L
    tempC = temp_const, 
    WSpd = wspd_const,
    salinity = salt_const, 
    I = fun_par_sin_model(`time, sec`) * 4.57,
    atmo.pressure = 1 
  ) %>%
  separate(DateTimeStamp, c('Date', 'Time'), sep = ' ') %>% 
  select(Date, Time, I, tempC, DO.meas, atmo.pressure, salinity, WSpd)  

# add DO saturated
data$DO.sat <- dosat_fun(data$tempC, data$salinity, data$atmo.pressure)

# Select dates
data$Date <- factor(data$Date, levels = unique(data$Date))
dates <- unique(data$Date)

# evaluate dates with complete record
n.records <- tapply(data$Date, INDEX=data$Date, FUN=length)
dates <- dates[n.records == (86400/interval)] # select only dates with full days

# iterate through each date to estimate metabolism ------------------------

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# setup log file
strt <- Sys.time()

# process
output <- foreach(d = dates, .packages = c('here', 'R2jags'), .export = c('interval')) %dopar% { 
  
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
  
  # Initial values, leave as NULL if no convergence issues
  inits <- function(){
    list(
      K = K.meas.mean.ts
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
  data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","DO.sat",
                    "K.est.n", "K.meas.mean.ts", "K.meas.sd.ts", "p.est.n", "theta.est.n")
  
  # Define monitoring variables (returned by jags)
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
  # insert results to table and write table
  result <- data.frame(Date=as.character(d), 
                       Time = data.sub$Time,
                       gppts = c(NA, metabfit$BUGSoutput$mean$gppts), # O2, g/m3/15 min
                       erpts = c(NA, metabfit$BUGSoutput$mean$erpts), # O2, g/m3/15 min 
                       gets = c(NA, metabfit$BUGSoutput$mean$gets), # O2, g/m3/15 min
                       dDO = c(NA, diff(metabfit$BUGSoutput$mean$DO.modelled)) # O2 g/m3/15 min
  )
  
  return(result)
  
}

stopCluster(cl)

# correct instantanous obs to daily, g to mmol
fwoxybasemetab <- do.call('rbind', output) %>% 
  na.omit() %>% 
  unite(DateTimeStamp, c('Date', 'Time'), sep = '_') %>% 
  mutate(
    DateTimeStamp = lubridate::ymd_hms(DateTimeStamp, tz = 'America/Jamaica'),
    Pg_vol = gppts * (86400 / interval) / 0.032, # O2 mmol/m3/d
    Rt_vol = erpts * (86400 / interval) / 0.032, # O2 mmol/m3/d
    D = -1 * gets * (86400 / interval) / 0.032, # O2 mmol/m3/d
    dDO = dDO * (86400 / interval) / 0.032 # O2 mmol/m3/d
  ) %>% 
  select(DateTimeStamp, Pg_vol, Rt_vol, D, dDO)
save(fwoxybasemetab, file = here('data/fwoxybasemetab.RData'), compress = 'xz')
