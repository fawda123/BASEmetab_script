# setup -------------------------------------------------------------------

library(R2jags)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(fwoxy)
library(here)

source(here('R/funcs.R'))

# number of seconds between observations
interval <- 900
nstepd <- 86400 / interval

# number of MCMC iterations 
n.iter <- 10000

# run metab_update if T
update.chains <- T

# number of MCMC chains to delete
n.burnin <- n.iter*0.5

# prep fwoxy output for BASEmetab
# DO.meas is mmol/m3
# 
# tempC is C, should be C
# WSpd is m/s, should be m/s
# salinity is ppt, should be ppt
# PAR is W/m2, should be umol/m2/s (multiply by 4.57 https://www.controlledenvironments.org/wp-content/uploads/sites/6/2017/06/Ch01.pdf)
# sc is unitless
# tempC is C
# DO.meas is mmol/m3
# DO.sat is mmol/m3
# U10 is m/s
load(file = here('data/APNERR2012dtd.RData'))
assign('data', APNERR2012dtd)
data <- data %>% 
  rename(
    PAR = I,
    U10 = WSpd
  ) %>% 
  mutate(
    DO.meas = DO.meas / 32 * 1000,
    DO.sat = fun_eqb_oxygen(tempC, salinity), 
    H = 1.852841, 
    sc = fun_schmidt_oxygen(tempC, salinity)
  ) %>% 
  select(Date, Time, H, PAR, tempC, sc, DO.meas, DO.sat, U10)

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
output <- foreach(d = dates, .packages = c('here', 'R2jags'), .export = c('nstepd')) %dopar% { 
  
  sink(here('log.txt'))
  cat('Log entry time', as.character(Sys.time()), '\n')
  cat(which(d == dates), ' of ', length(dates), '\n')
  print(Sys.time() - strt)
  sink()
  
  data.sub <- data[data$Date == d,]
  
  # Define data vectors
  num.measurements <- nrow(data.sub)
  DO.meas <- data.sub$DO.meas
  PAR <- data.sub$PAR
  DO.sat <- data.sub$DO.sat
  sc <- data.sub$sc
  H <- data.sub$H
  U10 <- data.sub$U10
  
  # Initial values, leave as NULL if no convergence issues
  inits <- NULL
  # inits <- function(){
  #   list(
  #     a = 0.2 / nstepd,
  #     r = 20 / nstepd,
  #     b = 0.251 / 400
  #   )
  
  # Different random seeds
  kern <- as.integer(runif(1000,min=1,max=10000))
  iters <- sample(kern,1)
  
  # Set 
  n.chains <- 3
  n.thin <- 10
  data.list <- list("num.measurements", "nstepd", "DO.meas", "PAR", "DO.sat", "sc", "H", "U10")
  
  # Define monitoring variables (returned by jags)
  params <- c("ats", "bts", "gppts", "erts", "gets", "DO.modelled")
  
  ## Call jags ##
  metabfit <- do.call(R2jags::jags.parallel, 
                      list(data = data.list, inits = inits, parameters.to.save = params, model.file = here::here("ebase_model.txt"),
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
                       DO.meas = data.sub$DO.meas,
                       DO.modelled = metabfit$BUGSoutput$mean$DO.modelled,
                       Time = data.sub$Time,
                       ats = c(NA, metabfit$BUGSoutput$mean$ats), # (mmol/m3/ts)/(W/m2)
                       bts = c(NA, metabfit$BUGSoutput$mean$bts), # ts/m
                       gppts = c(NA, metabfit$BUGSoutput$mean$gppts), # O2, mmol/m3/ts
                       erts = c(NA, metabfit$BUGSoutput$mean$erts), # O2, mmol/m3/ts 
                       gets = c(NA, metabfit$BUGSoutput$mean$gets), # O2, mmol/m3/ts
                       dDO = c(NA, diff(metabfit$BUGSoutput$mean$DO.modelled)) # O2 mmol/m3/ts
                       
  )
  
  return(result)
  
}

stopCluster(cl)

# correct instantanous obs to daily rates
resall <- do.call('rbind', output) %>% 
  na.omit() %>% 
  unite(DateTimeStamp, c('Date', 'Time'), sep = '_') %>% 
  mutate(
    DateTimeStamp = lubridate::ymd_hms(DateTimeStamp, tz = 'America/Jamaica'),
    a = ats * nstepd, # (mmol/m3/ts)/(W/m2) to (mmol/m3/d)/(W/m2)
    b = bts * 100 * 3600 / interval, # (m/d)/(m2/s2) to (cm/hr)/(m2/s2)
    Pg_vol = gppts * nstepd, # O2 mmol/m3/ts to O2 mmol/m3/d
    Rt_vol = erts * nstepd, # O2 mmol/m3/ts to O2 mmol/m3/d
    D = -1 * gets * nstepd, #  # O2 mmol/m3/ts to O2 mmol/m3/d
    dDO = dDO * nstepd #  # O2 mmol/m3/ts to O2 mmol/m3/d
  ) %>% 
  select(-ats, -bts, -gppts, -erts, -gets)

# summarize daily estimates, get NEM, convert to O2 to g/m3/d
outputebase <- resall %>% 
  mutate(
    Date = as.Date(DateTimeStamp)
  ) %>% 
  group_by(Date) %>% 
  summarise(
    GPP = mean(Pg_vol), 
    ER = mean(Rt_vol), 
    .groups = 'drop'
  ) %>% 
  mutate(
    NEP = GPP - ER,
    NEP = NEP * 0.032, 
    GPP = GPP * 0.032, 
    ER = ER * 0.032
  )

save(outputebase, file = here('data/outputebase.RData'))
