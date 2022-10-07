# setup -------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(fwoxy)
library(here)
library(EBASE)

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
                 wspd_in = wspd_const, plot = F)

data <- example %>% 
  mutate(
    DateTimeStamp = force_tz(as.POSIXct(`time, sec`, origin = Sys.Date(), tz = 'UTC'), tzone = 'America/Jamaica'),
    DO_obs = `oxy, mmol/m3` * 32 / 1000,
    Temp = temp_const, 
    Sal = salt_const, 
    PAR = fun_par_sin_model(`time, sec`),
    WSpd = wspd_const
  ) %>%
  select(DateTimeStamp, DO_obs, Temp, Sal, PAR, WSpd) %>% 
  .[-577,]

res <- ebase(data, interval = 900, ndays = 6, H = ht_const, progress = F, n.chains = 5, rprior = c(20, 5), bprior = c(0.251, 0.01), aprior = c(0.2, 0.1))

fwoxyebase <- res %>% 
  select(
    DateTimeStamp, 
    DO.meas = DO_obs, 
    DO.modelled = DO_mod, 
    dDO, 
    a, 
    b, 
    Pg_vol, 
    Rt_vol, 
    D
  )

save(fwoxyebase, file = here('data/fwoxyebase.RData'))
