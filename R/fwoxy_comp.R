library(fwoxy)
library(tidyverse)
library(lubridate)
# library(WtRegDO)
devtools::load_all('../WtRegDO')

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


tomod <- example %>% 
  mutate(
    DateTimeStamp = force_tz(as.POSIXct(`time, sec`, origin = Sys.Date(), tz = 'UTC'), tzone = 'America/Jamaica'),
    DO_obs = `oxy, mmol/m3` * 0.032, # to mg/L
    Temp = temp_const, 
    WSpd = wspd_const,
    Sal = salt_const, 
    ATemp = NA, 
    BP = 1013.25, 
    Tide = NA
  ) %>% 
  select(DateTimeStamp, Temp, Sal, DO_obs, WSpd, ATemp, BP, Tide)

# currently doesn't work
opmetab <- ecometab(tomod, DO_var = 'DO_obs', tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = 'instant', metab_units = 'grams') 
  

