library(tidyverse)
library(lubridate)

# for PAR conversion to W/m2
Jpmolph <-  0.2175e6 # 1 mol-photons = 0.2175e6 J for ave PAR wavelength of 550nm

# format 2012 AP data for use with BASEmetab, observed --------------------

# totpar is mmol/m2/15 min, should be Watts/m2 for EBASE
# sal is ppt, should be ppt
# DO is mg/l, should be mg/l
# temp is C, should be C
# BP is mb, should be atm
# WSpd is m/s, should be m/s
# last line is to fill NA (only 3 values) with last value
APNERR2012 <- read.csv('data/APNERR2012.csv') %>% 
  mutate(
    DateTimeStamp = as.character(DateTimeStamp), 
    totpar = ifelse(is.na(totpar), 0, totpar), 
    totpar = totpar * Jpmolph * 1e-3 / 15/ 60, # convert to W / m2
    BP = BP / 1013 # mb to atm
  ) %>% 
  separate(DateTimeStamp, c('Date', 'Time'), sep = ' ') %>% 
  select(Date, Time, PAR = totpar, tempC = Temp, DO.meas = DO_obs, atmo.pressure = BP, salinity = Sal, WSpd) %>% 
  fill(atmo.pressure, tempC)

# write.csv(APNERR2012, '../ebase/inst/APNERR2012.csv', row.names = F) # cant save here to no overwrite original file
save(APNERR2012, file = 'data/APNERR2012.RData', compress = 'xz')

# format 2012 AP data for use with BASEmetab, optimal detided -------------

# totpar is mmol/m2/15 min, should be Watts/m2 for EBASE
# sal is ppt, should be ppt
# DO is mg/l, should be mg/l
# temp is C, should be C
# BP is mb, should be atm
# WSpd is m/s, should be m/s
# last line is to fill NA (only 3 values) with last value

load(file = 'data/DO_APNERR2012_6_12_0.8.RData')
APNERR2012dtd <- DO_APNERR2012_6_12_0.8 %>% 
  mutate(
    DateTimeStamp = as.character(DateTimeStamp), 
    totpar = ifelse(is.na(totpar), 0, totpar), 
    totpar = totpar * Jpmolph * 1e-3 / 15/ 60, # convert to W / m2
    BP = BP / 1013 # mb to atm
  ) %>% 
  separate(DateTimeStamp, c('Date', 'Time'), sep = ' ') %>% 
  select(Date, Time, PAR = totpar, tempC = Temp, DO.meas = DO_nrm, atmo.pressure = BP, salinity = Sal, WSpd)

write.csv(APNERR2012dtd, file = 'data/APNERR2012dtd.csv', row.names = F)
save(APNERR2012dtd, file = 'data/APNERR2012dtd.RData', compress = 'xz')