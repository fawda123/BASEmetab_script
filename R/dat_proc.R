library(tidyverse)
library(lubridate)

# format 2012 AP data for use with BASEmetab, observed --------------------

# totpar is mmol/m2 total for 15 minute obs, umol/m2/s
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
    totpar = totpar * 1000 / (15 * 60), # convert to umol and per second
    BP = BP / 1013 # mb to atm
  ) %>% 
  separate(DateTimeStamp, c('Date', 'Time'), sep = ' ') %>% 
  select(Date, Time, I = totpar, tempC = ATemp, DO.meas = DO_obs, atmo.pressure = BP, salinity = Sal, WSpd) %>% 
  fill(atmo.pressure, tempC)

save(APNERR2012, file = 'data/APNERR2012.RData', compress = 'xz')


# format 2012 AP data for use with BASEmetab, optimal detided -------------

# totpar is mmol/m2 total for 15 minute obs, umol/m2/s
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
    totpar = totpar * 1000 / (15 * 60), # convert to umol and per second
    BP = BP / 1013 # mb to atm
  ) %>% 
  separate(DateTimeStamp, c('Date', 'Time'), sep = ' ') %>% 
  select(Date, Time, I = totpar, tempC = Temp, DO.meas = DO_nrm, atmo.pressure = BP, salinity = Sal, WSpd)

write.csv(APNERR2012dtd, file = 'data/APNERR2012dtd.csv', row.names = F)
save(APNERR2012dtd, file = 'data/APNERR2012dtd.RData', compress = 'xz')

# apa observed with cordat floored sun angle ------------------------------

# original apaobs created in appalachicola.Rmd

box::use(
  WtRegDO[ecometab, wtreg, meteval, met_day_fun, evalcor],
  here[...], 
  doParallel[registerDoParallel], 
  parallel[detectCores],
  dplyr[...],
  tidyr[...],
  SWMPr[import_local, qaqc, comb], 
  ggplot2[...], 
  tibble[enframe], 
  lubridate[...], 
  ggforce[facet_zoom],
  oce[...],
  foreach[...]
)

source('R/funcs.R')

data(apaobs)

apaobs <- apaobs %>%
  filter(year(DateTimeStamp) == 2012) %>% 
  select(-cordat)

# site metadata
locs <- SWMPr::stat_locs %>% 
  filter(station_code == 'apacp')
lat <- locs$latitude
long <- locs$longitude
tz <- attr(apaobs$DateTimeStamp, 'tzone')

# setup parallel backend
ncores <- detectCores()  
registerDoParallel(cores = ncores - 1)

# get evalcor results to check sun angle/tidal height correlations
cordatflr <- evalcorflr(apaobs, tz, lat, long, progress = T, plot = F)
apaobs$cordatflr <- cordatflr

# get regular evalcor results
cordat <- evalcor(apaobs, tz, lat, long, progress = T, plot = F)
apaobs$cordat <- cordat
apaevlcr <- apaobs

save(apaevlcr, file = here('data/apaevlcr.RData'))

