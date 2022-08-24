library(EBASE)
library(tidyverse)
library(lubridate)
library(doParallel)
library(here)

H <- 1.85
tz <- 'America/Jamaica'
dtrng1 <- ymd('2012-03-01', '2012-03-14')
dtrng2 <- ymd('2012-08-01', '2012-08-14')


# data prep -----------------------------------------------------------------------------------

load(file = 'data/DO_APNERR2012_6_12_0.8.RData')

# totpar is mmol/m2 total for 15 minute obs, umol/m2/s
# sal is ppt, should be ppt
# DO is mg/l, should be mg/l
# temp is C, should be C
# BP is mb, should be atm
# WSpd is m/s, should be m/s
# last line is to fill NA (only 3 values) with last value
ebdat <- DO_APNERR2012_6_12_0.8 %>% 
  mutate(
    totpar = ifelse(is.na(totpar), 0, totpar), 
    totpar = totpar * 1000 / (15 * 60), # convert to umol and per second
  ) %>% 
  select(DateTimeStamp, DO_obs, Temp, Sal, PAR = totpar, WSpd) 


# winter --------------------------------------------------------------------------------------

ebdat1 <- ebdat %>% 
  filter(as.Date(DateTimeStamp, tz = tz) >= dtrng1[1] & as.Date(DateTimeStamp, tz = tz) <= dtrng1[2])

stps <- 14

out1 <- vector('list', stps)
for(i in 1:stps){
  
  cat(i, '\t')
  H <- 1.85
  ebmetab <- ebase(ebdat1, H = H, interval = 900, ndays = i, n.chains = 5, progress = F)
  out1[[i]] <- ebmetab
  
}

ebasendays1 <- out1
save(ebasendays1, file = 'data/ebasendays1.RDAta')

# summer --------------------------------------------------------------------------------------

ebdat2 <- ebdat %>% 
  filter(as.Date(DateTimeStamp, tz = tz) >= dtrng2[1] & as.Date(DateTimeStamp, tz = tz) <= dtrng2[2])

stps <- 14
 
out2 <- vector('list', stps)
for(i in 1:stps){
  
  cat(i, '\t')
  H <- 1.85
  ebmetab <- ebase(ebdat2, H = H, interval = 900, ndays = i, n.chains = 5, progress = F)
  out2[[i]] <- ebmetab
  
}

ebasendays2 <- out2
save(ebasendays2, file = 'data/ebasendays2.RDAta')
