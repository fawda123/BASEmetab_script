library(R2jags)
library(foreach)
library(doParallel)
library(here)
library(WtRegDO)
library(EBASE)
# devtools::load_all('../EBASE', helpers = F)
library(tidyverse)

source(here('R/funcs.R'))

load(file = here('data/apacpdtd.RData'))

# depth at site
depth <-  apacpdtd$Tide %>% mean()

# Odum, obs ------------------------------------------------------------------------------------

tmp <- apacpdtd %>% 
  select(DateTimeStamp, Temp, Sal, DO_mgl = DO_obs, ATemp, BP, WSpd, Tide)

opmetab <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.7021, long = -84.8802, gasex = 'Wanninkhof', gasave = 'daily', metab_units = 'mmol', 
                    depth_val = NULL, depth_vec = depth, instant = T)

# all units to mmol O2 m-2 d-1
opmetaball <- opmetab %>% 
  group_by(metab_date) %>% 
  summarize(
    D = mean(D, na.rm = T), 
    P = mean(Pg, na.rm = T), 
    R = -1 * mean(Rt, na.rm = T), 
    NEM = mean(NEM, na.rm = T), 
    .groups = 'drop'
  ) %>% 
  rename(Date = metab_date) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum')

# EBASE, obs ------------------------------------------------------------------------------------

tmp <- apacpdtd %>% 
  select(DateTimeStamp, DO_obs, Temp, Sal, PAR, WSpd)

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# ebase, fix b
res <- ebase(tmp, interval = 900, Z = depth, progress = getwd(), n.chains = 4,
             bprior = c(0.251, 1e-6))

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

apaobscmp <- bind_rows(opmetaball, ebmetaball)

# save(apaobscmp, file = here('data/apaobscmp.RData'))

# Odum, dtd ------------------------------------------------------------------------------------

tmp <- apacpdtd %>% 
  select(DateTimeStamp, Temp, Sal, DO_mgl = DO_nrm, ATemp, BP, WSpd, Tide)

opmetab <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.7021, long = -84.8802, gasex = 'Wanninkhof', gasave = 'daily', metab_units = 'mmol', 
                    depth_val = NULL, depth_vec = depth, instant = T)

# all units to mmol O2 m-2 d-1
opmetaball <- opmetab %>% 
  group_by(metab_date) %>% 
  summarize(
    D = mean(D, na.rm = T), 
    P = mean(Pg, na.rm = T), 
    R = -1 * mean(Rt, na.rm = T), 
    NEM = mean(NEM, na.rm = T), 
    .groups = 'drop'
  ) %>% 
  rename(Date = metab_date) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum')

# EBASE, dtd ------------------------------------------------------------------------------------

tmp <- apacpdtd %>% 
  select(DateTimeStamp, DO_obs = DO_nrm, Temp, Sal, PAR, WSpd)

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# ebase, fix b
res <- ebase(tmp, interval = 900, Z = depth, progress = getwd(), n.chains = 4,
             bprior = c(0.251, 1e-6))

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

apadtdcmp <- bind_rows(opmetaball, ebmetaball)

# save(apadtdcmp, file = here('data/apadtdcmp.RData'))

# combine obs and dtd -------------------------------------------------------------------------

apacmp <- list(
    observed = apaobscmp, 
    detided = apadtdcmp
  ) %>% 
  enframe(name = 'dotyp') %>% 
  unnest('value')

save(apacmp, file = here('data/apacmp.RData'))
