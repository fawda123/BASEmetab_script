library(R2jags)
library(foreach)
library(doParallel)
library(here)
library(WtRegDO)
# library(EBASE)
devtools::load_all('../EBASE', helpers = F)
library(tidyverse)

source(here('R/funcs.R'))

load(file = here('data/apacpdtd.RData'))

# depth at site
depth <-  apacpdtd$Tide %>% mean()

tmp <- apacpdtd %>% 
  select(DateTimeStamp, DO_obs, Temp, Sal, PAR, WSpd, Tide)

# 1 day, fixed z ------------------------------------------------------------------------------

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# ebase, fix b
res <- ebase(tmp, interval = 900, Z = depth, progress = NULL, n.chains = 4,
             bprior = c(0.251, 1e-6))

stopCluster(cl)

apaobscmp1day <- res

save(apaobscmp1day, file = here('data/apaobscmp1day.RData'))

# 1 day, variable z ---------------------------------------------------------------------------

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# ebase, fix b
res <- ebase(tmp, interval = 900, Z = tmp$Tide, progress = NULL, n.chains = 4,
             bprior = c(0.251, 1e-6))

stopCluster(cl)

apaobscmp1dayvarz <- res

save(apaobscmp1dayvarz, file = here('data/apaobscmp1dayvarz.RData'))

# 7 days, fixed z -----------------------------------------------------------------------------

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# ebase, fix b
res <- ebase(tmp, interval = 900, Z = depth, ndays = 7, progress = NULL, n.chains = 4,
             bprior = c(0.251, 1e-6))

stopCluster(cl)

apaobscmp7day <- res

save(apaobscmp7day, file = here('data/apaobscmp7day.RData'))

# 7 days, variable z --------------------------------------------------------------------------

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# ebase, fix b
res <- ebase(tmp, interval = 900, Z = tmp$Tide, ndays = 7, progress = NULL, n.chains = 4,
             bprior = c(0.251, 1e-6))

stopCluster(cl)

apaobscmp7dayvarz <- res

save(apaobscmp7dayvarz, file = here('data/apaobscmp7dayvarz.RData'))

