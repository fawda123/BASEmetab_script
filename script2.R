
# Running original BASEmetab model ----------------------------------------

library(BASEmetab)
library(tidyverse)
library(R2jags)
library(foreach)
library(doParallel)

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

results.dir <- 'output'
data.dir <- results.dir

# 2012 four days
load(file = 'data/APNERR2012dtd.RData')
APNERR2012dtd <- APNERR2012dtd[6528:6911, ]
write.csv(APNERR2012dtd, 'output/APNERR2012dtd.csv', row.names = F)

# # 2012 all
# load(file = 'data/APNERR2012dtd.RData')
# assign('data', APNERR2012dtd)
# write.csv(APNERR2012dtd, 'output/APNERR2012dtd.csv', row.names = F)

#run model

# cat point bottom depth (mean of tidal vector plus 0.3m)
H <- 1.852841

# areal K, 0.80 is m/d mean wanninkhof for the year at apa 2012
Kareal <- 0.8040253

# volumetric
kvol <- Kareal / H

results <- bayesmetab(data.dir, results.dir, interval = 900, K.est = F, K.meas.mean = kvol, K.meas.sd = 1e-9, instant = F, update.chains = T)

# remove image files from output
file.remove(list.files(path = results.dir, pattern = '\\.jpg$', full.names = T))


