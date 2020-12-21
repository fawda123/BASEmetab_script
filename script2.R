
# Running original BASEmetab model ----------------------------------------

library(BASEmetab)
library(tidyverse)
library(R2jags)
library(foreach)
library(doParallel)

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 1)
registerDoParallel(cl)

results.dir <- 'output'
data.dir <- results.dir

# file.copy('C:/proj/BASEmetab_script/data/APNERR2012dtd.csv', '.')

# dat <- read.csv('C:/proj/BASEmetab_script/data/APNERR2012dtd.csv', stringsAsFactors = F) %>% 
#   .[8064:10079, ]
# write.csv(dat, file = 'dat.csv', row.names = F)

#run model

# cat point bottom depth (mean of tidal vector plus 0.3m)
H <- 1.852841

# areal K
Kareal <- 0.7993042

# volumetric
kvol <- Kareal / H

results <- bayesmetab(data.dir, results.dir, interval = 900, K.est = F, K.meas.mean = kvol, K.meas.sd = 1e-9, instant = F, update.chains = F)

# remove image files from output
file.remove(list.files(path = results.dir, pattern = '\\.jpg$', full.names = T))


