
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

results.dir <- '~/Desktop/Output/'
data.dir <- results.dir

setwd(data.dir)

# file.copy('C:/proj/BASEmetab_script/data/APNERR2012dtd.csv', '.')

# dat <- read.csv('C:/proj/BASEmetab_script/data/APNERR2012dtd.csv', stringsAsFactors = F) %>% 
#   .[8064:10079, ]
# write.csv(dat, file = 'dat.csv', row.names = F)

#run model
results <- bayesmetab(data.dir, results.dir, interval = 900, K.est = F, K.meas.mean = 0.3865506, K.meas.sd = 1e-6)

# results <- read.csv('~/Desktop/BASE_results_2020-11-15 135525.csv')
toplo <- results %>% 
  # filter(R.Rhat < 1.1) %>%
  # filter(PPP > 0.4) %>%
  mutate(ER.mean = -1 * ER.mean) %>% 
  select(Date, GPP.mean, ER.mean, NEP.mean) %>% 
  mutate(Date = lubridate::ymd(Date)) %>% 
  gather('var', 'val', -Date)

ggplot(toplo, aes(x = Date, y = val, colour = var)) + 
  geom_line() + 
  geom_point()
