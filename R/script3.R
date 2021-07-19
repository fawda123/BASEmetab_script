# setup -------------------------------------------------------------------

library(R2jags)
library(foreach)
library(doParallel)
library(WtRegDO)
library(here)

source(here('R/funcs.R'))

# number of seconds between observations
interval <- 900

# number of MCMC iterations 
n.iter <- 10000

# run metab_update if T
update.chains <- T

# number of MCMC chains to delete
n.burnin <- n.iter*0.5

# should p be estimated?
p.est <- F

# should theta be estimated?
theta.est <- F

# site depth (m)
depth <- 1.852841

# input dataset
load(file = here('data/APNERR2012dtd.RData'))
assign('data', APNERR2012dtd)
# data <- read.csv('data/Yallakool_example.csv')
# data <- read.csv('output/APNERR2012dtd.csv')

# add DO saturated
data$DO.sat <- dosat_fun(data$tempC, data$salinity, data$atmo.pressure)

# add Kw, wanninkhof is m/d, BASE model has k in d-1, divide by depth at site
data$K <- f_calcWanninkhof(data$tempC, data$salinity, data$WSpd)
data$Kinst <- data$K / depth / (86400 / interval)

# # initial values for inits
# # based on average A and ER from 2012 apa, A is unitless, R is mg O2 per day
# A.init <- 8e-6
# R.init <- 0.65

# Select dates
data$Date <- factor(data$Date, levels = unique(data$Date))
dates <- unique(data$Date)

# evaluate dates with complete record
n.records <- tapply(data$Date, INDEX=data$Date, FUN=length)
dates <- dates[n.records == (86400/interval)] # select only dates with full days

# iterate through each date to estimate metabolism ------------------------

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 3)
registerDoParallel(cl)

# setup log file
strt <- Sys.time()

# process
output <- foreach(d = dates, .packages = 'R2jags', .export = c('interval', 'depth')) %dopar% { 
  
  sink(here('log.txt'))
  cat('Log entry time', as.character(Sys.time()), '\n')
  cat(which(d == dates), ' of ', length(dates), '\n')
  print(Sys.time() - strt)
  sink()
  
  data.sub <- data[data$Date == d,]
  
  # Define data vectors
  num.measurements <- nrow(data.sub)
  tempC <- data.sub$tempC
  salinity <-data.sub$salinity
  atmo.pressure <- data.sub$atmo.pressure
  DO.meas <- data.sub$DO.meas
  PAR <- data.sub$I
  DO.sat <- data.sub$DO.sat
  Kinst <- data.sub$Kinst
  
  # Initial values, leave as NULL if no convergence issues
  inits <- NULL
  # inits <- function(){
  #   list(
  #     A = A.init,
  #     R = R.init / (86400 / interval)
  #   )
  # }
  
  # Different random seeds
  kern=as.integer(runif(1000,min=1,max=10000))
  iters=sample(kern,1)
  
  # Set 
  n.chains <- 3
  n.thin <- 10
  p.est.n <- as.numeric(p.est)
  theta.est.n <- as.numeric(theta.est)
  data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","DO.sat","p.est.n", "theta.est.n", "Kinst", "depth")#, "A.init", "R.init")  
  
  # Define monitoring variables (returned by jags)
  params <- c("Kday", "ER","GPP","NEP")
  
  ## Call jags ##
  metabfit <- do.call(R2jags::jags.parallel, 
                      list(data = data.list, inits = inits, parameters.to.save = params, model.file = here("BASE_metab_model.txt"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                           n.thin = n.thin, n.cluster = n.chains, DIC = TRUE,
                           jags.seed = 123, digits=5)
  )
  
  # update metab if no convergence
  metabfit <- metab_update(metabfit, update.chains, n.iter)
  
  # check final convergence
  srf <- metabfit$BUGSoutput$summary[,8]
  Rhat.test <- ifelse(any(srf > 1.1, na.rm = T) == TRUE, "Check convergence", "Fine")
  
  # insert results to table and write table
  result <- data.frame(Date=as.character(d), 
                       GPP = metabfit$BUGSoutput$mean$GPP, 
                       ER = metabfit$BUGSoutput$mean$ER, 
                       NEP = metabfit$BUGSoutput$mean$NEP,
                       Kday = metabfit$BUGSoutput$mean$Kday,
                       convergence = Rhat.test)
  
  return(result)
  
}

outputkints <- do.call('rbind', output)
save(outputkinst, file = here('data/outputkinst.RData'), compress = 'xz')
