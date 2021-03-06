# setup -------------------------------------------------------------------

library(R2jags)
library(foreach)
library(doParallel)
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

# initial values for inits
K.init <- 2
  
# should k be estimated with uninformative priors?
K.est <- F

# mean for the informed normal prior distribution if K.est = F
# 0.80 is m/d mean wanninkhof for the year, 1.85 is depth at the site, BASE model uses k as d-1
K.meas.mean <-   0.8040253 / 1.852841

# sd for the informed normal prior distribution if K.est = F
K.meas.sd <- 1e-9 

# should p be estimated?
p.est <- FALSE

# should theta be estimated?
theta.est <- FALSE 

# input dataset
load(file = here('data/APNERR2012dtd.RData'))
assign('data', APNERR2012dtd)
# data <- read.csv('data/Yallakool_example.csv')
# data <- read.csv('output/APNERR2012dtd.csv')

# add DO saturated
data$DO.sat <- dosat_fun(data$tempC, data$salinity, data$atmo.pressure)

# Select dates
data$Date <- factor(data$Date, levels = unique(data$Date))
dates <- unique(data$Date)

# evaluate dates with complete record
n.records <- tapply(data$Date, INDEX=data$Date, FUN=length)
dates <- dates[n.records == (86400/interval)] # select only dates with full days

# iterate through each date to estimate metabolism ------------------------

# setup parallel backend
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)

# setup log file
strt <- Sys.time()

# process
output <- foreach(d = dates, .packages = 'R2jags') %dopar% { 
  
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
  
  # Initial values
  inits <- function(){
    list(
      K = K.init / (86400 / interval)
      )
    }
  
  # Different random seeds
  kern=as.integer(runif(1000,min=1,max=10000))
  iters=sample(kern,1)
  
  # Set 
  n.chains <- 3
  n.thin <- 10
  p.est.n <- as.numeric(p.est)
  theta.est.n <- as.numeric(theta.est)
  K.est.n <- as.numeric(K.est)
  K.meas.mean.ts <- K.meas.mean / (86400/interval)
  K.meas.sd.ts <- K.meas.sd / (86400/interval)
  data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","DO.sat", "K.init",
                    "K.est.n", "K.meas.mean.ts", "K.meas.sd.ts", "p.est.n", "theta.est.n")  
  
  # Define monitoring variables
  params <- c("A","R","K","K.day","p","theta","tau","ER","GPP","NEP","PR","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled",
                         "gppts", "erpts", "kpts")
  
  ## Call jags ##
  metabfit <- do.call(R2jags::jags.parallel, 
                      list(data = data.list, inits = inits, parameters.to.save = params, model.file = here("BASE_metab_model_v2.3.txt"),
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
                       K = metabfit$BUGSoutput$mean$K.day, 
                       convergence = Rhat.test)
  
  return(result)
  
}

outputpar <- do.call('rbind', output)
save(outputpar, file = here('data/outputpar.RData', compress = 'xz'))
