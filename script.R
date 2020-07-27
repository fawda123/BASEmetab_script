#' bayesmetab
#'
#' Estimates single-station whole-stream metabolic rates from diel dissolved oxygen (DO) curves (see Grace et al. 2015).
#'
#' @param data.dir 	relative or absolute path to the folder containing csv input data files to be read.
#' @param results.dir 	relative or absolute path to the output folder where results (plots and tables) will be written.
#' @param interval 	Integer. The time interval in seconds (e.g. 10 minutes = 600 seconds)
#' @param n.iter 		Integer. Number of MCMC iterations (default = 20000)
#' @param n.burnin 	Integer. Number of iterations of MCMC chains to delete
#' @param update.chains 	Logical. Should the chains automatically update once if not converged? (default = TRUE)
#' @param extra.iter Numeric. Number of extra iterations to run if chains are not converged, as multiple of n.iter (default = 1 times)
#' @param smooth.DO 	Numeric. Proportion of high-frequency fluctuations to filter with fast Fourier transform (default = 0)
#' @param smooth.PAR 	Logical. Should PAR be smoothed with a moving average? (default = FALSE)
#' @param K.init 	Numeric. Initial value of chains for K (/day). Reasonable estimate aids convergence. (default value = 2)
#' @param K.est 		Logical. Should K be estimated with uninformative priors? (default = TRUE)
#' @param K.meas.mean 	Numeric. Mean for informed normal prior distribution when K.est = FALSE
#' @param K.meas.sd 	Numeric. Standard deviation for informed normal prior distribution when K.est = FALSE
#' @param p.est	Logical. Should p be estimated? (default = FALSE)
#' @param theta.est	Logical. Should theta be estimated? (default = FALSE)
#' @param instant 		Logical. Should a table of instantaneous rates be written? (default = FALSE)
#'
#' @return A dataframe and csv file of parameter estimates (mean, SD) and checks of model fit, plots of model fit (see Vignette for details https://github.com/dgiling/BASEmetab/blob/master/vignettes/BASEmetab.pdf).
#'
#'@references Grace et al. (2015) Fast processing of diel oxygen curves: estimating stream metabolism with BASE (BAyesian Single-station Estimation). Limnology and Oceanography: Methods, 13, 103-114.
#'
#' @author Darren Giling, Ralph Mac Nally
#' @examples
#'
#' ##Link to JAGS
#' library(R2jags)
#' 
#' ##View example data set.
#' #set path to example data.
#' data.dir <- system.file("extdata", package = "BASEmetab")
#' ex.data <- read.csv(file.path(data.dir, "Yallakool_example.csv"))
#' head(ex.data)
#' tail(ex.data)
#'
#' ##Run Example.
#'
#' #set output directory to Output folder in current working directory.
#' results.dir <- file.path(getwd(), "Output")
#' if (dir.exists(results.dir)){} else {
#' dir.create(results.dir)}
#'
#' #run model.
#' results <- bayesmetab(data.dir, results.dir, interval=600)
#'
#' @export
#' @import R2jags

bayesmetab <- function(data.dir, results.dir, interval, n.iter=20000, n.burnin=n.iter*0.5, K.init = 2, 
                       smooth.DO=0, smooth.PAR=FALSE, instant=FALSE, update.chains = TRUE, extra.iter=1,
                       K.est = TRUE, K.meas.mean = 0, K.meas.sd = 4, p.est=FALSE, theta.est=FALSE) 
{

data <- read.csv('Yallakool_example.csv', stringsAsFactors = F)
  
# Set up output tables
output.table<-data.frame(File=character(), Date=character(), 
                         GPP.mean=double(), GPP.sd=double(), GPP.median=double(),
                         ER.mean=double(), ER.sd=double(), ER.median=double(), 
                         NEP.mean=double(), NEP.sd=double(), NEP.median=double(), 
                         PR.mean=double(), PR.sd=double(), PR.median=double(),
                         K.mean=double(), K.sd=double(), K.median=double(),
                         theta.mean=double(), theta.sd=double(), theta.median=double(),
                         A.mean=double(), A.sd=double(), A.median=double(), 
                         p.mean=double(), p.sd=double(), p.median=double(),
                         R2=double(), PPP=double(), rmse=double(), rmse.relative=double(), mrl.fraction=double(), ER.K.cor=double(), 
                         convergence.check=double(), A.Rhat=double(), K.Rhat=double(), theta.Rhat=double(), p.Rhat=double(), R.Rhat=double(), GPP.Rhat=double(), 
                         DIC=double(), pD=double(),
                         totDailyLight=double(), aveDailyTemp=double(), 
                         interval= double(), smooth.DO=double() , smooth.PAR=logical(), n.iter= double(), n.burnin= double(),
                         stringsAsFactors=FALSE)
instant.rates<-data.frame(File=character(), Date=character(), interval=integer(), 
                          tempC=double(), I=double(), K.instant=double(), GPP.instant=double(), ER.instant=double(),
                          stringsAsFactors=FALSE)


## checks

# Select dates
data$Date <- factor(data$Date, levels = unique(data$Date))
dates <- unique(data$Date)

## Analyse days sequentially
for (d in dates){ 
  
  data.sub <- data[data$Date == d,]
  
  # Define data vectors
  num.measurements <- nrow(data.sub)
  tempC <- data.sub$tempC
  salinity <-data.sub$salinity
  atmo.pressure <- data.sub$atmo.pressure
  DO.meas <- data.sub$DO.meas
  PAR     <- data.sub$I
  
  # Initial values
  # Set these to something sensible if the model is becoming stuck in a bad parameter space
  # These values here are expressed per timestep, not per day. Divide desired initial K (/day) by the number of timesteps in a day, as shown in default below 
  inits <- function()	{	list(K = K.init / (86400/interval) ) }
  
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
  data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","salinity","atmo.pressure", "K.init", 
                    "K.est.n", "K.meas.mean.ts", "K.meas.sd.ts", "p.est.n", "theta.est.n")  
  
  # Define monitoring variables
  params=c("A","R","K","K.day","p","theta","tau","ER","GPP","NEP","PR","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled",
           "gppts", "erpts", "kpts")
  
  ## Call jags ##
  
  # Set debug = T below to inspect each file for model convergence 
  # (inspect the main parameters for convergence using bgr diagrams, history, density and autocorrelation)
  metabfit=NULL
  metabfit <- do.call(R2jags::jags.parallel,
                      list(data=data.list, inits=inits, parameters.to.save=params, model.file = file.path(system.file(package="BASEmetab"), "BASE_metab_model_v2.3.txt"),
                           n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                           n.thin = n.thin, n.cluster= n.chains, DIC = TRUE,
                           jags.seed = 123, digits=5))
  
  # print(metabfit, digits=2) # to inspect results of last metabfit
  
  ## diagnostic summaries
  # Rhat (srf) test
  srf<- metabfit$BUGSoutput$summary[,8]
  Rhat.test <- NULL
  Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
  
  # Check for convergence and update once if requested
  if(update.chains == TRUE) {
    if(Rhat.test == "Check convergence") {
      recompile(metabfit)
      metabfit <- update(metabfit, n.iter=n.iter*extra.iter) 
      
      # Rhat (srf) test - second round in case metabfit is updated
      srf<- metabfit$BUGSoutput$summary[,8]
      Rhat.test <- NULL
      Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
    }
  }
  
  # autocorr test
  metabfit.mcmc<-coda::as.mcmc(metabfit)
  ac.lag1 <- coda::autocorr.diag(metabfit.mcmc, lags = 1)
  auto.corr.test <- NULL
  auto.corr.test <- ifelse(any(abs(ac.lag1)>0.1, na.rm=T)==TRUE,"Check ac", "ac OK")
  
  PPP <- metabfit$BUGSoutput$summary["PPfit","mean"] # posterior predictive p-value
  
  DO.mod.means <- metabfit$BUGSoutput$mean$DO.modelled
  DO.mod.sd <- metabfit$BUGSoutput$sd$DO.modelled
  
  R2 = cor(DO.mod.means,DO.meas)^2
  rmse = sqrt(sum((metabfit$BUGSoutput$mean$DO.modelled-DO.meas)^2)/length(DO.meas))
  post.mean.dev <- metabfit$BUGSoutput$mean$deviance
  pD <- metabfit$BUGSoutput$pD
  DIC <- metabfit$BUGSoutput$DIC
  
  DO.lag<-DO.meas[2:length(DO.meas)]-DO.meas[1:(length(DO.meas)-1)]
  ptpvar <- sqrt((sum((DO.lag)^2)/(length(DO.meas)-1))) # point to point variation
  rmse.relative <- rmse / ptpvar
  
  diff<-metabfit$BUGSoutput$mean$DO.modelled-DO.meas
  mrl.max<-max(rle(sign(as.vector(diff)))$lengths)
  mrl.fraction<-max(rle(sign(as.vector(diff)))$lengths)/length(DO.meas) # proportion of largest run
  
  ER.K.cor <- cor(metabfit$BUGSoutput$sims.list$ER,metabfit$BUGSoutput$sims.list$K) # plot(metabfit$sims.list$ER ~ metabfit$sims.list$K)
  
  # insert results to table and write table
  result <- data.frame(File=as.character(fname), Date=as.character(d), 
                       metabfit$BUGSoutput$mean$GPP, metabfit$BUGSoutput$sd$GPP, metabfit$BUGSoutput$median$GPP,
                       metabfit$BUGSoutput$mean$ER, metabfit$BUGSoutput$sd$ER, metabfit$BUGSoutput$median$ER,
                       metabfit$BUGSoutput$mean$NEP, metabfit$BUGSoutput$sd$NEP, metabfit$BUGSoutput$median$NEP, 
                       metabfit$BUGSoutput$mean$PR, metabfit$BUGSoutput$sd$PR, metabfit$BUGSoutput$median$PR,
                       metabfit$BUGSoutput$mean$K.day, metabfit$BUGSoutput$sd$K.day, metabfit$BUGSoutput$median$K.day,  
                       metabfit$BUGSoutput$mean$theta, metabfit$BUGSoutput$sd$theta, metabfit$BUGSoutput$median$theta, 
                       metabfit$BUGSoutput$mean$A, metabfit$BUGSoutput$sd$A, metabfit$BUGSoutput$median$A,
                       metabfit$BUGSoutput$mean$p, metabfit$BUGSoutput$sd$p, metabfit$BUGSoutput$median$p,
                       R2, PPP, rmse, rmse.relative, mrl.fraction, ER.K.cor, 
                       Rhat.test, metabfit$BUGSoutput$summary["A",8] , metabfit$BUGSoutput$summary["K",8], metabfit$BUGSoutput$summary["theta",8], 
                       metabfit$BUGSoutput$summary["p",8], metabfit$BUGSoutput$summary["R",8], metabfit$BUGSoutput$summary["GPP",8],  
                       DIC, pD,
                       totDailyLight=sum(PAR), aveDailyTemp=mean(tempC),
                       interval= interval, smooth.DO=smooth.DO , smooth.PAR=smooth.PAR, n.iter= n.iter, n.burnin= n.burnin,
                       stringsAsFactors = FALSE)
  output.table[nrow(output.table)+1,] <- result
  
  # insert results to instantaneous table and write
  if(instant == TRUE) {
    instant.result <- data.frame(File=as.character(rep(fname,seconds/interval)), Date=as.character(rep(d,seconds/interval)),interval=1:(seconds/interval),
                                 tempC=tempC, I=PAR, 
                                 K.instant=as.vector(metabfit$BUGSoutput$mean$kpts),
                                 GPP.instant=as.vector(metabfit$BUGSoutput$mean$gppts),
                                 ER.instant=as.vector(metabfit$BUGSoutput$mean$erpts),
                                 stringsAsFactors = FALSE)
    instant.rates[(nrow(instant.rates)+1):(nrow(instant.rates)+(seconds/interval)),] <- instant.result
  }
  
}
  

