# calculate DO sat as in Grace et al. 
#
# tempC
# salinity 
# atmo.pressure
dosat_fun <- function(tempC, salinity, atmo.pressure){
  
  kelvin <- 273.15 + tempC
  
  # correction for salinity
  
  S1 <- 157570.1 / kelvin
  S2 <- -6.6423080E7 / (kelvin * kelvin)
  S3 <- 1.2438E10 / (kelvin ^ 3)
  S4 <-  -8.621949E11 / (kelvin ^ 4)
  sal.factor <- -1.0 * salinity * (0.017674 - 10.754 / kelvin + 2140.7 / (kelvin * kelvin))
  
  DOsalinity.corr <-exp(-139.34411+S1+S2+S3+S4+sal.factor)
  
  # correction for atmospheric pressure
  alpha <- 0.000975 - 0.00001426 * kelvin + 0.00000006436 * (kelvin ^ 2)		
  beta <- exp(11.8571 - 3840.7 / kelvin - 216961 / (kelvin ^ 2))
  gamma <- ((1 - beta / atmo.pressure) / (1 - beta)) * ((1 - alpha * atmo.pressure) / (1 - alpha))
  
  DO.sat <- DOsalinity.corr * atmo.pressure * gamma		
  
  return(DO.sat)
  
}

# update metabolism jags fit
#
# metabfit initial jags metabolism output
# update.chains logical to update, only if TRUE
# n.iter number of iterations 
metab_update <- function(metabfit, update.chains, n.iter){

  ## diagnostic summaries
  # Rhat (srf) test
  srf <- metabfit$BUGSoutput$summary[,8]
  Rhat.test <- NULL
  Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
  
  # Check for convergence and update once if requested
  if(update.chains == TRUE) {
    if(Rhat.test == "Check convergence") {
      recompile(metabfit)
      metabfit <- update(metabfit, n.iter=n.iter*1) 
      
      # Rhat (srf) test - second round in case metabfit is updated
      srf <- metabfit$BUGSoutput$summary[,8]
      Rhat.test <- NULL
      Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
    }
  }
  
  return(metabfit)
  
}

# function to make seconds in opmetab match fwoxy
getsec <- function(x){
  
  sec <- seconds(x)
  sec <- sec - (7.5 * 60)
  sec <- sec - min(sec)
  sec <- 1 + as.numeric(sec)
  
  return(sec)
  
}
