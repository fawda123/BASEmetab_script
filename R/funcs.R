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

# evalcor from wtregdo w/ sun angle floored at zero
evalcorflr <- function(dat_in, tz, lat, long, depth_val = 'Tide', daywin = 6, method = 'pearson', plot = TRUE, lims = c(-0.5, 0.5), progress = FALSE, harm = TRUE, chk_tide = FALSE, constituents = c('M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'MF', 'MM', 'SSA', 'M4', 'M6', 'S4', 'MS4')
){
  
  # sanity check
  chktz <- attr(dat_in$DateTimeStamp, 'tzone')
  if(tz != chktz)
    stop('dat_in timezone differs from tz argument')
  
  # check for duplicated rows
  chk <- duplicated(dat_in)
  if(any(chk))
    stop('Duplicated observations found, check rows: ', paste(which(chk), collapse = ', '))
  
  names(dat_in)[names(dat_in) %in% depth_val] <- 'Tide'
  
  # get decimal time
  tocor <- met_day_fun(dat_in, tz, lat, long)
  tocor$dec_time <- 365 * (lubridate::decimal_date(tocor$DateTimeStamp) - max(lubridate::year(tocor$DateTimeStamp)))
  
  # sun angle
  locs <- c(long, lat)
  utc_time <- as.POSIXlt(tocor$DateTimeStamp, tz = 'UTC')
  sun_angle <- oce::sunAngle(utc_time, locs[1], locs[2])$altitude
  sun_angle <- pmax(0, sun_angle)

  # get tidal predictions from harmonic regression
  if(harm){
    
    tide_mod <- oce::tidem(
      x = tocor[,'Tide'],
      t = tocor[,'DateTimeStamp'],
      constituents = constituents
    )
    tide_pred <- predict(tide_mod, newdata = tocor[,'DateTimeStamp'])
    
  }
  
  # return tidal model if TRUE
  if(chk_tide){
    
    out <- list(
      tide_mod = tide_mod,
      tide_pred = data.frame(DateTimeStamp = tocor$DateTimeStamp, tide_obs = tocor$Tide, tide_pred = tide_pred)
    )
    return(out)
    
  }
  
  # replace tidal observations with tidal predictions
  if(harm)
    tocor$Tide <- tide_pred
  
  # polar coords for tidal height, in degrees
  tocor$dTide <- with(tocor, c(diff(Tide)[1], diff(Tide)))
  tide_angle <- with(tocor, atan(dTide/c(diff(dec_time)[1], diff(dec_time))) * 180/pi)
  
  # for weights
  tocor$hour <- as.numeric(strftime(tocor$DateTimeStamp, '%H', tz = tz))
  tocor$hour <- tocor$hour + as.numeric(strftime(tocor$DateTimeStamp, '%M', tz = tz))/60
  
  #for counter
  strt <- Sys.time()
  
  # moving window correlations given daywin
  cor_out <- foreach(row = 1:nrow(tocor)) %dopar% {
    
    # progress
    if(progress){
      sink('log.txt')
      cat('Log entry time', as.character(Sys.time()), '\n')
      cat(row, ' of ', nrow(tocor), '\n')
      print(Sys.time() - strt)
      sink()
    }
    
    ref_in <- tocor[row, ]
    
    # return NA if no tide value at obs
    if(is.na(ref_in$Tide))
      return(NA)
    
    wts <- try({WtRegDO::wtfun(ref_in, tocor, wins = list(daywin, 1e6, 1e6))}, silent = TRUE)
    
    # return NA if error on wts, usually from insufficient data
    if(inherits(wts, 'try-error'))
      return(NA)
    
    gr_zero <- which(wts > 0)
    
    sun_in <- (sun_angle)[gr_zero]
    tide_in <- (tide_angle)[gr_zero]
    
    corout <- cor(sun_in, tide_in, method = method, use = 'complete.obs')
    
    return(corout)
    
  }
  
  cor_out <- unlist(cor_out)
  
  # create plot
  toplo <- data.frame(tocor, Correlation = cor_out)
  
  p <- ggplot(toplo, aes_string(x = 'DateTimeStamp', y = 'Correlation')) +
    geom_line() +
    scale_y_continuous(limits = lims) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    theme_bw()
  
  if(plot) return(p)
  
  return(cor_out)
  
}

# evalcor from wtregdo w/ output for one time step
evalcorstp <- function(dat_in, row, tz, lat, long, corv = TRUE, depth_val = 'Tide', daywin = 6, method = 'pearson', harm = TRUE, chk_tide = FALSE, constituents = c('M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'MF', 'MM', 'SSA', 'M4', 'M6', 'S4', 'MS4')
){
  
  # sanity check
  chktz <- attr(dat_in$DateTimeStamp, 'tzone')
  if(tz != chktz)
    stop('dat_in timezone differs from tz argument')
  
  # check for duplicated rows
  chk <- duplicated(dat_in)
  if(any(chk))
    stop('Duplicated observations found, check rows: ', paste(which(chk), collapse = ', '))
  
  names(dat_in)[names(dat_in) %in% depth_val] <- 'Tide'
  
  # get decimal time
  tocor <- met_day_fun(dat_in, tz, lat, long)
  tocor$dec_time <- 365 * (lubridate::decimal_date(tocor$DateTimeStamp) - max(lubridate::year(tocor$DateTimeStamp)))
  
  # sun angle
  locs <- c(long, lat)
  utc_time <- as.POSIXlt(tocor$DateTimeStamp, tz = 'UTC')
  sun_angle <- oce::sunAngle(utc_time, locs[1], locs[2])$altitude
  tocor$sun_angle <- sun_angle
  
  # get tidal predictions from harmonic regression
  if(harm){
    
    tide_mod <- oce::tidem(
      x = tocor[,'Tide'],
      t = tocor[,'DateTimeStamp'],
      constituents = constituents
    )
    tide_pred <- predict(tide_mod, newdata = tocor[,'DateTimeStamp'])
    
  }
  
  # return tidal model if TRUE
  if(chk_tide){
    
    out <- list(
      tide_mod = tide_mod,
      tide_pred = data.frame(DateTimeStamp = tocor$DateTimeStamp, tide_obs = tocor$Tide, tide_pred = tide_pred)
    )
    return(out)
    
  }
  
  # replace tidal observations with tidal predictions
  if(harm)
    tocor$Tide <- tide_pred
  
  # polar coords for tidal height, in degrees
  tocor$dTide <- with(tocor, c(diff(Tide)[1], diff(Tide)))
  tide_angle <- with(tocor, atan(dTide/c(diff(dec_time)[1], diff(dec_time))) * 180/pi)
  tocor$tide_angle <- tide_angle
  
  # for weights
  tocor$hour <- as.numeric(strftime(tocor$DateTimeStamp, '%H', tz = tz))
  tocor$hour <- tocor$hour + as.numeric(strftime(tocor$DateTimeStamp, '%M', tz = tz))/60
  
  ref_in <- tocor[row, ]
  
  # return NA if no tide value at obs
  if(is.na(ref_in$Tide))
    return(NA)
  
  wts <- try({WtRegDO::wtfun(ref_in, tocor, wins = list(daywin, 1e6, 1e6))}, silent = TRUE)
    
  # return NA if error on wts, usually from insufficient data
  if(inherits(wts, 'try-error'))
    return(NA)
    
  gr_zero <- which(wts > 0)

  if(!corv){
    out <- tocor[gr_zero, ]
    return(out)
  }   
  
  sun_in <- (sun_angle)[gr_zero]
  tide_in <- (tide_angle)[gr_zero]
  
  corout <- cor(sun_in, tide_in, method = method, use = 'complete.obs')
  
  return(corout)
  
}