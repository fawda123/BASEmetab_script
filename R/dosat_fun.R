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