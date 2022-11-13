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

# Fwoxy apa comparison to EBASE for different priors
priorcomp <- function(apasumdat, ind){
  
  met <- tibble(
    lbs = c('r2', 'rmse', 'aved'),
    lbspr = c('R^2', 'RMSE', 'Ave.\nDiff.'), 
    direc = c(-1, 1, 1)
  )
  
  toshw <- met$lbs[ind]
  leglb <- met$lbspr[ind]
  direc <- met$direc[ind]
  
  toplo <- apasumdat %>% 
    select(-out, -amean, -rmean, -bmean) %>% 
    unnest('ests') %>% 
    select(asd, rsd, bsd, var, matches(toshw)) %>%
    filter(!var %in% 'b') %>% 
    pivot_wider(names_from = 'var', values_from = !!toshw) %>% 
    mutate(
      ind = 1:nrow(.)
    )
  
  toplo1 <- toplo %>% 
    select(ind, asd, rsd, bsd) %>% 
    mutate(
      def = case_when(
        asd == 0.1 & rsd == 5 & bsd == 0.01 ~ '*', 
        T ~ ''
      ),
      asd = factor(asd, labels = c('L', 'M', 'H')), 
      rsd = factor(rsd, labels = c('L', 'M', 'H')), 
      bsd = factor(bsd, labels = c('L', 'M', 'H')), 
    ) %>% 
    pivot_longer(-c('ind', 'def'), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('asd', 'rsd', 'bsd'), 
                   labels = c('a', 'r', 'b'))
    ) 
  
  toplo2 <- toplo %>% 
    select(-asd, -rsd, -bsd) %>%  
    pivot_longer(-ind, names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('DO_mod', 'Pg_vol', 'Rt_vol', 'D', 'a'), 
                   labels = c('DO [mod]', 'Pg [vol]', 'Rt [vol]', 'D', 'a')
      )
    )
  
  p1 <- ggplot(toplo1, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    theme(
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size = 12),
      axis.ticks = element_blank(), 
      legend.position = 'left', 
      legend.title = element_blank()
    ) + 
    scale_fill_brewer(palette = 'Greys') + 
    scale_x_discrete(position = 'top', expand = c(0, 0)) + 
    scale_y_reverse(expand = c(0, 0), breaks = toplo1$ind, labels = toplo1$def) + 
    labs(
      y = NULL, 
      x = 'Variance of prior', 
      caption = '* EBASE default'
    )
  
  p2 <- ggplot(toplo2, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    theme(
      axis.text.x = element_text(face = 'italic', size = 12), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      legend.position = 'right'
    ) + 
    scale_fill_distiller(palette = 'YlOrRd', direction = direc) + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo2$var))) + 
    scale_y_reverse(expand = c(0, 0)) + 
    labs(
      y = NULL, 
      fill = parse(text = leglb),
      x = 'Parameter from EBASE vs Fwoxy'
    )
  
  out <- p1 + p2 + plot_layout(ncol = 2, widths = c(0.5, 1))
  
  return(out)
  
}