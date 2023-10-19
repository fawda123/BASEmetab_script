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

# Fwoxy apa comparison to EBASE for different priors sd only
priorcomp <- function(dat, ind){
  
  met <- tibble(
    lbs = c('r2', 'rmse'),
    lbspr = c('R^2', 'RMSE'), 
    direc = c(-1, 1)
  )

  toshw <- met$lbs[ind]
  leglb <- met$lbspr[ind]
  direc <- met$direc[ind]

  toplo <- dat %>% 
    select(-amean, -rmean, -bmean) %>% 
    unnest('ests') %>% 
    select(ndays, asd, rsd, bsd, var, matches(toshw)) %>%
    filter(!var %in% 'b') %>% 
    pivot_wider(names_from = 'var', values_from = !!toshw) %>% 
    mutate(
      ind = sort(rep(1: (nrow(.) / 2), times = 2)), 
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      )
    )
  
  toplo1 <- toplo %>% 
    select(ind, asd, rsd, bsd) %>% 
    unique() %>% 
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
    pivot_longer(-c(ind, ndays), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('DO_mod', 'P', 'R', 'D', 'a'), 
                   labels = c('DO [mod]', 'P', 'R', 'D', 'a')
      )
    )
  
  flmx <- 100
  if(toshw == 'rmse')
    flmx <- max(toplo2$val)
  
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
      legend.position = 'right', 
      strip.background = element_blank(), 
      strip.text = element_text(hjust = 0, size = 12, face = 'bold')
    ) + 
    facet_wrap(~ndays, ncol = 2) + 
    scale_fill_distiller(palette = 'YlOrRd', direction = direc, limits = c(0, flmx)) + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo2$var))) + 
    scale_y_reverse(expand = c(0, 0)) + 
    labs(
      y = NULL, 
      fill = parse(text = leglb),
      x = 'Parameter from EBASE vs Fwoxy'
    )

  out <- p1 + p2 + plot_layout(ncol = 2, widths = c(0.3, 1))
  
  return(out)
  
}

# Fwoxy apa comparison to EBASE for different priors, b mean and sd changes only
priorcompmean <- function(dat, ind){
  
  met <- tibble(
    lbs = c('r2', 'rmse', 'aved'),
    lbspr = c('R^2', 'RMSE', 'Ave.\nDiff.'), 
    direc = c(-1, 1, 1)
  )
  
  toshw <- met$lbs[ind]
  leglb <- met$lbspr[ind]
  direc <- met$direc[ind]

  toplo <- dat %>% 
    select(-amean, -asd, -rmean, -rsd) %>% 
    unnest('ests') %>% 
    select(ndays, bmean, bsd, var, matches(toshw)) %>%
    filter(!var %in% 'b') %>%
    pivot_wider(names_from = 'var', values_from = !!toshw) %>% 
    mutate(
      ind = sort(rep(1: (nrow(.) / 2), times = 2)), 
      ndays = case_when(
        ndays == 1 ~ paste(ndays, 'day'), 
        T ~ paste(ndays, 'days')
      )
    )
  
  toplo1 <- toplo %>% 
    select(ind, bmean, bsd) %>% 
    unique() %>% 
    mutate(
      def = case_when(
        bmean == 0.251 & bsd == 0.01 ~ '*', 
        T ~ ''
      ),
      bmean = factor(bmean, labels = c('L', 'M', 'H')),
      bsd = factor(bsd, labels = c('L', 'M', 'H'))
    ) %>% 
    pivot_longer(-c('ind', 'def'), names_to = 'var', values_to = 'val') %>% 
    mutate(
      var = factor(var, 
                   levels = c('bmean', 'bsd'), 
                   labels = c('mean', 'sd'))
    ) 
  
  toplo2 <- toplo %>% 
    select(-bmean, -bsd) %>%
    pivot_longer(-c(ind, ndays), names_to = 'var', values_to = 'val') %>% 
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
      x = 'b prior',
      caption = '* EBASE default'
    )
  
  p2 <- ggplot(toplo2, aes(y = ind, x = var, fill = val)) + 
    geom_tile(color = 'black') + 
    theme(
      axis.text.x = element_text(face = 'italic', size = 12), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      legend.position = 'right', 
      strip.background = element_blank(), 
      strip.text = element_text(hjust = 0, size = 12, face = 'bold')
    ) + 
    facet_wrap(~ndays, ncol = 2) + 
    scale_fill_distiller(palette = 'YlOrRd', direction = direc, limits = c(0, 100)) + 
    scale_x_discrete(position = 'top', expand = c(0, 0), labels = parse(text = levels(toplo2$var))) + 
    scale_y_reverse(expand = c(0, 0)) +
    labs(
      y = NULL, 
      fill = parse(text = leglb),
      x = 'Parameter from EBASE vs Fwoxy'
    )
  
  out <- p1 + p2 + plot_layout(ncol = 2, widths = c(0.3, 1))
  
  return(out)
  
}

# comparison of fwoxy and ebase results for selected prior at time step of ndays
optex <- function(apagrd, fwdatcmp, asdin, rsdin, bsdin, ndaysin){
  
  res <- apagrd %>% 
    filter(
      asd == asdin & rsd == rsdin & bsd == bsdin & ndays == ndaysin
    ) %>% 
    pull(out) %>% 
    .[[1]] %>% 
    .[[1]]
  cmp <- inner_join(fwdatcmp, res, by = c('Date', 'DateTimeStamp')) %>%
    select(-converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>%
    rename(
      DO_mod.x = DO_obs.x,
      DO_mod.y = DO_mod
    ) %>%
    pivot_longer(!all_of(c('DateTimeStamp', 'Date', 'grp')), names_to = 'var', values_to = 'val') %>%
    separate(var, c('var', 'mod'), sep = '\\.') %>%
    mutate(
      mod = case_when(
        mod == 'x' ~ 'Fwoxy',
        mod == 'y' ~ 'EBASE'
      )
    ) %>%
    pivot_wider(names_from = 'mod', values_from = 'val')

  toplo1 <- cmp %>% 
    filter(var %in% c('P', 'R', 'D')) %>% 
    group_by(grp, var) %>% 
    summarise(
      Fwoxy = mean(Fwoxy, na.rm = T), 
      EBASE = mean(EBASE, na.rm = T),
      Date  = min(Date),
      .groups = 'drop'
    ) %>% 
    pivot_longer(-c(Date, grp, var), names_to = 'model', values_to = 'est') %>% 
    mutate(
      var = factor(var, 
                   levels = c('P', 'R', 'D')
      )
    ) %>% 
    select(-grp)
  
  ylab <- expression(paste(O [2], ' (mmol ', m^-2, ' ', d^-1, ')'))
  
  p1 <- ggplot(toplo1, aes(x = Date, y = est, group = model, color = model)) + 
    geom_line() +
    geom_point() + 
    facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y', labeller = label_parsed) + 
    theme_minimal() + 
    theme(
      strip.placement = 'outside', 
      strip.background = element_blank(), 
      legend.position = 'top', 
      legend.title = element_blank(),
      strip.text = element_text(size = rel(1))
    ) + 
    labs(
      x = NULL, 
      y = ylab
    )
  
  labs <- c('DO[mod]~(mmol~O[2]~m^{3}~d^{-1})',
            'a~(mmol~m^{-3}~d^{-1})/(W~m^{-2})', 
            'b~(cm~hr^{-1})/(m^{2}~s^{-2})'
  )
  
  toplo2 <- cmp %>% 
    filter(var %in% c('DO_mod', 'a', 'b')) %>% 
    group_by(grp, var) %>%
    summarise(
      Fwoxy = mean(Fwoxy, na.rm = T),
      EBASE = mean(EBASE, na.rm = T),
      Date = min(Date),
      .groups = 'drop'
    ) %>%
    pivot_longer(-c(Date, grp, var), names_to = 'model', values_to = 'est') %>% 
    mutate(
      var = factor(var, 
                   levels = c('DO_mod', 'a', 'b'), 
                   labels = labs
      )
    ) %>% 
    select(-grp)
  
  p2 <- ggplot(toplo2, aes(x = Date, y = est, group = model, color = model)) + 
    geom_line() +
    geom_point() +
    facet_wrap(~var, ncol = 1, strip.position = 'left', scales = 'free_y', labeller = label_parsed) + 
    theme_minimal() + 
    theme(
      strip.placement = 'outside', 
      strip.background = element_blank(), 
      legend.position = 'top', 
      legend.title = element_blank(), 
      strip.text = element_text(size = rel(1))
    ) + 
    labs(
      x = NULL, 
      y = NULL
    )
  
  p1 + p2 + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = 'top')
  
}

# summary function for r2, rmse, and ave diff, used in fwoxyapacmp.R
sumfun <- function(x){
  
  lmmod <- lm(EBASE ~ Fwoxy, x) 
  
  # r2
  r2 <- lmmod %>% 
    summary() %>% 
    .$r.squared
  r2 <- 100 * r2
  
  # rmse
  rmse <- sqrt(mean((x$EBASE - x$Fwoxy)^2, na.rm = T))

  # mape
  mape <- abs((x$Fwoxy - x$EBASE) / x$Fwoxy)
  mape <- 100 * mean(mape, na.rm = TRUE)

  # mae
  mae <- mean(abs(x$Fwoxy - x$EBASE), na.rm = TRUE)

  # nash-sutcliffe efficiency
  # rss <- sum((x$Fwoxy - x$EBASE)^2)# sum(resid(lmmod)^2) # same as deviance(lmmod)
  # tss <- sum((x$Fwoxy - mean(x$Fwoxy, na.rm = T))^2)
  # cfd <- 1 - rss / tss
  nse <- NSE(x$EBASE, x$Fwoxy)
  
  out <- data.frame(r2 = r2, rmse = rmse, mape = mape, mae = mae, nse = nse)
  
  return(out)
  
}

# plot comparison of b for changing mean and sd
optexmean <- function(dat){
  
  toplo <- dat %>% 
    unnest('out') %>% 
    mutate(
      out = purrr::map(out, function(x){
        
        x %>% 
          group_by(grp) %>% 
          summarise(
            Date = min(Date),
            EBASE = mean(b), 
            .groups = 'drop'
          ) %>% 
          mutate(
            Fwoxy = 0.251
          )
        
      })
    ) %>% 
    unnest('out') %>% 
    select(ndays, bmean, bsd, Date, EBASE, Fwoxy) %>% 
    mutate(
      bmean = factor(bmean, labels = paste('mean:', c('L', 'M', 'H'))),
      bsd = factor(bsd, labels = paste('sd:', c('L', 'M', 'H')))
    ) %>% 
    pivot_longer(c(EBASE, Fwoxy), names_to = 'model', values_to = 'est') %>% 
    na.omit()
  
  toplo1 <- toplo %>% 
    filter(ndays == 1)
  
  toplo2 <- toplo %>% 
    filter(ndays == 7)
  
  thm <- theme_bw() + 
    theme(
      strip.placement = 'outside', 
      strip.background = element_blank(), 
      legend.position = 'top', 
      legend.title = element_blank(), 
      strip.text = element_text(size = rel(1)), 
      axis.text.x = element_text(size = 8)
    ) 
  
  p1 <- ggplot(toplo1, aes(x = Date, y = est, color = model)) +
    geom_line() +
    geom_point() +
    facet_grid(bmean ~ bsd) + 
    thm + 
    labs(
      x = NULL, 
      y = expression(paste('b (cm ', hr^-1, ') / ( ', m^2 ~ s^-2, ')')), 
      title = 'ndays = 1'
    )
  
  p2 <- ggplot(toplo2, aes(x = Date, y = est, color = model)) +
    geom_line() +
    geom_point() +
    facet_grid(bmean ~ bsd) + 
    thm +
    labs(
      x = NULL, 
      y = expression(paste('b (cm ', hr^-1, ') / ( ', m^2 ~ s^-2, ')')), 
      title = 'ndays = 7'
    )
  
  p <- p1 + p2 + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = 'bottom')
  
  return(p)
  
}

optexmeansum <- function(dat){
  
  toplo <- dat %>% 
    unnest('out') %>% 
    mutate(
      out = purrr::map(out, function(x){
        
        x %>% 
          group_by(grp) %>% 
          summarise(
            Date = min(Date),
            EBASE = mean(b), 
            .groups = 'drop'
          ) %>% 
          mutate(
            Fwoxy = 0.251
          )
        
      })
    ) %>% 
    unnest('out') %>% 
    select(ndays, bmean, bsd, Date, EBASE, Fwoxy) %>% 
    mutate(
      bmean = factor(bmean, labels = paste('mean:', c('L', 'M', 'H'))),
      bsd = factor(bsd, labels = paste('sd:', c('L', 'M', 'H')))
    ) %>% 
    pivot_longer(c(EBASE, Fwoxy), names_to = 'model', values_to = 'est') %>% 
    na.omit() %>% 
    group_by(ndays, bmean, bsd, model) %>% 
    summarise(
      avev = mean(est), 
      hiv = ifelse(model == 'EBASE', t.test(est)$conf.int[2], est),
      lov = ifelse(model == 'EBASE', t.test(est)$conf.int[1], est), 
      .groups = 'drop'
    ) %>% 
    unique()
  
  toplo1 <- toplo %>% 
    filter(ndays == 1)
  
  toplo2 <- toplo %>% 
    filter(ndays == 7)
  
  thm <- theme_bw() +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      strip.background = element_blank(), 
      axis.title.x = element_blank(), 
      legend.position  = 'bottom', 
      panel.grid.major.x = element_blank(), 
      panel.grid.minor.x = element_blank()
    )
  
  ddg <- position_dodge(width = 0.001)
  
  p1 <- ggplot(toplo1, aes(x = 1, y = avev, color = model)) + 
    geom_point(position = ddg) + 
    geom_errorbar(aes(ymin = lov, ymax = hiv), width = 0, position = ddg) +
    facet_wrap(~bsd + bmean, ncol = 9, strip.position = 'bottom') + 
    thm +
    labs(
      color = NULL,
      y = expression(paste('b (cm ', hr^-1, ') / ( ', m^2 ~ s^-2, ')')), 
      title = 'ndays = 1'
    )
  
  p2 <- ggplot(toplo2, aes(x = 1, y = avev, color = model)) + 
    geom_point(position = ddg) + 
    geom_errorbar(aes(ymin = lov, ymax = hiv), width = 0, position = ddg) +
    facet_wrap(~bsd + bmean, ncol = 9, strip.position = 'bottom') + 
    thm +
    labs(
      color = NULL,
      y = expression(paste('b (cm ', hr^-1, ') / ( ', m^2 ~ s^-2, ')')), 
      title = 'ndays = 7'
    )
  
  p <- p1 + p2 + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = 'bottom')
  
  return(p)
  
}


NSE <-function(sim, obs, ...) UseMethod("NSE")

NSE.default <- function (sim, obs, na.rm=TRUE, fun=NULL, ..., 
                         epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                         epsilon.value=NA){ 
  
  if ( is.na(match(class(sim), c("integer", "numeric", "ts", "zoo", "xts"))) |
       is.na(match(class(obs), c("integer", "numeric", "ts", "zoo", "xts")))
  ) stop("Invalid argument type: 'sim' & 'obs' have to be of class: c('integer', 'numeric', 'ts', 'zoo', 'xts')")      
  
  epsilon.type <- match.arg(epsilon.type)  
  
  # index of those elements that are present both in 'sim' and 'obs' (NON- NA values)
  vi <- valindex(sim, obs)
  
  if (length(vi) > 0) {	 
    # Filtering 'obs' and 'sim', selecting only those pairs of elements 
    # that are present both in 'x' and 'y' (NON- NA values)
    obs <- obs[vi]
    sim <- sim[vi]
    
    if (!is.null(fun)) {
      fun1 <- match.fun(fun)
      new  <- preproc(sim=sim, obs=obs, fun=fun1, ..., 
                      epsilon.type=epsilon.type, epsilon.value=epsilon.value)
      sim  <- new[["sim"]]
      obs  <- new[["obs"]]
    } # IF end     
    
    denominator <- sum( (obs - mean(obs))^2 )
    
    if (denominator != 0) {      
      NS <- 1 - ( sum( (obs - sim)^2 ) / denominator )     
    } else {
      NS <- NA
      warning("'sum((obs - mean(obs))^2)=0' => it is not possible to compute 'NSE'")  
    } 
  } else {
    NS <- NA
    warning("There are no pairs of 'sim' and 'obs' without missing values !")
  } # ELSE end
  
  return(NS)
  
} # 'NSE' end


NSE.matrix <- function(sim, obs, na.rm=TRUE, fun=NULL, ..., 
                       epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                       epsilon.value=NA){ 
  
  # Checking that 'sim' and 'obs' have the same dimensions
  if ( all.equal(dim(sim), dim(obs)) != TRUE )
    stop( paste("Invalid argument: dim(sim) != dim(obs) ( [", 
                paste(dim(sim), collapse=" "), "] != [", 
                paste(dim(obs), collapse=" "), "] )", sep="") )
  
  NS <- rep(NA, ncol(obs))       
  
  NS <- sapply(1:ncol(obs), function(i,x,y) { 
    NS[i] <- NSE.default( x[,i], y[,i], na.rm=na.rm, fun=fun, ..., 
                          epsilon.type=epsilon.type, epsilon.value=epsilon.value)
  }, x=sim, y=obs )    
  
  names(NS) <- colnames(obs)
  
  return(NS)
  
} # 'NSE.matrix' end


NSE.data.frame <- function(sim, obs, na.rm=TRUE, fun=NULL, ..., 
                           epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                           epsilon.value=NA){ 
  
  sim <- as.matrix(sim)
  obs <- as.matrix(obs)
  
  NSE.matrix(sim, obs, na.rm=na.rm, fun=fun, ..., 
             epsilon.type=epsilon.type, epsilon.value=epsilon.value)
  
} # 'NSE.data.frame' end


NSeff <-function(sim, obs, ...) UseMethod("NSE")


################################################################################
# Author: Mauricio Zambrano-Bigiarini                                          #
################################################################################
# Started: 22-Mar-2013                                                         #
# Updates: 12-Jul-2022 ; 13-Jul-2022                                           #
################################################################################
NSE.zoo <- function(sim, obs, na.rm=TRUE, fun=NULL, ..., 
                    epsilon.type=c("none", "Pushpalatha2012", "otherFactor", "otherValue"), 
                    epsilon.value=NA){ 
  
  sim <- zoo::coredata(sim)
  if (is.zoo(obs)) obs <- zoo::coredata(obs)
  
  if (is.matrix(sim) | is.data.frame(sim)) {
    NSE.matrix(sim, obs, na.rm=na.rm, fun=fun, epsilon.type=epsilon.type, 
               epsilon.value=epsilon.value, ...)
  } else NextMethod(sim, obs, na.rm=na.rm, fun=fun, epsilon.type=epsilon.type, 
                    epsilon.value=epsilon.value, ...)
  
} # 'NSE.zoo' end

valindex <- function(sim, obs, ...) UseMethod("valindex")

valindex.default <- function(sim, obs, ...) {  
  
  if ( length(obs) != length(sim) ) {
    stop( "Invalid argument: 'length(sim) != length(obs)' !! (", length(sim), "!=", length(obs), ") !!" )
  } else { 
    index <- which(!is.na(sim) & !is.na(obs))
    if (length(index)==0) warning("'sim' and 'obs' are empty or they do not have any common pair of elements with data !!")
    return( index  )
  } # ELSE end
  
} # 'valindex' END


################################################################################
# Author: Mauricio Zambrano-Bigiarini                                          #
################################################################################
# Started: 25-Jul-2011                                                         #
# Updates: 08-May-2012                                                         #
################################################################################
valindex.matrix <- function(sim, obs, ...) { 
  
  # Checking that 'sim' and 'obs' have the same dimensions
  if ( all.equal(dim(sim), dim(obs)) != TRUE ) {
    stop( paste("Invalid argument: dim(sim) != dim(obs) ( [", 
                paste(dim(sim), collapse=" "), "] != [", 
                paste(dim(obs), collapse=" "), "] )", sep="") )
  } else  
    return ( !is.na( sim) & !is.na(obs) )
  
} # 'valindex.matrix' END