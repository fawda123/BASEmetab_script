library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(here)
library(EBASE)
library(doParallel)
library(ggplot2)
library(patchwork)

fwdat <- read_csv(here("data/apafwoxy.csv")) 

# fwoxy for comparison
# convert areal to volumetric, extract b at right time step
fwdatcmp <- fwdat %>% 
  mutate(
    DateTimeStamp = dmy_hms(datet, tz = 'America/Jamaica'),
    Date = as.Date(DateTimeStamp, tz = 'America/Jamaica'),
    DO_obs = `oxy,mmol/m3`, 
    a = `aparam,(mmolO2/m2/d)/(W/m2)` / `ht,m`,
    Rt_vol = `er,mmol/m2/d` / `ht,m`,
    Pg_vol = `gpp,mmol/m2/d` / `ht,m`,
    D = -1 * `gasex,mmol/m2/d` / `ht,m`,
    b = 100 * 3600 * `kw,m/s` / `wspd2,m2/s2` / (`sc,dimensionless` / 600) ^ -0.5 # (m/s)/(m2/s2) to (cm/hr) / (m2/s2)
  ) %>% 
  select(Date, DateTimeStamp, DO_obs, a, b, Pg_vol, Rt_vol, D)

# fwoxy for input to ebase
fwdatinp <- fwdat %>% 
  mutate(
    datet = dmy_hms(datet, tz = 'America/Jamaica')
  ) %>% 
  select(
    DateTimeStamp = datet,
    DO_obs = `oxy,mmol/m3`, 
    Temp = `temp,degC`, 
    Sal = `salt,ppt`, 
    PAR = `par,W/m2`, 
    WSpd = `wspd2,m2/s2`, 
    H = `ht,m`
  ) %>% 
  mutate(
    DO_obs = DO_obs / 1000 * 32, # to mg/L
    WSpd = sqrt(WSpd)
  )

# # simple comparison ---------------------------------------------------------------------------
#
# # subset four days in June
# dat <- fwdatinp %>%
#   filter(month(fwdatinp$DateTimeStamp) == 6 & day(fwdatinp$DateTimeStamp) %in% 1:4)
# 
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# res <- ebase(dat, interval = 900, H = dat$H, progress = TRUE, n.chains = 4, interp = T)
# 
# stopCluster(cl)
# 
# fwdatinpcmp <- fwdatcmp %>% 
#   filter(Date <= max(res$Date) & Date >= min(res$Date))
# 
# cmp <- inner_join(fwdatcmp, res, by = c('Date', 'DateTimeStamp')) %>% 
#   select(-converge, -dDO, -DO_obs.y, -rsq, -matches('lo$|hi$')) %>% 
#   rename(
#     DO_mod.x = DO_obs.x, 
#     DO_mod.y = DO_mod
#   ) %>% 
#   pivot_longer(!all_of(c('DateTimeStamp', 'Date', 'grp')), names_to = 'var', values_to = 'val') %>% 
#   separate(var, c('var', 'mod'), sep = '\\.') %>% 
#   mutate(
#     mod = case_when(
#       mod == 'x' ~ 'Fwoxy', 
#       mod == 'y' ~ 'EBASE'
#     )
#   ) %>% 
#   pivot_wider(names_from = 'mod', values_from = 'val')
# 
# ggplot(cmp, aes(x = Fwoxy, y = EBASE, color = var)) + 
#   geom_point() +
#   facet_wrap(~var, ncol = 2, scales = 'free') +
#   theme_bw() + 
#   geom_abline(intercept = 0, slope = 1) +
#   theme(
#     legend.title = element_blank(),
#     strip.background = element_blank()
#   )

# gridded comparison --------------------------------------------------------------------------

# this takes about ten hours to run
grd <- crossing(
  amean = 0.2, #c(0, 0.2, 2),
  asd = c(0.01, 0.1, 1),
  rmean = 20, #c(0, 20, 200), 
  rsd = c(0.5, 5, 50), 
  bmean = 0.251, #c(0.0251, 0.251, 2.51), 
  bsd = c(0.001, 0.01, 0.1),
  ndays = c(1, 7),
  out = NA
)

str <- Sys.time()

# takes about 20 hours
for(i in 1:nrow(grd)){
  
  # counter
  cat(i, 'of', nrow(grd), '\n')
  print(Sys.time() - str)
  
  # get inputs
  selrow <- grd[i, ]
  aprior <- c(selrow$amean, selrow$asd)
  rprior <- c(selrow$rmean, selrow$rsd)
  bprior <- c(selrow$bmean, selrow$bsd)
  ndays <- c(selrow$ndays)
  
  # run model for inputs
  cl <- makeCluster(6)
  registerDoParallel(cl)
  
  # use interp for missing values
  res <- ebase(fwdatinp, interval = 900, H = fwdatinp$H, progress = TRUE, n.chains = 4, 
               aprior = aprior, rprior = rprior, bprior = bprior, ndays = ndays)
  
  stopCluster(cl)
  
  # append output to grd
  grd$out[[i]] <- list(res)
  
}

apagrd <- grd
save(apagrd, file = 'data/apagrd.RData', compress = 'xz')

# evaluate fit --------------------------------------------------------------------------------

load(file = here('data/apagrd.RData'))

# summary function for r2, rmse, and ave diff
sumfun <- function(x){
  
  # r2
  r2 <- lm(EBASE ~ Fwoxy, x) %>% 
    summary() %>% 
    .$r.squared
  r2 <- 100 * r2
  
  # rmse
  rmse <- sqrt(mean((x$EBASE - x$Fwoxy)^2, na.rm = T))
  
  # aved
  ts1 <- sum(x$EBASE, na.rm = TRUE)
  ts2 <- sum(x$Fwoxy, na.rm = TRUE)
  
  aved <- 100 * (ts1 - ts2)/((ts1 + ts2)/2)
  
  out <- data.frame(r2 = r2, rmse = rmse, aved = aved)
  
  return(out)
  
}

apasumdat <- apagrd %>% 
  mutate(
    ind = 1:nrow(.),
    ests = purrr::pmap(list(ind, out), function(ind, out){
      
      cat(ind, '\t')
      
      cmp <- inner_join(fwdatcmp, out[[1]], by = c('Date', 'DateTimeStamp')) %>%
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
      
      sumgrp <- cmp %>% 
        filter(var %in% c('a', 'Rt_vol', 'b')) %>% 
        group_by(grp, var) %>% 
        summarise(
          Fwoxy = mean(Fwoxy, na.rm = T), 
          EBASE = mean(EBASE, na.rm = T), 
          .groups = 'drop'
        ) 
      sumtms <- cmp %>% 
        filter(!var %in% c('a', 'Rt_vol', 'b'))
      
      sumcmp <- bind_rows(sumgrp, sumtms) %>% 
        group_by(var) %>% 
        nest() %>% 
        summarise(
          est = purrr::map(data, sumfun)
        ) %>% 
        unnest('est')
      
      return(sumcmp)
      
    })
  ) %>% 
  select(-out)

save(apasumdat, file = here('data/apasumdat.RData'), compress = 'xz')


