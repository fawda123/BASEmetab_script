library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(here)
library(EBASE)
library(doParallel)

fwdat <- read_csv(here("data/apafwoxy.csv")) 

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

# subset four days in June
dat <- fwdatinp %>%
  filter(month(fwdatinp$DateTimeStamp) == 6 & day(fwdatinp$DateTimeStamp) %in% 1:4)

cl <- makeCluster(4)
registerDoParallel(cl)

res <- ebase(dat, interval = 900, H = dat$H, progress = TRUE, n.chains = 4, interp = F)

stopCluster(cl)

fwdatinpcmp <- fwdatcmp %>% 
  filter(Date <= max(res$Date) & Date >= min(res$Date))

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

ggplot(cmp, aes(x = Fwoxy, y = EBASE, color = var)) + 
  geom_point() +
  facet_wrap(~var, ncol = 2, scales = 'free') +
  # coord_equal() + 
  # scale_color_manual(values = colors) + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = 1) +
  # labs(
  #   y = 'Odum, Flux, mmol/m3/d', 
  #   x = 'Fwoxy, flux, mmol/m3/d'
  # ) +
  theme(
    legend.title = element_blank(),
    strip.background = element_blank()
  )