library(fwoxy)
library(tidyverse)
library(lubridate)
# library(WtRegDO)
# devtools::load_all('../fwoxy')
devtools::load_all('../WtRegDO')

# Set model parameters
oxy_ic <- 250           # (mmol/m^3), initial oxygen concentration
a_param <- 0.2          # ((mmol/m^3)/day)/(W/m^2), light efficiency
er_param <- 20          # (mmol/m^3/day), ecosystem respiration

# Constant Forcings
ht_const <- 3           # m, height of the water column
salt_const <- 25        # ppt, salinity
temp_const <- 25        # deg C, water temperature
wspd_const <- 3         # m/s, wind speed at 10 m

example <- fwoxy(oxy_ic = oxy_ic, a_param = a_param, er_param = er_param, 
                 ht_in = ht_const, salt_in = salt_const, temp_in = temp_const,
                 wspd_in = wspd_const)

tomod <- example %>% 
  mutate(
    DateTimeStamp = force_tz(as.POSIXct(`time, sec`, origin = Sys.Date(), tz = 'UTC'), tzone = 'America/Jamaica'),
    DO_obs = `oxy, mmol/m3` * 0.032, # to mg/L
    Temp = temp_const, 
    WSpd = wspd_const,
    Sal = salt_const, 
    ATemp = NA, 
    BP = 1013.25, 
    Tide = NA
  ) %>% 
  select(DateTimeStamp, Temp, Sal, DO_obs, WSpd, ATemp, BP, Tide)

getsec <- function(x){
  
  sec <- seconds(x)
  sec <- sec - (7.5 * 60)
  sec <- sec - min(sec)
  sec <- 1 + as.numeric(sec)
  
  return(sec)
    
}

# estimate instantaneous Odum metabolism from Fwoxy output
opmetab <- ecometab(tomod, DO_var = 'DO_obs', tz = 'America/Jamaica', lat = 29.75, long = -85, depth_val = NULL, 
                    depth_vec = ht_const, gasex = 'Wanninkhof', gasave = 'instant', metab_units = 'mmol', instant = T) %>% 
  select(DateTimeStamp, dDO, D, Pg_vol, Rt_vol) %>% 
  mutate(
    secs = getsec(DateTimeStamp),
    D = -1 * D, 
    Rt_vol = -1 * Rt_vol, 
    typ = 'Odum'
  ) %>% 
  select(-DateTimeStamp) %>% 
  gather(var, val, -secs, -typ)

toplo1 <- example %>%
  select(
    secs = `time, sec`, 
    dDO = `troc, mmol/m3/d`, 
    D = `gasex, mmol/m3/d`,
    Pg_vol = `gpp, mmol/m3/d`,
    Rt_vol = `er, mmol/m3/d`
  ) %>% 
  mutate(
    typ = 'Fwoxy'
  ) %>% 
  gather(var, val, -secs, -typ) %>% 
  bind_rows(opmetab)

p1 <- ggplot(toplo1, aes(x = secs, y = val, color = var)) + 
  geom_line() +
  facet_wrap(~typ, ncol = 1) +
  theme_bw() + 
  labs(
    y = 'Flux, mmol/m3/d', 
    x = 'seconds'
  ) +
  theme(
    legend.title = element_blank(), 
    strip.background = element_blank()
  )
p1

toplo2 <- toplo1 %>% 
  spread(typ, val)
p2 <- ggplot(toplo2, aes(x = Fwoxy, y = Odum, color = var)) + 
  geom_point() +
  facet_wrap(~var, ncol = 2, scales = 'free') +
  theme_bw() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(
    y = 'Odum, Flux, mmol/m3/d', 
    x = 'Fwoxy, flux, mmol/m3/d'
  ) +
  theme(
    legend.title = element_blank(), 
    strip.background = element_blank()
  )
p2

p3 <- ggplot(toplo1, aes(x = secs, y = val, color = typ)) + 
  geom_line() +
  facet_wrap(~var, ncol = 1, scales = 'free') +
  theme_bw() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(
    y = 'Odum, Flux, mmol/m3/d', 
    x = 'Fwoxy, flux, mmol/m3/d'
  ) +
  theme(
    legend.title = element_blank(), 
    strip.background = element_blank()
  )
p3
