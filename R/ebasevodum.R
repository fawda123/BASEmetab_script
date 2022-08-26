library(WtRegDO)
library(EBASE)
library(tidyverse)
library(lubridate)
library(doParallel)
library(here)

H <-  1.85
tz <- 'America/Jamaica'
dtrng <- ymd('2012-07-01', '2012-07-07')

# data prep -----------------------------------------------------------------------------------

load(file = 'data/DO_APNERR2012_6_12_0.8.RData')

# DO as mg/L
# Temp (water) as deg C
# Salinity as ppt
# ATemp as deg C
# BP as mb
# WSpd as m/s
# Tide as m (but not used)
opdat <- DO_APNERR2012_6_12_0.8 %>%
  select(DateTimeStamp, DO_mgl = DO_obs, Temp, Sal, ATemp, BP, WSpd, Tide)

# totpar is mmol/m2 total for 15 minute obs, umol/m2/s
# sal is ppt, should be ppt
# DO is mg/l, should be mg/l
# temp is C, should be C
# BP is mb, should be atm
# WSpd is m/s, should be m/s
# last line is to fill NA (only 3 values) with last value
ebdat <- DO_APNERR2012_6_12_0.8 %>% 
  mutate(
    totpar = ifelse(is.na(totpar), 0, totpar), 
    totpar = totpar * 1000 / (15 * 60), # convert to umol and per second
  ) %>% 
  select(DateTimeStamp, DO_obs, Temp, Sal, PAR = totpar, WSpd) 

# estimate metabolism -------------------------------------------------------------------------

# odum
opmetab <- ecometab(opdat, tz = tz, lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = 'instant', metab_units = 'mmol', depth_vec = H, depth_val = NULL, replacemet = FALSE, instant = T)

# # ebase, takes about an hour
# cl <- makeCluster(10)
# registerDoParallel(cl)
# ebmetab <- ebase(ebdat, H = H, interval = 900, ndays = 7, n.chains = 5, progress = T)
# stopCluster(cl)
# 
# save(ebmetab, file = here('data/ebmetab.RData'))

load(file = here('data/ebmetab.RData'))

# gas exchange plots --------------------------------------------------------------------------

toplo1 <- opmetab %>% 
  select(DateTimeStamp, D_op = D) %>% 
  mutate(DateTimeStamp = DateTimeStamp + (60 * 7.5)) %>% 
  left_join(ebmetab, by = 'DateTimeStamp') %>% 
  select(DateTimeStamp, D_op, D_eb = D) %>% 
  mutate(
    month = month(DateTimeStamp, label = T, abbr = T)
  ) %>% 
  group_by(month) %>% 
  mutate(
    r2 = summary(lm(D_eb ~ D_op))$r.squared, 
    r2 = paste0(round(100 * r2, 0), '%')
  ) %>% 
  unite('month', month, r2, sep = ', ')

toplo2 <- toplo1 %>% 
  pivot_longer(-matches('DateTimeStamp|month'))

ylab <- expression(paste(O [2], ' (mmol ', m^-3, ' ', d^-1, ')'))

p1 <- ggplot(toplo2, aes(x = DateTimeStamp, y = value, color = name)) + 
  geom_line() + 
  theme_minimal() + 
  facet_wrap(. ~ month, scales = 'free', ncol = 2) + 
  theme(
    strip.background = element_blank(), 
    strip.text = element_blank(), 
    axis.text = element_text(size = 8), 
    legend.position = c(0.75, 0.05), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()
  ) + 
  scale_color_discrete(
    breaks = c('D_eb', 'D_op'),
    labels = c(expression(D[EBASE]), expression(D[Odum]))
  ) +
  scale_x_datetime(expand = c(0, 0)) + 
  labs(
    x = NULL, 
    y = ylab, 
    color = NULL
    )

p2 <- ggplot(toplo1, aes(x = D_op, y = D_eb)) + 
  geom_point(alpha = 0.5, show.legend = F) + 
  stat_smooth(formula = y ~ x, aes(alpha = ''), method = 'lm', se = F, color = 'red') +
  geom_abline(aes(color = '', intercept = 0, slope = 1), size = 0.8, linetype = 'dotted') +
  facet_wrap(~month, ncol = 4, scales = 'free') +
  scale_color_manual(values = 'red') +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    axis.text = element_text(size = 8), 
    legend.position = c(0.875, 0.12), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()
  ) +
  labs(
    x = expression(D[Odum]), 
    y = expression(D[EBASE]), 
    alpha = 'Regression',
    color = "1:1", 
    caption = 'R-squared values right of month'
  )

jpeg(here('docs/images/gasexts.jpg'), height = 5, width = 9, family = 'serif', units = 'in', res = 300)
print(p1)
dev.off()

jpeg(here('docs/images/gasexcmp.jpg'), height = 5, width = 9, family = 'serif', units = 'in', res = 300)
print(p2)
dev.off()

# weekly averages -----------------------------------------------------------------------------

toplo1 <- opmetab %>% 
  mutate(
    week = floor_date(DateTimeStamp, unit = 'weeks')
  ) %>% 
  select(week, Pg_vol, Rt_vol, NEM, D) %>% 
  tidyr::pivot_longer(cols = -matches('week'), values_to = 'opval') 
  
toplo2 <- ebmetab %>% 
  mutate(
    week = floor_date(DateTimeStamp, unit = 'weeks'),
    Rt_vol = -1 * Rt_vol,
    NEM = Pg_vol + Rt_vol
  ) %>% 
  select(week, Pg_vol, Rt_vol, NEM, D) %>%  
  tidyr::pivot_longer(cols = -matches('week'), values_to = 'ebval')
  
toploa <- full_join(toplo1, toplo2, by = c('week', 'name')) %>% 
  group_by(week, name) %>% 
  summarise(
    opvalue = mean(opval, na.rm = T), 
    ebvalue = mean(ebval, na.rm = T), 
    .groups = 'drop'
  ) %>% 
  mutate(
    name = factor(name, levels = c('NEM', 'Pg_vol', 'Rt_vol', 'D'), labels = c('NEM', 'Pg', 'Rt', 'D'))
  ) %>% 
  group_by(name) %>% 
  mutate(
    r2 = summary(lm(ebvalue ~ opvalue))$r.squared, 
    r2 = paste0(round(100 * r2, 0), '%')
  ) %>% 
  unite('name', name, r2, sep = ', ')

toplob <- toploa %>% 
  pivot_longer(matches('opvalue|ebvalue'), names_to = 'var') %>% 
  mutate(var = factor(var, levels = c('ebvalue', 'opvalue'), labels = c('EBASE', 'Odum')))

p1 <- ggplot(toploa, aes(x = opvalue, y = ebvalue)) + 
  geom_point() + 
  facet_wrap(~name, scales = 'free') + 
  stat_smooth(formula = y ~ x, aes(alpha = ''), method = 'lm', se = F, color = 'red') +
  geom_abline(aes(color = '', intercept = 0, slope = 1), size = 0.8, linetype = 'dotted') +
  scale_color_manual(values = 'red') +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    axis.text = element_text(size = 8), 
    legend.position = 'top', 
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()
  ) + 
  labs(
    x = 'Odum', 
    y = 'EBASE', 
    alpha = 'Regression',
    color = "1:1", 
    caption = 'R-squared values right of parameter'
  )

p2 <- ggplot(toplob, aes(x = week, y = value, group = var, color = var)) + 
  geom_line() +
  geom_point() + 
  facet_wrap(~name, scales = 'free', ncol = 2) +
  theme_minimal() +
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(hjust = 0), 
    axis.text = element_text(size = 8), 
    legend.position = 'top', 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()
  ) + 
  scale_x_datetime(expand = c(0, 0)) + 
  labs(
    x = NULL, 
    y = ylab, 
    color = NULL,
    caption = 'R-squared values right of parameter'
  )

jpeg(here('docs/images/weeklycmp.jpg'), height = 5, width = 7, family = 'serif', units = 'in', res = 300)
print(p1)
dev.off()

jpeg(here('docs/images/weeklyts.jpg'), height = 5, width = 7, family = 'serif', units = 'in', res = 300)
print(p2)
dev.off()

# hourly averages by quarter ------------------------------------------------------------------

toplo1 <- opmetab %>% 
  mutate(
    quarter = quarter(DateTimeStamp), 
    hour = hour(DateTimeStamp)
  ) %>% 
  select(quarter, hour, Pg_vol, Rt_vol, NEM, D) %>% 
  tidyr::pivot_longer(cols = -matches('quarter|hour'), values_to = 'opval') %>% 
  group_by(quarter, hour, name) %>% 
  summarise(
    opval = mean(opval, na.rm = T), 
    .groups = 'drop'
  )

toplo2 <- ebmetab %>% 
  mutate(
    quarter = quarter(DateTimeStamp), 
    hour = hour(DateTimeStamp),
    Rt_vol = -1 * Rt_vol,
    NEM = Pg_vol + Rt_vol
  ) %>% 
  select(quarter, hour, Pg_vol, Rt_vol, NEM, D) %>%  
  tidyr::pivot_longer(cols = -matches('quarter|hour'), values_to = 'ebval') %>% 
  group_by(quarter, hour, name) %>% 
  summarise(
    ebval = mean(ebval, na.rm = T), 
    .groups = 'drop'
  )

toploa <- full_join(toplo1, toplo2, by = c('quarter', 'hour', 'name')) %>% 
  mutate(
    name = factor(name, levels = c('NEM', 'Pg_vol', 'Rt_vol', 'D'), labels = c('NEM', 'Pg', 'Rt', 'D')), 
    quarter = factor(quarter, levels = c('1', '2', '3', '4'), labels = c('JFM', 'AMJ', 'JAS', 'OND'))
  ) %>% 
  pivot_longer(matches('ebval|opval'), names_to = 'var') %>% 
  mutate(
    var = factor(var, levels = c('ebval', 'opval'), labels = c('EBASE', 'Odum'))
  )

p1 <- ggplot(toploa, aes(x = hour, y = value, group = var, color = var)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(name ~ quarter, scales = 'free_y') +
  theme_minimal() +
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(hjust = 0), 
    axis.text = element_text(size = 8), 
    legend.position = 'top', 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.minor.y = element_blank()
  ) + 
  scale_x_continuous(expand = c(0, 0)) + 
  labs(
    x = 'Hour of day', 
    y = ylab, 
    color = NULL
  )

jpeg(here('docs/images/hourts.jpg'), height = 5, width = 9, family = 'serif', units = 'in', res = 300)
print(p1)
dev.off()
