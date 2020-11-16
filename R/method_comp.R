library(tidyverse)
library(lubridate)
library(WtRegDO)
library(patchwork)

# compare Bayes with Odum -------------------------------------------------

# outputpar is from script.R
data(outputpar)

bsmetab <- outputpar %>% 
  select(
    Pg = GPP, 
    Rt = ER, 
    NEM = NEP,
    Date
  ) %>% 
  mutate(
    Rt = -1 * Rt,
    Date = ymd(Date) 
  ) %>% 
  gather('var' ,'val', -Date) %>% 
  mutate(typ = 'Bayes')

load(file = 'data/DO_APNERR2012_6_12_0.8.RData')
tmp <- DO_APNERR2012_6_12_0.8 %>%
  select(DateTimeStamp, DO_mgl = DO_nrm, Temp, Sal, ATemp, BP, WSpd, totpar, Tide)

opmetab <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = T) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum')

toplo <- bind_rows(bsmetab, opmetab)

ggplot(toplo, aes(x = Date, y = val, colour = var, fill = var)) + 
  geom_line() + 
  geom_point() + 
  theme_minimal() + 
  facet_wrap(~typ, ncol = 1) + 
  theme(
    strip.background = element_blank()
  ) + 
  scale_y_continuous(limits= c(-300, 300))

# compare different gas exchange for Odum ---------------------------------

load(file = 'data/DO_APNERR2012_6_12_0.8.RData')
tmp <- DO_APNERR2012_6_12_0.8 %>%
  select(DateTimeStamp, DO_mgl = DO_nrm, Temp, Sal, ATemp, BP, WSpd, totpar, Tide)

opmetab1 <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Thiebault', gasave = F) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum - Thiebault continuous')

opmetab2 <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Thiebault', gasave = T) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum - Thiebault daily average')

opmetab3 <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = F) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum - Wanninkhof continuous')

opmetab4 <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = T) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'Odum - Wanninkhof daily average')

toplo <- bind_rows(opmetab1, opmetab2, opmetab3, opmetab4)

ggplot(toplo, aes(x = Date, y = val, colour = var, fill = var)) + 
  geom_line() + 
  geom_point() + 
  theme_minimal() + 
  facet_wrap(~typ, ncol = 1, scales = 'free') + 
  theme(
    strip.background = element_blank()
  )# + 
# scale_y_continuous(limits= c(-150, 150))


# compare original BASEmetb with Odum -------------------------------------

# BASE results
bsres <- read.csv('results/BASE_results_2020-11-16 101039.csv', stringsAsFactors = F) %>% 
  # filter(R.Rhat < 1.1) %>%
  # filter(PPP > 0.4) %>%
  mutate(ER.mean = -1 * ER.mean) %>% 
  select(Date, Pg = GPP.mean, Rt = ER.mean, NEM = NEP.mean) %>% 
  mutate(Date = lubridate::ymd(Date)) %>% 
  gather('var', 'val', -Date) %>% 
  mutate(typ = 'BASE')

# Odum results
load(file = 'data/DO_APNERR2012_6_12_0.8.RData')
tmp <- DO_APNERR2012_6_12_0.8 %>%
  select(DateTimeStamp, DO_mgl = DO_nrm, Temp, Sal, ATemp, BP, WSpd, totpar, Tide)

# get metab as g/m2/d, convert to mg/L/d (multiply by  to get mg, divide by H to get m3, divide by 1000 to get L)
opres <- ecometab(tmp, tz = 'America/Jamaica', lat = 29.75, long = -85, gasex = 'Wanninkhof', gasave = T, metab_units = 'grams') %>% 
  gather('var', 'val', -Date) %>% 
  mutate(val = val / 2.08) %>% 
  mutate(typ = 'Odum')

cols <- list('#F8766D', '#00BA38', '#619CFF') 
names(cols) <- c('NEM', 'Pg', 'Rt')
alph <- 0.8

toplo1 <- bind_rows(bsres, opres) %>% 
  mutate(var = factor(var, levels = c('NEM', 'Pg', 'Rt')))

p1 <- ggplot(toplo1, aes(x = Date, y = val, colour = var)) + 
  geom_line(alpha = alph) + 
  geom_point(alpha = alph) + 
  theme_minimal() + 
  facet_wrap(~typ, ncol = 1) + 
  scale_x_date(date_labels = '%b', date_breaks = '1 month') +
  # scale_colour_manual(values = cols) +
  theme(
    strip.background = element_blank(), 
    legend.title = element_blank(), 
    axis.title.x = element_blank(), 
    axis.ticks.x = element_line(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.x = element_blank()
  ) + 
  labs(
    y = expression(paste('mg ', O [2], ' ', L^-1, d^-1))
  )
p1

toplo2 <- toplo1 %>% 
  pivot_wider(names_from = typ, values_from = val)

sel <- 'NEM'
p2a <- ggplot(toplo2[toplo2$var == sel, ], aes(x = BASE, y = Odum, colour = var, fill = var)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(alpha = alph) + 
  geom_smooth(method = 'lm', se = F, colour = 'grey') +
  scale_colour_manual(values = cols) +
  theme_minimal() + 
  theme(
    legend.position = 'none', 
    axis.title.x = element_blank()
  ) + 
  labs(
    subtitle = sel
  )
sel <- 'Pg'
p2b <- ggplot(toplo2[toplo2$var == sel, ], aes(x = BASE, y = Odum, colour = var, fill = var)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(alpha = alph) + 
  geom_smooth(method = 'lm', se = F, colour = 'grey') +
  scale_colour_manual(values = cols, drop = F) +
  theme_minimal() + 
  theme(
    legend.position = 'none', 
    axis.title.y = element_blank()
  ) + 
  labs(
    subtitle = sel
  )
sel <- 'Rt'
p2c <- ggplot(toplo2[toplo2$var == sel, ], aes(x = BASE, y = Odum, colour = var, fill = var)) + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_point(alpha = alph) + 
  geom_smooth(method = 'lm', se = F, colour = 'grey') +
  scale_colour_manual(values = cols, drop = F) +
  theme_minimal() + 
  theme(
    legend.position = 'none', 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  ) + 
  labs(
    subtitle = sel
  )

p2 <- p2a + p2b + p2c + plot_layout(ncol = 3)

p <- p1 + p2 + plot_layout(ncol = 1)
png('~/Desktop/methodcomp.png', height = 7, width = 8, units = 'in', res = 300)
print(p)
dev.off()
