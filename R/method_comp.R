library(tidyverse)
library(lubridate)
library(WtRegDO)

# compare Bayes with Odum -------------------------------------------------

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
  facet_wrap(~typ, ncol = 1, scales = 'free') + 
  theme(
    strip.background = element_blank()
  ) + 
  scale_y_continuous(limits= c(-150, 150))


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

