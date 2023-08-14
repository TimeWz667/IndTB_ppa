library(tidyverse)
library(tidybayes)

theme_set(theme_bw())



cas <- read_csv(here::here("docs", "tabs", "cascade.csv"))

locations <- read_csv(here::here("data", "locations.csv"))


g_inc <- cas %>% 
  mutate(
    Location = reorder(Location, Inc)
  ) %>% 
  ggplot() +
  stat_halfeye(aes(x = Inc, y = Location)) +
  scale_x_continuous("Incidence, per 100 000", labels = scales::number_format(scale = 1e5)) +
  expand_limits(x = 0)


ggsave(g_inc, filename = here::here("docs", "figs", "g_inc.png"), width = 7, height = 6)



g_delay_pat <- cas %>% 
  mutate(
    Location = reorder(Location, DelayPat)
  ) %>% 
  ggplot() +
  stat_halfeye(aes(x = DelayPat, y = Location)) +
  scale_x_continuous("Patient delay, months", labels = scales::number_format(scale = 30)) +
  expand_limits(x = 0)

g_delay_pat



g_delay_sys <- cas %>% 
  mutate(
    Location = reorder(Location, DelaySys)
  ) %>% 
  ggplot() +
  stat_halfeye(aes(x = DelaySys, y = Location)) +
  scale_x_continuous("System delay, months", labels = scales::number_format(scale = 30)) +
  expand_limits(x = 0)

g_delay_sys


g_delay <- ggpubr::ggarrange(g_delay_pat, g_delay_sys)


ggsave(g_delay, filename = here::here("docs", "figs", "g_delay.png"), width = 10, height = 6)


