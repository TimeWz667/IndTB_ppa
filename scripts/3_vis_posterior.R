library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


post <- read_csv(here::here("out", "cs_11", "post.csv"))


g_ppm <- post %>% 
  ggplot() + 
  stat_halfeye(aes(x = ppm), alpha = 0.7) + 
  scale_x_continuous("Engaged private visits / all private visits, %", labels = scales::percent_format()) +
  scale_y_continuous("") + 
  expand_limits(x = c(0, 1)) +
  theme(axis.text.y = element_blank())

g_ppm

g_pdx_pri <- post %>% 
  ggplot() + 
  stat_halfeye(aes(x = pdx_eng), alpha = 0.7) + 
  scale_x_continuous("TB case detected per private visits, %", labels = scales::percent_format()) +
  scale_y_continuous("") + 
  expand_limits(x = c(0, 1)) +
  theme(axis.text.y = element_blank())

g_pdx_pri


g_dur_pri <- post %>% 
  ggplot() + 
  stat_halfeye(aes(x = dur_pri), alpha = 0.7) + 
  scale_x_continuous("Treatment duration, on private drug, month", labels = scales::number_format(scale = 12)) +
  scale_y_continuous("") + 
  expand_limits(x = 0) +
  theme(axis.text.y = element_blank())


g_dur_pri

g_ppv_pri <- post %>% 
  ggplot() + 
  stat_halfeye(aes(x = ppv_pri), alpha = 0.7) + 
  scale_x_continuous("Positive Predictive Value, in the private sector, %", labels = scales::percent_format()) +
  scale_y_continuous("") + 
  expand_limits(x = c(0, 1)) +
  theme(axis.text.y = element_blank())


g_bind <- ggpubr::ggarrange(g_ppm, g_pdx_pri, g_ppv_pri, g_dur_pri, ncol=1)

g_bind

ggsave(g_bind, filename = here::here("out", "cs_11", "g_posterior.png"), width = 5, height = 5)









