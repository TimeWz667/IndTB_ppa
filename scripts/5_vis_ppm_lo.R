library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))



ds <- dir(here::here("out", "sens_ppm_lo"))
ds <- ds[startsWith(ds, "post_")]


post <- bind_rows(lapply(ds, function(file) {
  read_csv(here::here("out", "sens_ppm_lo", file)) %>% 
    select(ppm, dur_pri, ppv_pri, p_pri_on_pub, 
           tp_pri_drug, tp_pri_drug_time, tp_pri_txi) %>% 
    mutate(Scenario = file)
})) %>% 
  extract(Scenario, "PPM_lb", "post_(\\S+).csv", convert=T)



g_ppm <- post %>% 
  ggplot() + 
  stat_halfeye(aes(x = ppm, fill = as.factor(PPM_lb)), alpha = 0.7) + 
  scale_x_continuous("Unengaged private Txi / all private Txi, %", labels = scales::percent_format()) +
  scale_y_continuous("") + 
  scale_fill_discrete("LB of PPM") +
  expand_limits(x = c(0, 1)) +
  theme(axis.text.y = element_blank())


g_ppv <- post %>% 
  ggplot() + 
  stat_halfeye(aes(x = ppv_pri, fill = as.factor(PPM_lb)), alpha = 0.7) + 
  scale_x_continuous("Positive Predictive Value, in the private sector, %", labels = scales::percent_format()) +
  scale_y_continuous("") + 
  scale_fill_discrete("LB of PPM") +
  expand_limits(x = c(0, 1)) +
  theme(axis.text.y = element_blank())



g_dur <- post %>% 
  ggplot() + 
  stat_halfeye(aes(x = dur_pri, fill = as.factor(PPM_lb)), alpha = 0.7) + 
  scale_x_continuous("Treatment duration, on private drug, month", labels = scales::number_format(scale = 12)) +
  scale_y_continuous("") + 
  scale_fill_discrete("LB of PPM") +
  expand_limits(x = c(0, 1)) +
  theme(axis.text.y = element_blank())


g_bind <- ggpubr::ggarrange(g_ppv, g_dur, legend = "bottom", nrow=1, common.legend = T)
ggsave(g_bind, filename = here::here("docs", "figs", "g_sens_ppmlb.png"), width = 9, height = 5)


