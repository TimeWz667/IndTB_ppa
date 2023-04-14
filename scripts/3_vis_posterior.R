library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


for (sc in c("cs_11", "dx_11", "tx_11", "tx_01", "tx_10", "tx_00")) {
  
  post <- read_csv(here::here("out", sc, "post.csv"))
  
  
  g_ppm <- post %>% 
    ggplot() + 
    stat_halfeye(aes(x = ppm), alpha = 0.7) + 
    scale_x_continuous("Engaged private visits / all private visits, %", labels = scales::percent_format()) +
    scale_y_continuous("") + 
    expand_limits(x = c(0, 1)) +
    theme(axis.text.y = element_blank())
  
  g_ppm
  
  # g_pdx_pri <- post %>% 
  #   ggplot() + 
  #   stat_halfeye(aes(x = pdx_eng), alpha = 0.7) + 
  #   scale_x_continuous("TB case detected per private visits, %", labels = scales::percent_format()) +
  #   scale_y_continuous("") + 
  #   expand_limits(x = c(0, 1)) +
  #   theme(axis.text.y = element_blank())
  # 
  # g_pdx_pri
  
  
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
  
  
  g_bind <- ggpubr::ggarrange(g_ppm, g_ppv_pri, g_dur_pri, ncol=1)
  
  g_bind
  
  ggsave(g_bind, filename = here::here("out", sc, "g_posterior.png"), width = 5, height = 5)

}





post <- bind_rows(lapply(c("tx_00", "tx_01", "tx_11"), function(folder) {
  read_csv(here::here("out", folder, "post.csv")) %>% 
    select(ppm, dur_pri, ppv_pri, p_pri_on_pub, 
           tp_pri_drug, tp_pri_drug_time, tp_pri_txi) %>% 
    mutate(Scenario = folder)
})) %>% 
  mutate(
    Scenario = case_when(
      Scenario == "tx_00" ~ "Prior distribution",
      Scenario == "tx_01" ~ "Drug sales data alone",
      T ~ "Drug sales + prevalence survey data"
    )
  )



g_ppv <- post %>% 
  filter(Scenario != "Prior distribution") %>% 
  ggplot() +
  geom_rect(xmin=0.2, xmax=0.85, ymin=0, ymax=1 / (0.85 - 0.2), alpha = 0.01, 
            colour = "black", linetype = 2,
            aes(fill = "Prior distribution")) +
  geom_density(aes(x = ppv_pri, fill = Scenario), alpha = 0.4, position = "identity") +
  scale_x_continuous("Positive predictive value, private sector, %", labels = scales::percent) +
  scale_y_continuous("Probability density") +
  scale_fill_discrete("") +
  guides(fill = guide_legend(reverse = T))




g_dur <- post %>% 
  filter(Scenario != "Prior distribution") %>% 
  ggplot() +
  geom_rect(xmin=1 / 12, xmax=2, ymin=0, ymax=1 / (2 - 0.1), alpha = 0.01, 
            colour = "black", linetype = 2,
            aes(fill = "Prior distribution")) +
  geom_density(aes(x = dur_pri, fill = Scenario), alpha = 0.4, position = "identity") +
  scale_x_continuous("Treatment period, private providers, month", labels = scales::number_format(scale = 12)) +
  scale_y_continuous("Probability density") +
  scale_fill_discrete("") +
  guides(fill = guide_legend(reverse = T))


g_ppv
g_dur



g_post <- ggpubr::ggarrange(g_ppv, g_dur, common.legend = T, legend = "bottom")



g_cross <- post %>% 
  filter(Scenario != "Prior distribution") %>% 
  ggplot() +
  geom_rect(xmin=1 / 12, xmax=2, ymin=0.2, ymax=0.85, alpha = 0.01, 
            linetype = 2, fill = NA,
            aes(colour = "Prior distribution")) +
  geom_point(aes(x = dur_pri, y = ppv_pri, colour = Scenario), alpha = 0.4, position = "identity") +
  scale_x_continuous("Treatment period, private providers, month", labels = scales::number_format(scale = 12))  +
  scale_y_continuous("Positive predictive value, private sector, %", labels = scales::percent) +
  scale_color_discrete("")


g_cross

ggsave(g_ppv, filename = here::here("docs", "figs", "g_post_ppv.png"), width = 6, height = 4)

ggsave(g_dur, filename = here::here("docs", "figs", "g_post_dur.png"), width = 6, height = 4)

ggsave(g_post, filename = here::here("docs", "figs", "g_post.png"), width = 9, height = 6)

ggsave(g_cross, filename = here::here("docs", "figs", "g_post_cross.png"), width = 9, height = 6)



