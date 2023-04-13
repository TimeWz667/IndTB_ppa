library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


post <- bind_rows(lapply(c("tx_00", "tx_01", "tx_11", "dx_11", "cs_11"), function(folder) {
  read_csv(here::here("out", folder, "post.csv")) %>% 
    select(ppm, dur_pri, ppv_pri, p_pri_on_pub, p_under,
           tp_pri_drug, tp_pri_drug_time, tp_pri_txi) %>% 
    mutate(Scenario = folder)
})) %>% 
  mutate(
    Scenario = factor(Scenario, c("tx_00", "tx_01", "tx_11", "dx_11", "cs_11"))
  )




post %>% 
  ggplot() +
  stat_halfeye(aes(x = ppm, y = Scenario)) 


post %>% 
  ggplot() +
  stat_halfeye(aes(x = tp_pri_txi, y = Scenario)) 



g_p_under <- post %>% 
  filter(Scenario %in% c("tx_11", "dx_11", "cs_11")) %>% 
  ggplot() +
  stat_halfeye(aes(x = p_under, y = Scenario)) +
  scale_x_continuous("Proportion cases not notified, %", labels = scales::percent) +
  scale_y_discrete("Modelling from ", labels = c(tx_11="Treatment", dx_11="Diagnosis", cs_11="Care-seeking")) + 
  expand_limits(x=0)
  
g_p_under


ggsave(g_p_under, filename=here::here("docs", "figs", "incr_p_under.png"), width=4, height=5)



post %>% 
  ggplot() +
  geom_point(aes(x = dur_pri, y = p_under)) +
  facet_wrap(.~Scenario)



post %>% 
  ggplot() +
  geom_point(aes(x = dur_pri, y = ppv_pri)) +
  facet_wrap(.~Scenario)


