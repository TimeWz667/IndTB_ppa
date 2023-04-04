library(tidyverse)
library(tidybayes)

theme_set(theme_bw())


post <- bind_rows(lapply(c("dx_00", "dx_01", "dx_10", "dx_11"), function(folder) {
  read_csv(here::here("out", folder, "post.csv")) %>% 
    mutate(Scenario = folder)
}))



post %>% 
  select(Scenario, ppv_pri, pdx_pri, dur_pri, p_pri_on_pub, p_under) %>% 
  pivot_longer(-Scenario, names_to = "Index") %>% 
  ggplot() +
  stat_halfeyeh(aes(x = value, y = Index, colour = Scenario, fill = Scenario), alpha = 0.3)



post %>% 
  select(Scenario, ppv_pri, pdx_pri, dur_pri, p_pri_on_pub, p_under) %>% 
  pivot_longer(-Scenario, names_to = "Index") %>% 
  group_by(Scenario, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  ggplot() + 
  geom_pointrange(aes(x = M, xmin = L, xmax = U, y = Index, colour = Scenario), 
                  position = position_dodge(0.3))




post <- bind_rows(lapply(c("tx_00", "tx_01", "tx_10", "tx_11"), function(folder) {
  read_csv(here::here("out", folder, "post.csv")) %>% 
    mutate(Scenario = folder)
}))



post %>% 
  select(Scenario, ppv_pri, dur_pri, p_pri_on_pub, p_under) %>% 
  pivot_longer(-Scenario, names_to = "Index") %>% 
  ggplot() +
  stat_halfeye(aes(x = value, y = Index, colour = Scenario, fill = Scenario), alpha = 0.3)



post %>% 
  select(Scenario, ppv_pri, dur_pri, p_pri_on_pub, p_under) %>% 
  pivot_longer(-Scenario, names_to = "Index") %>% 
  group_by(Scenario, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  ) %>% 
  ggplot() + 
  geom_pointrange(aes(x = M, xmin = L, xmax = U, y = Index, colour = Scenario), 
                  position = position_dodge(0.3))




post %>% 
  filter(Scenario %in% c("tx_01", "tx_10", "tx_11")) %>% 
  select(Scenario, tp_pri_txi, tp_pri_drug_time) %>% 
  pivot_longer(-Scenario, names_to = "Index") %>% 
  group_by(Scenario, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  ) %>% 
  ggplot() + 
  geom_pointrange(aes(x = M, xmin = L, xmax = U, y = Index, colour = Scenario), 
                  position = position_dodge(0.3))


tp_pri_txi



post %>% 
  filter(Scenario %in% c("tx_01", "tx_10", "tx_11")) %>% 
  select(Scenario, ppv_pri, dur_pri) %>% 
  ggplot() +
  geom_point(aes(x = ppv_pri, y = dur_pri, colour = Scenario)) +
  facet_grid(.~Scenario)


