library(tidyverse)
library(tidybayes)

theme_set(theme_bw())


post <- bind_rows(lapply(c("dx_00", "dx_01", "dx_10", "dx_11", "tx_00", "tx_01", "tx_10", "tx_11"), function(folder) {
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



post %>% 
  filter(Scenario %in% c("tx_01", "dx_11")) %>% 
  select(Scenario, ppv_pri, dur_pri) %>% 
  mutate(
    Scenario = ifelse(Scenario == "tx_01", "Without TBPS data", "With TBPS data")
  ) %>% 
  ggplot() +
  geom_point(aes(y = ppv_pri, x = dur_pri), alpha = 0.4) +
  scale_y_continuous("PPV, private diagnosis, %", labels = scales::percent_format()) + 
  scale_x_continuous("Treatment duration, treated privately, months", labels = scales::number_format(scale = 12)) +
  facet_grid(.~Scenario)


g_txi_pri2 <- post %>% 
  filter(Scenario %in% c("tx_01", "dx_11")) %>% 
  select(Scenario, tp_pri_txi) %>% 
  mutate(
    Scenario = ifelse(Scenario == "tx_01", "Drug sales data alone", "Drug sales + prevalence survey data"),
    Scenario = factor(Scenario, c("Drug sales data alone", "Drug sales + prevalence survey data"))
  ) %>% 
  ggplot() +
  geom_density(aes(x = tp_pri_txi, fill = Scenario), alpha = 0.4) +
  scale_x_continuous("Annual number of patients with TB initiating treatment in private sector (million)", 
                     labels = scales::number_format(scale = 1e-6), limits = c(0, 3e6)) +
  scale_y_continuous("Probability density") +
  scale_fill_discrete("") +
  theme(legend.position = c(1, 1), legend.justification = c(1.05, 1.2), axis.text.y = element_blank())


g_txi_pri2

ggsave(g_txi_pri2, filename = here::here("docs", "figs", "txi_pri2.png"), width = 6, height = 4)




sub <- post %>% 
  filter(Scenario %in% c("tx_01", "tx_11","dx_11")) %>% 
  mutate(
    Scenario = case_when(
      Scenario == "tx_01" ~ "Without TBPS", 
      Scenario == "tx_11" ~ "With TBPS:Tx",
      T ~ "With TBPS: Tx+CS"
    ),
    Scenario = factor(Scenario, c("Without TBPS", "With TBPS:Tx", "With TBPS: Tx+CS"))
  )



sub %>% 
  ggplot() +
  geom_density(aes(x = ppv_pri, fill = Scenario), alpha = 0.4)

sub %>% 
  ggplot() +
  geom_density(aes(x = tp_pri_drug, fill = Scenario), alpha = 0.4)

sub %>% 
  ggplot() +
  geom_point(aes(x = txi_pri, y = txi_eng, colour = Scenario), alpha = 0.1)


g_ppm <- sub %>% 
  ggplot() +
  geom_density(aes(x = ppm, fill = Scenario), alpha = 0.4) +
  scale_x_continuous("PPM (treatment initiation in engaged private / all private), %", labels = scales::percent)

g_under <- sub %>% 
  ggplot() +
  geom_density(aes(x = p_under, fill = Scenario), alpha = 0.4) +
  scale_x_continuous("Under reporting (unengaged private / overall case detection), %", labels = scales::percent)


ggsave(g_ppm, filename = here::here("docs", "figs", "g_ppm.png"), width = 5, height = 3)
ggsave(g_under, filename = here::here("docs", "figs", "g_under.png"), width = 5, height = 3)




sub %>% 
  ggplot() +
  geom_density(aes(x = dur_pri, fill = Scenario), alpha = 0.4)


sub %>% 
  ggplot() +
  geom_density(aes(x = p_pub, fill = Scenario), alpha = 0.4)




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
  filter(Scenario %in% c("tx_01", "dx_11", "tx_11")) %>% 
  select(Scenario, ppv_pri, dur_pri, p_under, tp_pri_drug_time, tp_pri_txi) %>% 
  pivot_longer(-Scenario, names_to = "Index") %>% 
  ggplot() +
  geom_density(aes(x = value)) +
  facet_wrap(Index~Scenario, scales = "free_x", ncol=3)




post %>% 
  filter(Scenario %in% c("tx_01", "dx_11", "tx_11")) %>% 
  select(Scenario, tp_pri_txi, tp_pri_drug, tp_pri_drug_time) %>% 
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
  filter(Scenario %in% c("tx_01", "tx_11")) %>% 
  select(Scenario, tp_pri_txi) %>% 
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
  select(Scenario, ppv_pri, dur_pri) %>% 
  ggplot() +
  geom_point(aes(x = ppv_pri, y = dur_pri, colour = Scenario)) +
  facet_grid(.~Scenario)




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




