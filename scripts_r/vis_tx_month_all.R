library(tidyverse)
library(tidybayes)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))


# Drug sale data

targets <- read_csv(here::here("data", "cascade", "d_cascade_2019.csv")) %>% 
  select(Region, State, PrevUt, PrevTx, CNR_pub, CNR_eng) %>% 
  mutate(PrevUt = PrevUt, PrevTx = PrevTx) %>% 
  pivot_longer(c(PrevUt, PrevTx, CNR_pub, CNR_eng), names_to = "Index") %>% 
  rename(Location = State)


targets


drug <- read_csv(here::here("data", "drug_sale.csv")) %>% 
  extract(Tx_month_pri, c("TxMo_M", "TxMo_L", "TxMo_U"), "(\\S+) \\((\\S+), (\\S+)\\)", convert = T) %>% 
  extract(Ratio_private_RNTCP, c("Ratio_M", "Ratio_L", "Ratio_U"), "(\\S+) \\((\\S+), (\\S+)\\)", convert = T) %>% 
  mutate(across(starts_with("TxMo"), function(x) x * 1e-5))


targets


sims <- read_csv(here::here("out", "tx_time1.csv"))
sims




sims %>% 
  ggplot() +
  geom_point(aes(x = ppv_pri, y = dur_tx_pri)) +
  facet_grid(Target~Location)


g_margin <- sims %>% 
  mutate(
    txmo_pri = 12 * (prev_pos_pri),
    txmo_engpri = 12 * (prev_pos_pri + prev_pos_eng)
  ) %>% 
  inner_join(drug %>% rename(Location = State)) %>% 
  select(Location, Target, ppv_pri, txmo_pri, txmo_engpri, dur_tx_pri, TxMo_M, TxMo_L, TxMo_U) %>% 
  filter(txmo_pri <　TxMo_U & txmo_engpri > TxMo_L) %>% 
  ggplot() +
  geom_point(aes(x = ppv_pri, y = dur_tx_pri)) +
  facet_grid(Location~Target)



sims %>% 
  mutate(
    txmo_pri = 12 * (prev_pos_pri),
    txmo_engpri = 12 * (prev_pos_pri + prev_pos_eng)
  ) %>% 
  inner_join(drug %>% rename(Location = State)) %>% 
  select(Location, Target, ppv_pri, txmo_pri, txmo_engpri, dur_tx_pri, TxMo_M, TxMo_L, TxMo_U) %>% 
  filter(txmo_pri <　TxMo_U & txmo_engpri > TxMo_L) %>% 
  ggplot() +
  geom_point(aes(x = ppv_pri, y = dur_tx_pri)) +
  facet_grid(Location~Target)



ggsave(g_margin, filename = here::here("docs", "g_margin.png"), width = 6, height = 35)



# sims %>% 
#   group_by(Target, Location) %>% 
#   summarise(across(starts_with("ratio"), list(
#     M = median,
#     L = function(x) quantile(x, 0.25),
#     U = function(x) quantile(x, 0.75)
#   ))) %>% 
#   inner_join(drug %>% rename(Location = State)) %>% 
#   ggplot() +
#   geom_pointrange(aes(x = ratio_pe_p_M, y = Ratio_M, ymin = Ratio_L, ymax = Ratio_U, colour = Location)) +
#   geom_pointrange(aes(x = ratio_pe_p_M, xmin = ratio_pe_p_L, xmax = ratio_pe_p_U, y = Ratio_M, colour = Location)) +
#   geom_pointrange(aes(x = ratio_p_ep_M, y = Ratio_M, ymin = Ratio_L, ymax = Ratio_U, colour = Location)) +
#   geom_pointrange(aes(x = ratio_p_ep_M, xmin = ratio_p_ep_L, xmax = ratio_p_ep_U, y = Ratio_M, colour = Location)) +
#   geom_abline(slope = 1) + 
#   facet_grid(.~Target)




# sims %>% 
#   mutate(
#     txmo_pri = 12 * (prev_pos_pri),
#     txmo_engpri = 12 * (prev_pos_pri + prev_pos_eng)
#   ) %>% 
#   group_by(Target, Location) %>% 
#   summarise(across(starts_with("txmo_"), list(
#     M = median,
#     L = function(x) quantile(x, 0.25),
#     U = function(x) quantile(x, 0.75)
#   ))) %>% 
#   inner_join(drug %>% rename(Location = State)) %>% 
#   ggplot() +
#   geom_pointrange(aes(x = txmo_pri_M, y = TxMo_M, ymin = TxMo_L, ymax = TxMo_U, colour = Location)) +
#   geom_pointrange(aes(x = txmo_pri_M, xmin = txmo_pri_L, xmax = txmo_pri_U, y = TxMo_M, colour = Location)) +
#   geom_pointrange(aes(x = txmo_engpri_M, y = TxMo_M, ymin = TxMo_L, ymax = TxMo_U, colour = Location)) +
#   geom_pointrange(aes(x = txmo_engpri_M, xmin = txmo_engpri_L, xmax = txmo_engpri_U, y = TxMo_M, colour = Location)) +
#   geom_abline(slope = 1) + 
#   expand_limits(x = c(0, 0.1), y = c(0, 0.1)) +
#   facet_grid(.~Target)



g_txmo <- sims %>% 
  mutate(
    txmo_pri = 12 * (prev_pos_pri),
    txmo_engpri = 12 * (prev_pos_pri + prev_pos_eng)
  ) %>% 
  inner_join(drug %>% rename(Location = State)) %>% 
  filter(Target == "CNR+Prev") %>% 
  ggplot() + 
  geom_rect(aes(xmin = TxMo_L, xmax = TxMo_U, ymin = 0, ymax = 1), alpha = 0.01) +
  geom_vline(aes(xintercept = TxMo_M)) +
  stat_halfeye(aes(x = txmo_engpri, fill = "Engaged + Unengaged, Private "), alpha = 0.4) +
  stat_halfeye(aes(x = txmo_pri, fill = "Unengaged Private"), alpha = 0.4) +
  scale_x_continuous("Drug months, per thousand person year", 
                     labels = scales::number_format(scale = 1e5)) +
  scale_y_continuous("Density") + 
  scale_fill_discrete("") +
  facet_wrap(Location~Target, scales = "free", nrow = 5) +
  theme(legend.position = "bottom")
  

g_txmo


sims %>% 
  mutate(
    txmo_pri = 12 * (prev_pos_pri),
    txmo_engpri = 12 * (prev_pos_pri + prev_pos_eng)
  ) %>% 
  inner_join(drug %>% rename(Location = State)) %>% 
  filter(Target == "CNR+Prev") %>% 
  select(Location, Target, txmo_pri, txmo_engpri, prev_pos) %>% 
  pivot_longer(-c(Location, Target), names_to = "Index") %>% 
  group_by(Location, Target, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  )



g_cali <- sims %>% 
  select(Location, Target, PrevUt = prev_ut, PrevTx = prev_pos, CNR_pub = cnr_pub, CNR_eng = cnr_eng) %>% 
  pivot_longer(-c(Location, Target), names_to = "Index") %>% 
  group_by(Location, Target, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.025),
    U = quantile(value, 0.975)
  ) %>% 
  ggplot() + 
  geom_pointrange(aes(x = Index, y = M, ymin = L, ymax = U)) +
  geom_point(data = targets, aes(x = Index, y = value, colour = Index), size = 5, alpha = 0.2) + 
  scale_y_continuous("per 100k", labels = scales::number_format(scale = 1e5)) +
  facet_grid(Location~Target) +
  guides(colour = guide_none())

g_cali





ggsave(g_txmo, filename = here::here("docs", "g_txmo.png"), width = 9, height = 12)

ggsave(g_cali, filename = here::here("docs", "g_cali.png"), width = 6, height = 35)







