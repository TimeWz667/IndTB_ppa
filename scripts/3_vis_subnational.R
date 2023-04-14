library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


locs <- local({
  locs <- dir(here::here("out", "tx_11_subnational"))
  locs <- locs[startsWith(locs, "post_")]
  locs <- gsub("post_", "", locs)
  locs <- gsub(".csv", "", locs)
  locs
  
})


post <- bind_rows(lapply(locs, function(loc) {
  read_csv(here::here("out", "tx_11_subnational", "post_" + glue::as_glue(loc) +".csv")) %>% 
    select(ppm, dur_pri, ppv_pri, p_pri_on_pub, p_under,
           tp_pri_drug, tp_pri_drug_time, tp_pri_txi) %>% 
    mutate(State = loc)
}))


stats <- post %>% 
  pivot_longer(-State, names_to = "Index") %>% 
  group_by(State, Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.25),
    U = quantile(value, 0.75)
  ) %>% 
  ungroup()


post %>% 
  mutate(
    State = reorder(State, ppm)
  ) %>% 
  ggplot() +
  stat_halfeye(aes(x = ppm, y = State))





post %>% 
  mutate(
    State = reorder(State, dur_pri)
  ) %>% 
  ggplot() +
  stat_halfeye(aes(x = dur_pri, y = State))



post %>% 
  ggplot() +
  geom_point(aes(x = dur_pri, y = ppv_pri)) +
  facet_wrap(.~State)



stats %>% 
  filter(Index == "ppm") %>% 
  mutate(
    State = reorder(State, M)
  ) %>% 
  ggplot() +
  geom_pointinterval(aes(x = M, xmin = L, xmax = U, y = State))


g_dur_pri <- stats %>% 
  filter(Index == "dur_pri") %>% 
  mutate(
    State = reorder(State, M)
  ) %>% 
  ggplot() +
  geom_pointinterval(aes(x = M, xmin = L, xmax = U, y = State)) + 
  scale_x_continuous("Treatment duration, on private drug, month", labels = scales::number_format(scale = 12)) +
  expand_limits(x = 0)


g_dur_pri


g_ppv_pri <- stats %>% 
  filter(Index == "ppv_pri") %>% 
  mutate(
    State = reorder(State, M)
  ) %>% 
  ggplot() +
  geom_pointinterval(aes(x = M, xmin = L, xmax = U, y = State)) + 
  scale_x_continuous("Positive Predictive Value, in the private sector, %", labels = scales::percent_format()) + 
  expand_limits(x = c(0, 1))

g_ppv_pri



ggsave(g_ppv_pri, filename = here::here("docs", "figs", "g_post_ppv_sn.png"), width = 6, height = 7)

ggsave(g_dur_pri, filename = here::here("docs", "figs", "g_post_dur_sn.png"), width = 6, height = 7)




g_drug <- stats %>% 
  filter(Index == "tp_pri_drug") %>% 
  mutate(
    State = reorder(State, M)
  ) %>% 
  ggplot() +
  geom_pointinterval(aes(x = M, xmin = L, xmax = U, y = State)) + 
  scale_x_continuous("Txi with private drugs, millions", labels=scales::number_format(scale = 1e-6)) + 
  expand_limits(x = c(0, 1))

g_drug


g_drugT <- stats %>% 
  filter(Index == "tp_pri_drug_time") %>% 
  mutate(
    State = reorder(State, M)
  ) %>% 
  ggplot() +
  geom_pointinterval(aes(x = M, xmin = L, xmax = U, y = State)) + 
  scale_x_continuous("Caseloads with private drugs, millions", labels=scales::number_format(scale = 1e-6)) + 
  expand_limits(x = c(0, 1))

g_drugT



g_txi<- stats %>% 
  filter(Index == "tp_pri_txi") %>% 
  mutate(
    State = reorder(State, M)
  ) %>% 
  ggplot() +
  geom_pointinterval(aes(x = M, xmin = L, xmax = U, y = State)) + 
  scale_x_continuous("Txi  in the unengaged private, millions", labels=scales::number_format(scale = 1e-6)) + 
  expand_limits(x = c(0, 1))

g_txi




ggsave(g_drug, filename = here::here("docs", "figs", "g_post_drug_sn.png"), width = 6, height = 7)
ggsave(g_drugT, filename = here::here("docs", "figs", "g_post_drugt_sn.png"), width = 6, height = 7)
ggsave(g_txi, filename = here::here("docs", "figs", "g_post_txi_sn.png"), width = 6, height = 7)



stats %>% 
  filter(Index %in% c("tp_pri_txi", "tp_pri_drug_time", "tp_pri_drug")) %>% 
  mutate(
    State = reorder(State, M)
  ) %>% 
  ggplot() +
  geom_pointinterval(aes(x = M, xmin = L, xmax = U, y = State, colour = Index), position = position_dodge(-0.3)) + 
  scale_x_log10("Caseloads, millions", labels=scales::number_format(scale = 1e-6)) + 
  expand_limits(x = c(0, 1))





