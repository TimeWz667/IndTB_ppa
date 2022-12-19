library(tidyverse)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))


summariser <- list(
  M = median,
  L = function(x) quantile(x, 0.05),
  U = function(x) quantile(x, 0.95)
)


ds <- dir("out/sim/ss/")


sims <- bind_rows(lapply(ds, function(d) {
  read_csv(here::here("out", "sim", "ss", d))
})) %>% 
  select(-Key, -N) %>% 
  group_by(Time, Source) %>%
  summarise(across(everything(), summariser)) %>% 
  ungroup() 


sims %>% 
  extract(Source, c("Type", "Location"), "pars_(s1|s2|s3)_(\\S+)") %>% 
  pivot_longer(-c(Time, Type, Location)) %>% 
  extract(name, c("Index", "Stats"), "(\\S+)_(M|L|U)") %>% 
  pivot_wider(names_from = Stats, values_from = value)


sims %>% write_csv(here::here("docs", "tabs", "sim_tofit.csv"))


sims_pri <- bind_rows(lapply(ds, function(d) {
  read_csv(here::here("out", "sim", "ss", d))
})) %>% 
  filter(Time == 2020) %>% 
  mutate(CarePri = TxEng + TxPri) %>% 
  extract(Source, c("Type", "Location"), "pars_s(1-3)_(\\S+)") %>% 
  select(Key, Type, Location, CarePri, TxPri, TxEng)



targets <- read_csv(here::here("data", "targets.csv")) %>% rename(Location = State)


g_fit_prev <- sims %>% 
  filter(startsWith(Index, "Prev")) %>% 
  filter(Time == 2020) %>% 
  filter(Type == "onds") %>% 
  ggplot() + 
  geom_pointrange(aes(x = Index, y = M, ymin = L, ymax = U)) +
  geom_point(data = targets %>% filter(Index %in% c("PrevA", "PrevS", "PrevC")),
             aes(x = Index, y = Observed), colour = 'red', alpha = 0.3, size = 4) + 
  scale_y_continuous("Prevalence, per 100 000", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0) +
  facet_wrap(.~ Location, ncol = 6)

g_fit_prev
ggsave(g_fit_prev, filename = here::here("docs", "figs", "g_fit_onds_prev.png"), width = 12, height = 10)



g_fit_prev <- sims %>% 
  filter(startsWith(Index, "Prev")) %>% 
  filter(Time == 2020) %>% 
  filter(Type == "nods") %>% 
  ggplot() + 
  geom_pointrange(aes(x = Index, y = M, ymin = L, ymax = U)) +
  geom_point(data = targets %>% filter(Index %in% c("PrevA", "PrevS", "PrevC")),
             aes(x = Index, y = Observed), colour = 'red', alpha = 0.3, size = 4) + 
  scale_y_continuous("Prevalence, per 100 000", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0) +
  facet_wrap(.~ Location, ncol = 6)

g_fit_prev
ggsave(g_fit_prev, filename = here::here("docs", "figs", "g_fit_nods_prev.png"), width = 12, height = 10)



for (type in c("onds", "nods")) {
  for (cnr in c("CNR_Pub", "CNR_Eng", "CNR_Acf")) {
    g_fit_cnr <- sims %>% 
      filter(Index == cnr) %>%
      filter(Type == type) %>% 
      ggplot() +
      geom_ribbon(aes(x = Time, ymin = L, ymax = U), alpha = 0.5) + 
      geom_line(aes(x = Time, y = M)) +
      geom_point(data = targets %>% filter(Index == cnr), aes(x = Year, y = Observed), colour = "red") +
      scale_y_continuous("Case notification, per 100 000", labels = scales::number_format(scale = 1e5)) +
      scale_x_continuous("Year", breaks = c(2015, 2020, 2025)) + 
      expand_limits(y = 0) +
      facet_wrap(.~ Location, ncol = 6) +
      theme(axis.text.x = element_text(angle = 40))  
    
    
    ggsave(g_fit_cnr, filename = here::here("docs", "figs", glue::as_glue("g_fit_") + type + "_" + cnr + ".png"), width = 12, height = 10)
  }
}



g_inc <- sims %>% 
  filter(Index == "IncR") %>% 
  filter(Time == 2022) %>% 
  mutate(
    Location = factor(Location),
    Location = reorder(Location, M)
  ) %>% 
  ggplot() + 
  geom_pointrange(aes(x = M, xmin = L, xmax = U, y = Location, colour = Type), position = position_dodge(-0.4)) +
  scale_x_continuous("Incidence rate, 2020, per 100 000", labels = scales::number_format(scale = 1e5)) +
  scale_colour_discrete("Drug sale data", labels = c(nods="No", onds = "Yes")) +
  expand_limits(x = 0) +
  theme(legend.position = c(1, 0), legend.justification = c(1.05, -0.05))

ggsave(g_inc, filename = here::here("docs", "figs", "g_inc.png"), width = 6, height = 5)


tar_drug <- targets %>% filter(startsWith(Index, "OnPriDrug")) %>% pivot_wider(names_from = "Index", values_from = "Observed")
tar_care <- targets %>% filter(Index == "PrevTxPri")


sims_pri %>% 
  ggplot() +
  geom_point(aes(x = CarePri, y = TxPri, colour = Type), alpha = 0.2) +
  facet_wrap(.~ Location, ncol = 6)


sims %>% 
  filter(Index %in% c("CarePri", "TxPri")) %>% 
  filter(Time == 2020) %>% 
  pivot_wider(names_from = "Index", values_from = c("M", "L", "U")) %>% 
  ggplot() +
  geom_pointrange(aes(x = M_TxPri, y = M_CarePri, ymin = L_CarePri, ymax = U_CarePri, colour = Type)) +
  geom_linerange(aes(x = M_TxPri, xmin = L_TxPri, xmax = U_TxPri, y = M_CarePri, colour = Type)) +
  geom_vline(data = tar_drug, aes(xintercept = OnPriDrugM)) +
  geom_vline(data = tar_drug, aes(xintercept = OnPriDrugL), linetype = 2) +
  geom_vline(data = tar_drug, aes(xintercept = OnPriDrugU), linetype = 2) +
  geom_hline(data = tar_care, aes(yintercept = Observed)) +
  facet_wrap(.~ Location, ncol = 6)
  







