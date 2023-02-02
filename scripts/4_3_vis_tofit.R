library(tidyverse)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))


load(here::here("data", "dat_cas.rdata"))
tofit <- read_csv(here::here("docs", "tabs", "to_fit.csv"))



sim_prev <- tofit %>% 
  filter(name %in% c("pr_a", "pr_s", "pr_c")) %>% 
  select(State, name, Type = Scenario, M = `50%`, L = `25%`, U = `75%`)


d_prev <- dat_tbps %>% 
  mutate(
    N = N_Asym + N_NotCS + N_NotDet,
    pr_a = N_Asym / N,
    pr_s = N_NotCS / N,
    pr_c = N_NotDet / N
  ) %>% 
  select(State, N, starts_with("pr_")) %>% 
  pivot_longer(starts_with("pr")) %>% 
  mutate(
    L = qbinom(0.25, size=N, prob=value) / N,
    U = qbinom(0.75, size=N, prob=value) / N,
    Type = "Data"
  ) %>% 
  select(State, name, Type,  M = value, L, U)


g_prev <- bind_rows(sim_prev, d_prev) %>% 
  mutate(
    Type = factor(Type, c("Data", "s0", "s1", "s2", "s3")),
    name = factor(name, c("pr_a", "pr_s", "pr_c"))
  ) %>% 
  ggplot() +
  geom_pointrange(aes(x = name, y = M, ymin = L, ymax = U, colour = Type), position = position_dodge(width = 0.6)) +
  facet_wrap(.~State, ncol = 4) +
  scale_y_continuous("Proportion, %", labels = scales::percent) + 
  scale_x_discrete("Stage", labels = c(pr_a="Subclinical", pr_s="Pre-care seeking", pr_c="Pre-case detection")) + 
  scale_colour_discrete("", labels = c(Data="Target data", s0="~.", s1="~(1)", s2="~(2)", s3="~(1)+(2)")) +
  expand_limits(y = c(0, 1)) +
  labs(caption = "(1) private ATT usage from drug sale data\n(2) prevalent TB on private ATT") +
  theme(axis.text.x = element_text(angle = -35, hjust=0))

g_prev


ggsave(g_prev, filename = here::here("docs", "figs", "g_fit_prev.png"), width = 10, height = 12) 



