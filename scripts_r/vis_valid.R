library(tidyverse)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))



sims <- read_csv(here::here("out", "rate_valid", "valid.csv"))



sel_ind <- c(
  IncR = "Incidence",
  Prev = "Prevalence",
  CNR_Pub = "Case notif.\npub.",
  CNR_Eng = "Case notif.\nengaged pri."
)


stats <- sims %>% 
  pivot_longer(-Key, names_to = "Index") %>% 
  group_by(Index) %>% 
  summarise(
    M = median(value),
    L = quantile(value, 0.1),
    U = quantile(value, 0.9)
  ) %>% 
  filter(Index %in% names(sel_ind))


dat <- read_csv(here::here("data",  "cascade", "d_cascade_2019.csv")) %>% 
  filter(State == "India") %>% 
  mutate(Prev = PrevUt * 1e-5, IncR = 193 * 1e-5) %>% 
  rename(CNR_Pub = CNR_pub, CNR_Eng = CNR_eng) %>% 
  select(Prev, CNR_Pub, CNR_Eng, IncR) %>% 
  pivot_longer(everything(), names_to = "Index")


g_valid <- stats %>% 
  left_join(dat) %>% 
  mutate(
    Index = factor(Index, names(sel_ind))
  ) %>% 
  ggplot() +
  geom_point(aes(x = Index, y = value, colour = "Data"), size = 4, alpha = 0.8) +
  geom_pointrange(aes(x = Index, y = M, ymax = U, ymin = L)) +
  scale_y_continuous("per 100 000", labels = scales::number_format(scale = 1e5)) +
  scale_x_discrete("", labels = sel_ind) +
  scale_colour_discrete("") +
  expand_limits(y = 0) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))


g_valid

bd <- sims %>% 
  crossing(DurPri = seq(0, 24, 0.25), PPV_Pri = seq(0.01, 1, 0.01)) %>% 
  group_by(DurPri, PPV_Pri) %>%
  summarise(
    dt = mean((Det_Pri / PPV_Pri * PrTxi_pri * DurPri + Dur_Eng * CNR_Eng) * 1310 > 17)
  )


# 


g_cross <- bd %>% 
  ggplot() + 
  geom_raster(aes(x = PPV_Pri, y = DurPri, fill = dt)) +
  scale_x_continuous("PPV in unengaged private services, %", labels = scales::percent) +
  scale_y_continuous("Treatment duration in unengaged private services, months", breaks = seq(0, 24, 6)) +
  scale_fill_distiller("Pr(patient-month estimate > drug-sale data)", labels = scales::percent) +
  theme(legend.position = "bottom") +
  labs(subtitle = "Patient-month from drug-sale data: 17 million")



g_cross 

ggsave(g_valid, filename = here::here("out", "rate_valid", "g_valid.png"), width = 5, height = 4)

ggsave(g_cross, filename = here::here("out", "rate_valid", "g_drug_time.png"), width = 5, height = 4)

