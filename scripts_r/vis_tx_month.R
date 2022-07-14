library(tidyverse)
library(jsonlite)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))

# Drug sale data
sel_loc <- c("India", "Kerala", "Uttar Pradesh", "Gujarat")


dat <- read_csv(here::here("data", "cascade", "d_cascade_2019.csv")) %>% 
  filter(State %in% sel_loc) %>% 
  select(Location = State, N = Pop)

drug_sale <- tibble(
  Location = sel_loc,
  PT_M = c(1350.0, 473.0, 2400.0, 1320.0)
)



ds <- dir(here::here("out", "cascade_2019"))
ds <- ds[startsWith(ds, "Cascade_")]


cas <- bind_rows(lapply(ds, function(d) {
  print(d)
  css <- read_json(here::here("out", "cascade_2019", d))
  
  tibble(
    Location = gsub(".json", "", gsub("Cascade_", "", d)), 
    DetR_Pri = sapply(css, function(cs) cs$DetR[[3]]),
    Prop_Pri = sapply(css, function(cs) cs$DetR[[3]] / sum(unlist(cs$DetR))),
    TxI_Pri = sapply(css, function(cs) cs$P_TxI[[3]]),
    PrevTx_Eng = sapply(css, function(cs) cs$PrevTx[[2]]),
    PrevTx_Pri = sapply(css, function(cs) cs$PrevTx[[3]]),
    Inc = sapply(css, function(cs) cs$PrevUt * (cs$R_Die_Asym + cs$R_Onset + cs$R_SelfCure))
  )
}))



g_dense <- cas %>% 
  filter(Location %in% sel_loc) %>% 
  group_by(Location) %>% 
  mutate(Location = factor(Location, sel_loc)) %>% 
  ggplot() + 
  geom_density(aes(x = (PrevTx_Pri + PrevTx_Eng) * 12), fill = "grey8", alpha = 0.3) +
  geom_vline(data = drug_sale %>% mutate(PT_M = PT_M * 1e-5, Location = factor(Location, sel_loc)), 
             aes(xintercept = PT_M)) + 
  scale_y_continuous("Density", labels = NULL) +
  scale_x_continuous("Patient-month, private sector, per 10 000 population", labels = scales::number_format(scale = 1e5)) +
  expand_limits(x = 0) +
  facet_grid(.~Location) +
  theme(axis.text.x = element_text(angle = - 30))


cas %>% 
  filter(Location %in% sel_loc) %>% 
  group_by(Location) %>% 
  mutate(Location = factor(Location, sel_loc)) %>% 
  ggplot() + 
  geom_density(aes(x = Inc), fill = "grey8", alpha = 0.3) +
  scale_y_continuous("Density", labels = NULL) +
  scale_x_continuous("Incidence, per 10 000 population", labels = scales::number_format(scale = 1e5)) +
  expand_limits(x = 0) +
  facet_grid(.~Location) +
  theme(axis.text.x = element_text(angle = - 30))


g_inc <-cas %>% 
  filter(Location %in% sel_loc) %>%
  group_by(Location) %>% 
  summarise(
    M = median(Inc), 
    L = quantile(Inc, 0.25),
    U = quantile(Inc, 0.75)
  ) %>% 
  mutate(Location = factor(Location, sel_loc)) %>% 
  ggplot() + 
  geom_pointrange(aes(x = Location, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Incidence, per 10 000 population", labels = scales::number_format(scale = 1e5)) +
  expand_limits(y = 0)



g_prop <-cas %>% 
  filter(Location %in% sel_loc) %>%
  group_by(Location) %>% 
  summarise(
    M = median(Prop_Pri), 
    L = quantile(Prop_Pri, 0.25),
    U = quantile(Prop_Pri, 0.75)
  ) %>% 
  mutate(Location = factor(Location, sel_loc)) %>% 
  ggplot() + 
  geom_pointrange(aes(x = Location, y = M, ymin = L, ymax = U)) +
  scale_y_continuous("Proportion private treated TB, %", labels = scales::percent) +
  expand_limits(y = c(0, 1))
  
g_prop



bd <- cas %>% 
  filter(Location %in% sel_loc) %>% 
  crossing(PPV_Pri = seq(0.05, 1, 0.05), DurTx_Pri = seq(2, 12, 0.25)) %>% 
  left_join(drug_sale) %>%
  mutate(
    PT_M = PT_M * 1e-5,
    DrugTime = DetR_Pri / PPV_Pri * TxI_Pri * DurTx_Pri
  ) %>%  
  group_by(Location, DurTx_Pri, PPV_Pri) %>%
  summarise(
    dt = mean(DrugTime > PT_M)
  ) 

 


g_cross <- bd %>% 
  mutate(Location = factor(Location, sel_loc)) %>% 
  ggplot() + 
  geom_raster(aes(x = PPV_Pri, y = DurTx_Pri, fill = dt)) +
  scale_x_continuous("PPV in unengaged private services, %", labels = scales::percent) +
  scale_y_continuous("Treatment duration, unengaged private, months", breaks = seq(0, 24, 6)) +
  scale_fill_distiller("Pr(patient-month estimate from TBPS > drug-sale data)", labels = scales::percent) +
  theme(legend.position = "bottom") + 
  facet_wrap(.~Location, nrow = 1)

g_cross


ggsave(g_prop, filename = here::here("docs", "prop_pri.pdf"), width = 5, height = 3)


ggsave(g_cross, filename = here::here("docs", "drug_month.pdf"), width = 11, height = 4.5)


ggsave(g_dense, filename = here::here("docs", "drug_month_dense.pdf"), width = 11, height = 4)


