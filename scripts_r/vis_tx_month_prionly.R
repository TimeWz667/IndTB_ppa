library(tidyverse)


dat <- read_csv(here::here("data", "cascade", "d_cascade_2019.csv")) %>% 
  select(Region, State, PrevUt, PrevTx, PrevTxPri, starts_with("DrugTime_"), starts_with("CNR"))


dat <- as.list(dat[dat$State == "India", ])

dat
# 
# dat %>% 
#   ggplot() +
#   geom_point(aes(x = PrevTxPri * 12, y = TxMo_M)) +
#   geom_pointrange(aes(x = PrevTxPri * 12, y = TxMo_M, ymin = TxMo_L, ymax = TxMo_U)) +
#   geom_abline(slope = 1) +
#   expand_limits(x = 0, y = 0)





sim0 <- bind_rows(jsonlite::read_json(here::here("out", "cds0_pars", "fit_India.json")))

sim1 <- bind_rows(jsonlite::read_json(here::here("out", "cds1_pars", "fit_India.json")))

sim <- bind_rows(
  sim0 %>% mutate(Scenario = "No drug-sale data"),
  sim1 %>% mutate(Scenario = "With drug-sale data")
) 


g1 <- sim %>% 
  ggplot() + 
  geom_point(aes(x = PpvPri, y = PrDetPri, colour = Scenario)) +
  scale_y_continuous("Cases detected by unengaged private system %, true TB", label = scales::percent) +
  scale_x_continuous("PPV %, unengaged private", label = scales::percent) +
  expand_limits(x = 0:1, y = 0:1)

g1
# sim %>% 
#   ggplot() +
#   geom_density(aes(x = PrOnPriDrug, fill = Scenario), alpha = 0.3)


g2 <- sim %>% 
  ggplot() + 
  geom_density(aes(x = OnPriDrug, fill = Scenario), alpha = 0.3) + 
  geom_segment(x=dat$DrugTime_M, xend=dat$DrugTime_M, y=0, yend=Inf, linetype = 2) +
  # geom_segment(x=dsd[2], xend=dsd[2], y=0, yend=Inf, linetype = 3) +
  # geom_text(x=dat$DrugTime_M, y=6500, label = "2019 drug-sale data", hjust = -0.1) +
  scale_x_continuous("Patient months of private TB treatment per 100k",
                     label = scales::number_format(scale = 12e5))

g2


d <- tibble(
  Sector = c("Public", "Private"),
  value = c(dat$CNR_pub, dat$CNR_eng)
)


sim %>% 
  select(Scenario, Public = CnrPub, Private = CnrEng) %>% 
  pivot_longer(-Scenario, names_to = "Sector") %>% 
  ggplot() +
  geom_density(aes(x=value, fill = Scenario), alpha = 0.3) +
  geom_segment(data = d, aes(x=value, xend=value, y=0, yend=Inf)) +
  facet_grid(.~Sector)
  



sim %>% 
  # filter(OnPriDrug < dat$DrugTime_U & OnPriDrug > dat$DrugTime_L) %>% 
  ggplot() +
  geom_density(aes(x = DrugTimeEng + DrugTimePri, fill = Scenario), alpha = 0.3) +
  geom_segment(aes(x=dat$PrevTxPri, xend=dat$PrevTxPri, y=0, yend=Inf))


g3 <- sim %>% 
  mutate(Fil = PrOnPriDrug > 0.7) %>% 
  ggplot() +
  geom_point(aes(x = DrugTimeEng + DrugTimePri, y = OnPriDrug, colour = Scenario), alpha = 0.3) +
  # geom_density_2d(aes(x = DrugTimeEng + DrugTimePri, y = OnPriDrug, colour = Scenario)) +
  geom_vline(xintercept = dat$PrevTxPri) +
  geom_hline(yintercept = dat$DrugTime_M) +
  scale_x_continuous("Prevalent cases on private TB service per 100k", label=scales::number_format(scale = 1e5)) +
  scale_y_continuous("Patient months of private TB treatment per 100k", label=scales::number_format(scale = 12e5)) +
  scale_colour_discrete("Target: TBPS + CNR") +
  expand_limits(x = 0, y = 0)
  

g4 <- sim %>% 
  mutate(Fil = ifelse(PrOnPriDrug > 0.7, "> 70%", "<70%")) %>% 
  ggplot() +
  geom_point(aes(x = DrugTimeEng + DrugTimePri, y = OnPriDrug, colour = Scenario), alpha = 0.3) +
  # geom_density_2d(aes(x = DrugTimeEng + DrugTimePri, y = OnPriDrug, colour = Scenario)) +
  geom_vline(xintercept = dat$PrevTxPri) +
  geom_hline(yintercept = dat$DrugTime_M) +
  scale_x_continuous("Prevalent cases on private TB service per 100k", label=scales::number_format(scale = 1e5)) +
  scale_y_continuous("Patient months of private TB treatment per 100k", label=scales::number_format(scale = 12e5)) +
  scale_colour_discrete("Target: TBPS + CNR") +
  facet_grid(.~Fil) +
  expand_limits(x = 0, y = 0)



ggsave(g1, filename = here::here("docs", "gds1.png"), width = 8, height = 6.5)
ggsave(g2, filename = here::here("docs", "gds2.png"), width = 8, height = 6.5)
ggsave(g3, filename = here::here("docs", "gds3.png"), width = 8, height = 6.5)
ggsave(g4, filename = here::here("docs", "gds4.png"), width = 12, height = 6.5)




# g1 <- sim %>% 
#   ggplot() +
#   geom_histogram(aes(x = pridrug_pri * 12)) +
#   geom_vline(data = dat %>% filter(State == "India"), aes(xintercept = TxMo_M)) +
#   geom_vline(data = dat %>% filter(State == "India"), aes(xintercept = TxMo_L), linetype = 2) +
#   geom_vline(data = dat %>% filter(State == "India"), aes(xintercept = TxMo_U), linetype = 2) +
#   scale_x_continuous("Months with private drugs, months per 100k", label = scales::number_format(scale = 1e5)) + 
#   scale_y_continuous("Density") +
#   expand_limits(x = 0) +
#   theme(axis.text.y = element_blank())


g1 <- sim %>% 
  ggplot() +
  geom_density(aes(x = pridrug_pri * 12 + pridrug_eng /0.8 * 0.5 * 12, fill="Model projection"), alpha = 0.3) +
  geom_segment(x=dsd[1], xend=dsd[1], y=0, yend=Inf, linetype = 2) +
  # geom_segment(x=dsd[2], xend=dsd[2], y=0, yend=Inf, linetype = 3) +
  geom_text(x=dsd[1], y=700, label = "2019 drug-sale data", hjust = -0.1) +
  # geom_text(x=dsd[2], y=700, label = "2020 drug-sale data", hjust = -0.1) +
  scale_x_continuous("Patient months of private TB treatment per 100k, 2019", label = scales::number_format(scale = 1e5)) +
  scale_colour_brewer(palette = "Paired") +
  scale_y_continuous("Density") + 
  expand_limits(x = c(0, 0.008), y = 700) +
  theme(axis.text.y = element_blank(), legend.position = "none")




g1

g2 <- sim %>% 
  ggplot() +
  geom_point(aes(x = ppv_pri, y = dur_pri, colour = p_pridrug)) + 
  scale_x_continuous("PPV, unengaged private, %", label = scales::percent) +
  scale_y_continuous("Treatment duration, unengaged private, year")

g2


ggsave(g1, filename = here::here("docs", "drug time.png"), width = 5, height = 4)
ggsave(g2, filename = here::here("docs", "ppv_dur.png"), width = 5, height = 4)


