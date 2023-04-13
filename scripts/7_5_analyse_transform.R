library(tidyverse)



loc <- "India"
loc <- "Andhra_Pradesh"

sims <- read_csv(here::here("out", "sub_cas", "sims_baseline_" + glue::as_glue(loc) + ".csv")) %>% 
  mutate(
    Model = ifelse(Model == "Juan", "Eq0", Model),
    Model = factor(Model, c("Eq", "Eq0", "AllRe", "Free"))
  )


sims %>%
  filter(t == 2020) %>% 
  ggplot() +
  geom_histogram(aes(x = det_pub, fill = Model), alpha = 0.6, position = position_identity())


sims %>%
  filter(t == 2020) %>% 
  ggplot() +
  geom_histogram(aes(x = det_eng, fill = Model), alpha = 0.6, position = position_identity())



sims %>%
  filter(t == 2020) %>% 
  ggplot() +
  geom_histogram(aes(x = det_pri, fill = Model), alpha = 0.6, position = position_identity())


sims %>%
  filter(t == 2020) %>% 
  ggplot() +
  geom_histogram(aes(x = del_pat, fill = Model), alpha = 0.6, position = position_identity())

sims %>%
  filter(t == 2020) %>% 
  ggplot() +
  geom_histogram(aes(x = del_sys, fill = Model), alpha = 0.6, position = position_identity())



sims %>%
  ggplot() +
  geom_histogram(aes(x = p_dx0_pub, fill = Model), alpha = 0.6, position = position_identity())


sims %>%
  ggplot() +
  geom_histogram(aes(x = p_dx1_pub, fill = Model), alpha = 0.6, position = position_identity())


sims %>%
  ggplot() +
  geom_histogram(aes(x = r_recsi, fill = Model), alpha = 0.6, position = position_identity())


sims %>%
  ggplot() +
  geom_histogram(aes(x = r_csi, fill = Model), alpha = 0.6, position = position_identity())




sims <- read_csv(here::here("out", "sub_cas", "sims_" + glue::as_glue(loc) + ".csv")) %>% 
  mutate(
    Model = ifelse(Baseline == "Zero0", "AllRe", Baseline)
  )


sims %>%
  ggplot() +
  geom_histogram(aes(x = p_dx1_pub, fill = Model), alpha = 0.6, position = position_identity()) +
  facet_grid(Scenario~.)


sims %>%
  ggplot() +
  geom_histogram(aes(x = p_dx0_pub, fill = Model), alpha = 0.6, position = position_identity()) +
  facet_grid(Scenario~.)


sims %>%
  ggplot() +
  geom_histogram(aes(x = Delay0, fill = Model), alpha = 0.6, position = position_identity()) +
  facet_grid(Scenario~.)


sims %>%
  ggplot() +
  geom_histogram(aes(x = CDR0, fill = Model), alpha = 0.6, position = position_identity()) +
  facet_grid(Scenario~.)


sims %>%
  mutate(diffDelay = Delay0 - Delay1) %>% 
  ggplot() +
  geom_histogram(aes(x = diffDelay * 12, fill = Model), alpha = 0.6, position = position_identity()) +
  facet_grid(Scenario~.)


sims %>%
  mutate(diffCDR = CDR0 - CDR1) %>% 
  ggplot() +
  geom_histogram(aes(x = diffCDR * 100, fill = Model), alpha = 0.6, position = position_identity()) +
  facet_grid(Scenario~.)










