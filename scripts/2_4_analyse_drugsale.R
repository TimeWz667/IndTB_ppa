library(tidyverse)



loc <- "India"

names.pars <- c("inc0", "p_ppm", "ppv_pri", "dur_pri", 
                "p_pri_on_pub",
                "r_onset", "r_csi", "r_recsi", "r_acf", 
                "p_dx1_pub", "p_dx0_pub", "p_dx1_pri", "p_dx1_pri")



tab0 <- local({
  load(here::here("out", "sub_csd", "post_csd_d0_" + glue::as_glue(loc) + ".rdata"))
  as_tibble(rstan::extract(post, pars=names.pars)) %>% mutate(Location = loc, Scenario = "No drug sale")
})



tab1 <- local({
  load(here::here("out", "sub_csd", "post_csd_d1_" + glue::as_glue(loc) + ".rdata"))
  as_tibble(rstan::extract(post, pars=names.pars)) %>% mutate(Location = loc, Scenario = "With drug sale")
})


tab2 <- local({
  load(here::here("out", "sub_csd", "post_csd_d2_" + glue::as_glue(loc) + ".rdata"))
  as_tibble(rstan::extract(post, pars=names.pars)) %>% mutate(Location = loc, Scenario = "No Private")
})


tab <- bind_rows(tab0, tab1, tab2)


tab %>% 
  ggplot() +
  geom_density(aes(x = inc0, fill = Scenario), alpha = 0.5)




d <- as_tibble(dat)


d %>% 
  mutate(Prev = (Asym + Sym + CS) / N, CNR = (NotiPub + NotiEng) / Pop) %>% 
  select(Prev, CNR)


tab %>% 
  ggplot() +
  geom_point(aes(x = p_pri_on_pub, y = inc0, colour = Scenario))


