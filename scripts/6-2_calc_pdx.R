library(tidyverse)
library(jsonlite)



p_dx1_pub = 0.7
rat_dx01 = 0.4


load(here::here("data", "dat_cas.rdata"))


loc <- glue::as_glue("India")


sel_cs <- read_csv(here::here("out",  "sub_cs", "post_cs_s1_" + loc + ".csv")) %>% 
  mutate(Key = sample(1:n(), n(), rep=F)) %>% 
  arrange(Key) %>% 
  select(Key, adr, prv0, 
         r_death_a, r_death_s, r_sc,
         r_sym, r_csi = r_aware,
         r_det, r_det_all, 
         det_pub, det_eng, det_pri, ppv_pub, ppv_eng, ppv_pri, p_pri_on_pub,
         r_acf, sens_acf, spec_acf,
         dur_pub, dur_pri,
         p_txi_pri=txi_pri)


sel_cs %>% 
  summarise(
    pr_det_pub = mean(det_pub)
  )



d_tbps <- read_csv(here::here("data", "TBPS", "TBPS_ASC_Region.csv"))

d_tbps %>% 
  mutate(pr_pub = N_Pub / (N_Pub + N_Pri), pr_pub = N_Pub / (N_Pub + N_Pri)) %>% 
  select(State, pr_pub, Pr_Pub_CSI0, Pr_Pub_CSI1, Pr_Pub_CSI)




x_pub <- 0.483
x_pri <- 1 - x_pub
txi_pub <- 0.7
pdx_pub <- 0.66
pdx_pri <- x_pub * pdx_pub * (1 - txi_pub) / txi_pub / x_pri 

pdx_pri

pdx = x_pub * pdx_pub + x_pri * pdx_pri
pdx


(1 / 2.7 - x_pub * pdx_pub) / x_pri
(1 / 1.9 - x_pub * pdx_pub) / x_pri
(1 / 12.3 - x_pub * pdx_pub) / x_pri








r_csi <- 1 / 0.25

3 * pdx / pdx_pub









