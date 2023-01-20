library(tidyverse)




locations <- dir(here::here("out", "pars"))
locations <- gsub("pars_", "", locations)
locations <- gsub(".csv", "", locations)


locations


loc <- "India"

res <- bind_rows(lapply(locations, function(loc) {
  pars <- read_csv(here::here("out", "pars", "pars_" + glue::as_glue(loc) + ".csv")) %>% 
    select(adr, starts_with("r_death_"), starts_with("ppv_"), 
           starts_with("p_"), starts_with("r_"), starts_with("dur_"), ends_with("_acf")) %>% 
    mutate(
      p_to_pub = p_txs_pub + p_txl_pub + p_txd_pub,
      p_txs_pub = p_txs_pub / p_to_pub,
      p_txl_pub = p_txl_pub / p_to_pub, 
      p_txd_pub = p_txd_pub / p_to_pub,
      p_to_eng = p_txs_eng + p_txl_eng + p_txd_eng,
      p_txs_eng = p_txs_eng / p_to_eng,
      p_txl_eng = p_txl_eng / p_to_eng, 
      p_txd_eng = p_txd_eng / p_to_eng,
      p_txs_pri = p_txs_eng,
      p_txl_pri = p_txl_eng,
      p_txd_pri = p_txd_eng
    ) %>% 
    select(-p_txi0, -p_txi1, -r_det, -p_to_pub, -p_to_eng)
  
  
  pars %>% 
    mutate(
      p_txi0_pub = p_ent_pub * p_dx0_pub * p_txi_pub,
      p_txi0_eng = p_ent_eng * p_dx0_eng * p_txi_eng,
      p_txi0_pri = p_ent_pri * p_dx0_pri * p_txi_pri,
      p_txi0 = p_txi0_pub + p_txi0_eng + p_txi0_pri,
      
      p_txi1_pub = p_ent_pub * p_dx1_pub * p_txi_pub,
      p_txi1_eng = p_ent_eng * p_dx1_eng * p_txi_eng,
      p_txi1_pri = p_ent_pri * p_dx1_pri * p_txi_pri,
      p_txi1 = p_txi1_pub + p_txi1_eng + p_txi1_pri,
      
      r_acf_tp = r_acf * sens_acf,
      
      r_drop_a = r_death_a + r_sc,
      r_prog_a = r_sym,
      r_det_a = r_acf_tp,
      r_out_a = r_drop_a + r_prog_a + r_det_a,
      Dur_A = 1 / r_out_a,
      ProgA = r_prog_a / r_out_a, DropA = r_drop_a / r_out_a, DetA = r_det_a / r_out_a,
      DetA_Pub = DetA,
      
      r_drop_s = r_death_s + r_sc,
      r_prog_s = r_csi * (1 - p_txi0),
      r_det_s = r_acf_tp + r_csi * p_txi0,
      r_out_s = r_drop_s + r_prog_s + r_det_s,
      Dur_S = 1 / r_out_s,
      ProgS = r_prog_s / r_out_s, DropS = r_drop_s / r_out_s, DetS = r_det_s / r_out_s,
      DetS_Pub = (r_acf_tp + r_csi * p_txi0_pub) / r_out_s,
      DetS_Eng = r_csi * p_txi0_eng / r_out_s,
      DetS_Pri = r_csi * p_txi0_pri / r_out_s,
      
      r_drop_c = r_death_s + r_sc,
      r_det_c = r_acf_tp + r_recsi * p_txi1,
      r_out_c = r_drop_c + r_det_c,
      Dur_C = 1 / r_out_c,
      DropC = r_drop_c / r_out_c, DetC = r_det_c / r_out_c,
      DetC_Pub = (r_acf_tp + r_recsi * p_txi1_pub) / r_out_c,
      DetC_Eng = r_recsi * p_txi1_eng / r_out_c,
      DetC_Pri = r_recsi * p_txi1_pri / r_out_c
    ) %>% 
    mutate(
      End_SC = DropA * r_sc / r_drop_a,
      End_SC = End_SC + ProgA * DropS * r_sc / r_drop_s,
      End_SC = End_SC + ProgA * ProgS * DropC * r_sc / r_drop_c,
      End_Dead = DropA * r_death_a / r_drop_a,
      End_Dead = End_Dead + ProgA * DropS * r_death_s / r_drop_s,
      End_Dead = End_Dead + ProgA * ProgS * DropC * r_death_s / r_drop_c,
      
      pd_pub = DetA_Pub + ProgA * DetS_Pub + ProgA * ProgS * DetC_Pub,
      pd_eng = ProgA * DetS_Eng + ProgA * ProgS * DetC_Eng,
      pd_pri = ProgA * DetS_Pri + ProgA * ProgS * DetC_Pri, 
      
      pd_a = DetA,
      pd_s = ProgA * DetS,
      pd_c = ProgA * ProgS * DetC,
      
      End_TC = pd_pub * p_txs_pub + pd_eng * p_txs_eng + pd_pri * p_txs_pri,
      End_TD = pd_pub * p_txd_pub + pd_eng * p_txd_eng + pd_pri * p_txd_pri,
      End_TL = pd_pub * p_txl_pub + pd_eng * p_txl_eng + pd_pri * p_txl_pri,
      
      Gap_A = (1 - DropA),
      Gap_S = (1 - DropS),
      Gap_C = (1 - DropC),
      Gap_T = End_TC / (End_TC + End_TD + End_TL),
      
      Cas_A = Gap_A,
      Cas_S = Cas_A * Gap_S,
      Cas_C = Cas_S * Gap_C,
      Cas_T = Cas_C * Gap_T,
      
      Delay_Pat = Dur_S,
      Delay_Sym = Dur_C * DetS / (DetS + ProgS),
      Delay_Tot = Delay_Pat + Delay_Sym,
      Location = loc
    ) %>% 
    select(Location, starts_with(c("Dur_", "Delay_", "End_", "Gap_", "Cas_"), ignore.case=F))
  
  
}))



delays <- res %>% 
  group_by(Location) %>% 
  summarise(across(starts_with("Delay_"), list(
    Mu = mean,
    M = median,
    L = function(x) quantile(x, 0.025),
    U = function(x) quantile(x, 0.975)
  ))) %>% 
  mutate(
    across(starts_with("Delay"), function(x) x * 12)
  )

g_pat <- delays %>% 
  ggplot() + 
  geom_pointrange(aes(y = Location, x = Delay_Pat_M, xmin = Delay_Pat_L, xmax = Delay_Pat_U)) +
  scale_x_continuous("Patient Delay, mo.", breaks = c(0, 3, 6, 12))


g_sym <- delays %>% 
  ggplot() + 
  geom_pointrange(aes(y = Location, x = Delay_Sym_M, xmin = Delay_Sym_L, xmax = Delay_Sym_U)) +
  scale_x_continuous("System Delay, mo.", breaks = c(0, 3, 6, 12))



ggsave(g_pat, filename = here::here("docs", "figs", "g_pat.png"), width = 6, height = 10)
ggsave(g_sym, filename = here::here("docs", "figs", "g_sys.png"), width = 6, height = 10)






