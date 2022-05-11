calc_rate <- function(df, exo, settings) {
  require(tidyverse)
  
  with(settings, {
    df %>% 
      mutate(
        lam = (Prev / Prev_Std) ^ 2,
        k = Prev / lam
      ) %>%
      merge(tibble(Key = exo$Key)) %>% 
      mutate(
        Prev = rpois(n(), lam) * k
      ) %>% 
      select(-lam, -k, -Prev_Std) %>% 
      left_join(exo) %>% 
      arrange(Location, Key) %>% 
      relocate(Location, Key) %>%  
      mutate(
        TxI_pri = TxI_pri,
        TxR_pub = CNR_pub * TxI_pub,
        TxR_eng = CNR_eng * TxI_eng,
        r_succ_pub = r_succ_pub,
        r_succ_eng = r_succ_eng,
        r_succ_pri = r_succ_pri,
        r_ltfu_pub = r_succ_pub / TxSucc_pub * TxLTFU_pub,
        r_dead_pub = r_succ_pub / TxSucc_pub * TxDead_pub,
        r_ltfu_eng = r_succ_eng / TxSucc_eng * TxLTFU_eng,
        r_dead_eng = r_succ_eng / TxSucc_eng * TxDead_eng,
        r_ltfu_pri = r_succ_pri / TxSucc_eng * TxLTFU_eng,
        r_dead_pri = r_succ_pri / TxSucc_eng * TxDead_eng,
        dur_tx_pub = 1 / (r_succ_pub + r_ltfu_pub + r_dead_pub),
        dur_tx_eng = 1 / (r_succ_eng + r_ltfu_eng + r_dead_eng),
        dur_tx_pri = 1 / (r_succ_pri + r_ltfu_pri + r_dead_pri),
        time_tx_pub = TxR_pub * dur_tx_pub,
        time_tx_eng = TxR_eng * dur_tx_eng,
        time_tx_pri = TxPriPub * time_tx_pub / TxTime_pri2pub - time_tx_eng,
        TxR_pri = time_tx_pri / dur_tx_eng,
        TxR = TxR_pub + TxR_eng + TxR_pri,
        DetR_pub = CNR_pub,
        DetR_eng = CNR_eng,
        DetR_pri = TxR_pri / TxI_pri,
        DetR = DetR_pub + DetR_eng + DetR_pri,
        Prev_A = Prev * Pr_Asym,
        Prev_Sym = Prev * (1 - Pr_Asym),
        Prev_S = Prev_Sym * (1 - Pr_CS_Sym) * (Pr_NA_NC),
        Prev_C =Prev_Sym * (1 - Pr_CS_Sym) * (1 - Pr_NA_NC),
        Prev_E = Prev_Sym * Pr_CS_Sym,
        r_sc = r_sc,
        r_die_asym = r_die_untx * rr_die_asym,
        r_die_sym = r_die_untx,
        r_mu_asym = r_sc + r_die_asym,
        r_mu_sym = r_sc + r_die_sym,
        r_det = DetR / Prev_E,
        r_cs = Prev_E / Prev_C * (r_mu_sym + r_det),
        r_aware = Prev_C / Prev_S * (r_mu_sym + r_cs),
        r_onset = Prev_S / Prev_A * (r_mu_sym + r_aware)
      ) %>% 
      select(Location, Key, starts_with("DetR"), starts_with("TxR"), starts_with("TxI_"), 
             Prev, Prev_A, Prev_S, Prev_C, Prev_E, starts_with("r_"))
  })
}


calc_rate_shared_r_onset <- function(df) {
  df %>% 
    group_by(Key) %>% 
    mutate(
      r_onset = r_onset[Location == "India"],
      Prev = Prev_A + Prev_S + Prev_C + Prev_E,
      Prev_SCE0 = Prev_S + Prev_C + Prev_E,
      Prev_A = (DetR + r_mu_sym * Prev) / (r_onset + r_mu_sym),
      Prev_S = Prev_S / Prev_SCE0 * (Prev - Prev_A),
      Prev_C = Prev_C / Prev_SCE0 * (Prev - Prev_A),
      Prev_E = Prev_E / Prev_SCE0 * (Prev - Prev_A),
      r_aware = r_onset * Prev_A / Prev_S - r_mu_sym,
      r_cs = r_aware * Prev_S / Prev_C - r_mu_sym,
      r_det = r_cs * Prev_C / Prev_E - r_mu_sym, 
      Det2 = Prev_E * r_det,
      PrAsym = Prev_A / Prev,
      PNratio = Prev / DetR
    ) %>% 
    filter(r_aware > 0 & r_cs > 0 & r_det > 0) %>% 
    select(-Prev_SCE0) %>% 
    ungroup()
}



