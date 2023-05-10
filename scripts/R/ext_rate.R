ext_rate <- function(df, exo, settings) {
  require(tidyverse)
  
  with(settings, {
    df %>% 
      merge(tibble(Key = exo$Key)) %>% 
      mutate(
        PrevUt = rbinom(n(), prob = PrevUt / 1e5, size = N) / N,
        PrevTx = rbinom(n(), prob = PrevTx / 1e5, size = N) / N
      ) %>% 
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
        r_ltfu_pri = r_succ_pri / TxSucc_eng * TxLTFU_eng * rr_ltfu_pri,
        r_dead_pri = r_succ_pri / TxSucc_eng * TxDead_eng,
        
        Prev_Tx_Pub = TxR_pub / (r_succ_pub + r_ltfu_pub + r_dead_pub),
        Prev_Tx_Eng = TxR_eng / (r_succ_eng + r_ltfu_eng + r_dead_eng),
        Prev_Tx_Pri = PrevTx - Prev_Tx_Pub - Prev_Tx_Eng,
        Prev_Tx_Pri = ifelse(Prev_Tx_Pri < 0, 0, Prev_Tx_Pri),
        Prev_Tx_Pri = Prev_Tx_Pri * ppv_pri,
        TxR_pri = Prev_Tx_Pri * (r_succ_pri + r_ltfu_pri + r_dead_pri),
        
        DetR_pub = CNR_pub,
        DetR_eng = CNR_eng,
        DetR_pri = TxR_pri / TxI_pri,
        DetR = DetR_pub + DetR_eng + DetR_pri,
        Prev_A = PrevUt * Pr_Asym,
        Prev_S = PrevUt * Pr_NotAware,
        Prev_C = PrevUt * Pr_NotCS,
        Prev_E = PrevUt * Pr_NotDet,
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
      select(Location, Key, Pop=N, starts_with("DetR"), starts_with("TxR"), starts_with("TxI_"), 
             PrevUt, PrevTx, Prev_A, Prev_S, Prev_C, Prev_E, starts_with("r_"), starts_with("Pr_Pub_"))
  })
}
