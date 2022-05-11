calc_cascade <- function(df) {
  require(tidyverse)
  
  df %>% 
    mutate(
      r_g0 = r_onset + r_sc + r_die_asym,
      Gap0_prog = r_onset / r_g0,
      Gap0_sc = r_sc / r_g0,
      Gap0_dead = r_die_asym / r_g0,
      
      r_g1 = r_aware + r_sc + r_die_sym,
      Gap1_prog = r_aware / r_g1, 
      Gap1_sc = r_sc / r_g1,
      Gap1_dead = r_die_sym / r_g1,
      
      r_g2 = r_cs + r_sc + r_die_sym,
      Gap2_prog = r_cs / r_g2, 
      Gap2_sc = r_sc / r_g2,
      Gap2_dead = r_die_sym / r_g2,
      
      r_g3 = r_det + r_sc + r_die_sym,
      Gap3_prog = r_det / r_g3, 
      Gap3_sc = r_sc / r_g3,
      Gap3_dead = r_die_sym / r_g3,
      
      Gap4_prog = TxR / DetR,
      Gap4_ltfu = 1 - Gap4_prog,
      
      ent_pub = DetR_pub * TxI_pub,
      ent_eng = DetR_eng * TxI_eng,
      ent_pri = DetR_pri * TxI_pri,
      across(c(ent_pub, ent_eng, ent_pri), function(x) x / TxR),
      
      r_g5_pub = r_succ_pub + r_ltfu_pub + r_dead_pub,    
      Gap5_succ_pub = r_succ_pub / r_g5_pub,
      Gap5_ltfu_pub = r_ltfu_pub / r_g5_pub,
      Gap5_dead_pub = r_dead_pub / r_g5_pub,
      
      r_g5_eng = r_succ_eng + r_ltfu_eng + r_dead_eng,    
      Gap5_succ_eng = r_succ_eng / r_g5_eng,
      Gap5_ltfu_eng = r_ltfu_eng / r_g5_eng,
      Gap5_dead_eng = r_dead_eng / r_g5_eng,
      
      r_g5_pri = r_succ_pri + r_ltfu_pri + r_dead_pri,    
      Gap5_succ_pri = r_succ_pri / r_g5_pri,
      Gap5_ltfu_pri = r_ltfu_pri / r_g5_pri,
      Gap5_dead_pri = r_dead_pri / r_g5_pri,
      
      Gap5_succ = ent_pub * Gap5_succ_pub + ent_eng * Gap5_succ_eng + ent_pri * Gap5_succ_pri,
      Gap5_ltfu = ent_pub * Gap5_ltfu_pub + ent_eng * Gap5_ltfu_eng + ent_pri * Gap5_ltfu_pri,
      Gap5_dead = ent_pub * Gap5_dead_pub + ent_eng * Gap5_dead_eng + ent_pri * Gap5_dead_pri,
      
      Cas0 = Gap0_prog,
      Cas1 = Cas0 * Gap1_prog,
      Cas2 = Cas1 * Gap2_prog,
      Cas3 = Cas2 * Gap3_prog,
      Cas4 = Cas3 * Gap4_prog,
      Cas5 = Cas4 * Gap5_succ,
      
      Drop0_dead = Gap0_dead, 
      Drop0_sc = Gap0_sc,
      Drop1_dead = Cas0 * Gap1_dead, 
      Drop1_sc = Cas0 * Gap1_sc,
      Drop2_dead = Cas1 * Gap2_dead, 
      Drop2_sc = Cas1 * Gap2_sc,
      Drop3_dead = Cas2 * Gap3_dead, 
      Drop3_sc = Cas2 * Gap3_sc,
      Drop4_ltfu = Cas3 * Gap4_ltfu, 
      Drop5_ltfu = Cas4 * Gap5_ltfu,
      Drop5_dead = Cas4 * Gap5_dead,
      
      DurAsym = 1 / r_g0,
      DurNA = 1 / r_g1,
      DurNC = 1 / r_g2,
      DurCS = 1 / r_g3,
      
      TTSym = DurAsym,
      TTAware = TTSym + DurNA,
      TTCare = TTAware + DurNC,
      TTDet = TTCare + DurCS,
      
      DelayPatient = DurNA + DurNC,
      DelaySystem = DurCS,
      DelayTotal = DelayPatient + DelaySystem,
      
      Burden = TxR * Pop
    ) %>% 
    select(
      Location, Key,
      Burden,
      starts_with("Gap"), starts_with("Cas"), starts_with("Dur"), starts_with("TT"),
      starts_with("Delay"), starts_with("Drop"))
  
}