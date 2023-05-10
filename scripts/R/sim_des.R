sample_des <- function(i, inp) {
  endpts <- c("SelfCure", "DieUt", "DieTx", "Det_pub", "Det_eng", "Det_pri")
  
  with(inp, {
    state <- "Start"
    pathway <- state
    ti <- 0
    tis <- ti
    visit <- 0
    delay_pat <- NA
    delay_sys <- NA
    sector0 <- NA
    sector1 <- NA
    
    
    trans <- c(SelfCure = rexp(1, r_sc), 
               DieUt = rexp(1, r_die_sym),
               Aware = rexp(1, r_aware))
    
    state <- names(which.min(trans))
    ti <- min(trans) + ti
    tis <- c(tis, ti)
    pathway <- c(pathway, state)
    
    
    if (!(state %in% endpts)) {
      trans <- c(SelfCure = rexp(1, r_sc), 
                 DieUt = rexp(1, r_die_sym),
                 Seeking = rexp(1, r_csi))
      
      state <- names(which.min(trans))
      ti <- min(trans) + ti
      tis <- c(tis, ti)
      
      if (state == "Seeking") {
        state <- sample(names(entry), 1, prob = entry)
        sector0 <- state
        if (runif(1) < p_diag0[state]) {
          sector1 <- state
          state <- paste0("Det_", state) 
        }
        delay_pat <- ti
      }
      pathway <- c(pathway, state)
    }
    
    
    while (!(state %in% endpts)) {
      ent <- trm_shift[state, ]
      pdi <- p_diag1[state, ]
      
      trans <- c(SelfCure = rexp(1, r_sc), 
                 DieUt = rexp(1, r_die_sym),
                 Seeking = rexp(1, r_recsi))
      
      state <- names(which.min(trans))
      ti <- min(trans) + ti
      tis <- c(tis, ti)
      
      if (state == "Seeking") {
        state <- sample(names(ent), 1, prob = ent)
        
        if (runif(1) < pdi[state]) {
          sector1 <- state
          state <- paste0("Det_", state) 
        }
      }
      pathway <- c(pathway, state)
    }
    
    if (startsWith(state, "Det")) {
      delay_sys = ti - delay_pat
    }
    
    
    list(
      ID = i, 
      Pathway = paste0(pathway, collapse = ":"),
      Tis = paste0(scales::number(tis, .001), collapse = ":"),
      DelayPatient = delay_pat,
      DelaySystem = delay_sys,
      End = state,
      SectorStart = sector0,
      SectorEnd = sector1
    )
  })
}
