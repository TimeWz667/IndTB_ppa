extract_shifting <- function(amp_pri, pdi_pri, p_eng, adj_pub, mat0, det) {
  require(tidyverse)
  
  mat <- mat0[c("pub", "pub", "pri", "UC"), c("pub", "pri", "pri", "pub_det")]
  colnames(mat)[2] <- "eng"
  rownames(mat)[2] <- "eng"
  mat <- cbind(mat, eng_det = 0, pri_det = 0)
  mat[, c("eng", "pri")] <- mat[, c("eng", "pri")] * amp_pri
  
  # mat["pub", "pub_det"] <- mat["pub", "pub_det"] * x[1]
  
  p_diag_pub <- mat[, "pub_det"] / (mat[, "pub_det"] + mat[, "pub"])
  odds <- p_diag_pub / (1 - p_diag_pub) * adj_pub
  p_diag_pub <- odds / (1 + odds)
  
  p_diag <- matrix(0, 4, 3)
  colnames(p_diag) <- c("pub", "eng", "pri")
  rownames(p_diag) <- c("pub", "eng", "pri", "UC")
  
  p_diag[, "pub"] <- p_diag_pub[c("pub", "pub", "pri", "UC")]
  p_diag[, "eng"] <- p_diag[, "pub"]
  
  odds <- p_diag[, "pub"] / (1 - p_diag[, "pub"]) * pdi_pri
  p_diag[, "pri"] <- odds / (1 + odds)
  
  pub <- mat[, "pub"] + mat[, "pub_det"]
  mat[, "pub"] <- pub * (1 - p_diag[, "pub"])
  mat[, "pub_det"] <- pub * p_diag[, "pub"]
  
  k <- p_eng * (1 - p_diag[, "eng"]) + (1 - p_eng) * (1 - p_diag[, "pri"])
  mat[, "eng"] <- mat[, "eng"] * p_eng * (1 - p_diag[, "eng"]) / k
  mat[, "pri"] <- mat[, "pri"] * (1 - p_eng) * (1 - p_diag[, "pri"]) / k
  mat[, "eng_det"] <- mat[, "eng"] * (p_diag[, "eng"]) / (1 - p_diag[, "eng"])
  mat[, "pri_det"] <- mat[, "pri"] * (p_diag[, "pri"]) / (1 - p_diag[, "pri"])
  mat <- mat / rowSums(mat)
  
  
  tr_mat <- mat
  
  m0 <- list(vis = 1)
  m0 <- c(m0, tr_mat["UC", ])
  
  ms <- list(m0)
  
  for (i in 2:20) {
    m <- list(vis = i)
    m <- c(m, colSums(tr_mat[c("pub", "eng", "pri"), ] * c(m0$pub, m0$eng, m0$pri)))
    
    ms[[length(ms) + 1]] <- m
    m0 <- m
  }
  
  ms <- bind_rows(ms)
  
  return(list(
    tr_mat = tr_mat,
    p_diag = p_diag,
    dy = ms
  ))
}


find_r_recsi <- function(d, sim_shifting) {
  pdi <- sim_shifting$p_diag
  shift <- sim_shifting$tr_mat[, 1:3] + sim_shifting$tr_mat[, 4:6]
  
  entry <- shift["UC", ]
  pdi0 <- pdi["UC", ]
  
  trm <- shift[1:3, ]
  pdi1 <- pdi[1:3, ]
  
  det0 <- d$r_cs * entry * pdi0 * d$Prev_C
  fn0 <- d$r_cs * entry * (1 -pdi0) * d$Prev_C
  
  det1 <- c(d$DetR_pub, d$DetR_eng, d$DetR_pri) - det0
  
  fn <- function(r_recsi, d, trm, pdi1, fn0) {
    prev_e <- solve(a = (r_recsi * t(trm * (1 - pdi1)) - diag(d$r_mu_sym + r_recsi, 3)), b = -fn0)
    (sum(prev_e) - d$Prev_E) ^ 2
  }
  
  
  r_recsi <- nlminb(1, fn, lower = 0, d = d, trm = trm, pdi1 = pdi1, fn0 = fn0)$par
  r_recsi
}




