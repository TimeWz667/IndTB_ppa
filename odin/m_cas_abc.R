deriv(PrevA) <- inc - (r_onset + r_sc + r_death_a + r_acf * sens_acf) * PrevA

deriv(PrevS) <- r_onset * PrevA - (r_csi + r_sc + r_death_s + r_acf * sens_acf) * PrevS

deriv(PrevC) <- fn0 - sum(det1) - (r_sc + r_death_s + r_acf * sens_acf) * PrevC

deriv(PrevTx) <- sum(txi) - PrevTx / dur_tx

N <- PrevA + PrevS + PrevC + PrevTx

output(N) <- TRUE


dur_a <- 1 / (r_onset + r_sc + r_death_a + r_acf * sens_acf)
dur_s <- 1 / (r_csi + r_sc + r_death_s + r_acf * sens_acf)
dur_c <- 1 / (r_det1 + r_sc + r_death_s + r_acf * sens_acf)


p_det_a <- (r_acf * sens_acf) * dur_a
p_prog_a <- r_onset * dur_a
p_drop_a <- 1 - p_det_a - p_prog_a

p_det_s <- (sum(det0) + r_acf * sens_acf) * dur_s
p_prog_s <- fn0 / PrevS  * dur_s
p_drop_s <- 1 - p_det_s - p_prog_s

r_det1 <- sum(det1) / PrevC

p_det_c <- (r_det1 + r_acf * sens_acf) * dur_c
p_drop_c <- 1 - p_det_c


output(dur_a) <- TRUE
output(dur_s) <- TRUE
output(dur_c) <- TRUE

output(p_drop_a) <- TRUE
output(p_drop_s) <- TRUE
output(p_drop_c) <- TRUE


del_pat <- dur_s 
del_sys <- dur_c * p_prog_s * p_det_c / (p_prog_s * p_det_c + p_det_a)
output(del_pat) <- TRUE
output(del_sys) <- TRUE
output(del_tot) <- del_pat + del_sys


initial(PrevA) <- prv_a
initial(PrevS) <- prv_s
initial(PrevC) <- prv_c
initial(PrevTx) <- prv_t

prv_a <- user(0)
prv_s <- user(0)
prv_c <- user(0)
prv_t <- user(0)

prv <- PrevA + PrevS + PrevC

output(PrA) <- PrevA / prv
output(PrS) <- PrevS / prv
output(PrC) <- PrevC / prv


det0[] <- r_csi * PrevS * entry[i] * p_dx0[i]
dim(det0) <- 3


det1[] <- r_recsi * PrevC * entry[i] * p_dx1[i]
dim(det1) <- 3

fn0 <- r_csi * PrevS - sum(det0)
  
txi[] <- (det0[i] + det1[i]) / ppv[i]
dim(txi) <- 3


det_acf <- r_acf * (sens_acf * (PrevA + PrevS + PrevC) + spec_acf * (1 - N))

output(det_acf) <- TRUE
output(det_pub) <- det0[1] + det1[1]
output(det_eng) <- det0[2] + det1[2]
output(det_pri) <- det0[3] + det1[3]
output(det_all) <- sum(det0) + sum(det1)


# inc <- r_death_a * PrevA + r_death_s * (PrevS + PrevC) + sum(det0) + sum(det1) - (r_sc + r_acf * sens_acf + adr) * prv
inc <- (r_onset + r_sc + r_death_a + r_acf * sens_acf - adr) * PrevA

output(cdr) <- (r_acf * sens_acf * prv + det0[1] + det1[1] + det0[2] + det1[2]) / inc

r_sc <- user(0.2)
r_death_a <- user(0)
r_death_s <- user(0.1)

r_onset <- user(1)
r_csi <- user(1)
r_recsi <- user(1)
dur_tx <- user(0.5)

r_acf <- user(0.0001)
sens_acf <- user(0.8)
spec_acf <- user(0.995)


intv <- user(0)
intv_t <- if (t >= 2023) 1 - intv else 1
output(intv_t) <- TRUE

p_dx0[1] <- 1 - (1 - p_dx0_pub) * intv_t
p_dx0[2] <- 1 - (1 - p_dx0_eng) * intv_t
p_dx0[3] <- p_dx0_pri
dim(p_dx0) <- 3

p_dx1[1] <- 1 - (1 - p_dx1_pub) * intv_t
p_dx1[2] <- 1 - (1 - p_dx1_eng) * intv_t
p_dx1[3] <- p_dx1_pri
dim(p_dx1) <- 3

p_dx0_pub <- user(0.5)
p_dx0_eng <- user(0.5)
p_dx0_pri <- user(0.5)
p_dx1_pub <- user(0.5)
p_dx1_eng <- user(0.5)
p_dx1_pri <- user(0.5)


ppv[1] <- ppv_pub
ppv[2] <- ppv_eng
ppv[3] <- ppv_pri
dim(ppv) <- 3

ppv_pub <- user(1)
ppv_eng <- user(1)
ppv_pri <- user(1)


entry[1] <- p_pub
entry[2] <- p_eng
entry[3] <- p_pri
dim(entry) <- 3

p_pub <- user(0.3)
p_eng <- user(0.3)
p_pri <- user(0.4)

adr <- user(0.01)

