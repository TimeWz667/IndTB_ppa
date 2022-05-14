
deriv(PrevA) <- inc - (r_onset + r_sc + r_die_asym) * PrevA 

deriv(PrevS) <- r_onset * PrevA - (r_aware + r_sc + r_die_sym) * PrevS

deriv(PrevC) <- r_aware * PrevS - (r_csi + r_sc + r_die_sym) * PrevC

deriv(PrevE[]) <- (fn0[i] + sum(fn1[, i])) - (r_recsi + r_sc + r_die_sym) * PrevE[i]
dim(PrevE) <- 3

deriv(PrevTx[]) <- (det0[i] + sum(det1[, i])) * p_txi[i] - (r_ltfu[i] + r_succ[i] + r_die_tx[i]) * PrevTx[i]
dim(PrevTx) <- 3



inc <- (r_sc + r_die_asym) * PrevA + (r_sc + r_die_sym) * (PrevS + PrevC + sum(PrevE)) + sum(det0) + sum(det1)


output(Prev) <- PrevA + PrevS + PrevC + sum(PrevE) + sum(PrevTx)

initial(PrevA) <- prev0
initial(PrevS) <- 0
initial(PrevC) <- 0
initial(PrevE[]) <- 0  
initial(PrevTx[]) <- 0

prev0 <- user()


det0[] <- r_csi * entry[i] * p_diag0[i] * PrevC
dim(det0) <- 3

fn0[] <- r_csi * entry[i] * (1 - p_diag0[i]) * PrevC
dim(fn0) <- 3

det1[, ] <- r_recsi * trm_shift[i, j] * p_diag1[i, j] * PrevE[i]
dim(det1) <- c(3, 3)

fn1[, ] <- r_recsi * trm_shift[i, j] * (1 - p_diag1[i, j]) * PrevE[i]
dim(fn1) <- c(3, 3)


lost[] <- (det0[i] + sum(det1[, i])) * (1 - p_txi[i]) + r_ltfu[i] * PrevTx[i]
dim(lost) <- 3


output(det_pub) <- sum(det1[, 1]) + det0[1]
output(det_eng) <- sum(det1[, 2]) + det0[2]
output(det_pri) <- sum(det1[, 3]) + det0[3]
output(det_all) <- sum(det1) + sum(det0)


r_sc <- user(0.2)
r_die_asym <- user()
r_die_sym <- user()

r_onset <- user()
r_aware <- user()

r_csi <- user()

entry[] <- user()
dim(entry) <- 3

p_diag0[] <- user()
dim(p_diag0) <-3


r_recsi <- user(2)

trm_shift[, ] <- user()
dim(trm_shift) <- c(3, 3)

p_diag1[, ] <- user()
dim(p_diag1) <- c(3, 3)


p_txi[] <- user()
dim(p_txi) <- 3

r_succ[] <- user()
dim(r_succ) <- 3
r_ltfu[] <- user()
dim(r_ltfu) <- 3
r_die_tx[] <- user()
dim(r_die_tx) <- 3
