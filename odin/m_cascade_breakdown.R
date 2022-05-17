
deriv(PrevA) <- - (r_onset + r_sc + r_die_asym) * PrevA 

deriv(PrevS) <- r_onset * PrevA - (r_aware + r_sc + r_die_sym) * PrevS

deriv(PrevC) <- r_aware * PrevS - (r_csi + r_sc + r_die_sym) * PrevC

deriv(PrevE[]) <- (fn0[i] + sum(fn1[, i])) - (r_recsi + r_sc + r_die_sym) * PrevE[i]
dim(PrevE) <- 3

deriv(PrevTx[]) <- (det0[i] + sum(det1[, i])) * p_txi[i] - (r_ltfu[i] + r_succ[i] + r_die_tx[i]) * PrevTx[i]
dim(PrevTx) <- 3


deriv(Cas0) <- r_onset * PrevA 
deriv(Cas1) <- r_aware * PrevS
deriv(Cas2) <- r_csi * PrevC 
deriv(Cas3_0) <- sum(det0)
deriv(Cas3_1) <- sum(det1)
deriv(Cas3) <- sum(det0) + sum(det1)

deriv(Cas4) <- sum(det0) + sum(det1) - sum(lost_txi)
deriv(Cas5) <- r_succ[1] * PrevTx[1] + r_succ[2] * PrevTx[2] + r_succ[3] * PrevTx[3]


deriv(Drop0_dead) <- r_die_asym * PrevA
deriv(Drop0_sc) <- r_sc * PrevA
deriv(Drop1_dead) <- r_die_sym * PrevS
deriv(Drop1_sc) <- r_sc * PrevS
deriv(Drop2_dead) <- r_die_sym * PrevC
deriv(Drop2_sc) <- r_sc * PrevC
deriv(Drop3_dead) <- r_die_sym * sum(PrevE)
deriv(Drop3_sc) <- r_sc * sum(PrevE)
deriv(Drop4_ltfu) <- sum(lost_txi)
deriv(Drop5_dead) <- r_die_tx[1] * PrevTx[1] + r_die_tx[2] * PrevTx[2] + r_die_tx[3] * PrevTx[3]
deriv(Drop5_ltfu) <- sum(lost_txd)


Drop_dead <- Drop0_dead + Drop1_dead + Drop2_dead + Drop3_dead + Drop5_dead
Drop_sc <- Drop0_sc + Drop1_sc + Drop2_sc + Drop3_sc
Drop_ltfu <- Drop4_ltfu + Drop5_ltfu

output(Drop_dead) <- TRUE
output(Drop_sc) <- TRUE
output(Drop_ltfu) <- TRUE
output(End) <- Drop_dead + Drop_sc + Drop_ltfu + Cas5

initial(PrevA) <- 1
initial(PrevS) <- 0
initial(PrevC) <- 0
initial(PrevE[]) <- 0  
initial(PrevTx[]) <- 0


initial(Cas0) <- 0
initial(Cas1) <- 0
initial(Cas2) <- 0
initial(Cas3_0) <- 0
initial(Cas3_1) <- 0
initial(Cas3) <- 0
initial(Cas4) <- 0
initial(Cas5) <- 0

initial(Drop0_dead) <- 0
initial(Drop0_sc) <- 0
initial(Drop1_dead) <- 0
initial(Drop1_sc) <- 0
initial(Drop2_dead) <- 0
initial(Drop2_sc) <- 0
initial(Drop3_dead) <- 0
initial(Drop3_sc) <- 0
initial(Drop4_ltfu) <- 0
initial(Drop5_dead) <- 0
initial(Drop5_ltfu) <- 0


det0[] <- r_csi * entry[i] * p_diag0[i] * PrevC
dim(det0) <- 3

fn0[] <- r_csi * entry[i] * (1 - p_diag0[i]) * PrevC
dim(fn0) <- 3

det1[, ] <- r_recsi * trm_shift[i, j] * p_diag1[i, j] * PrevE[i]
dim(det1) <- c(3, 3)

fn1[, ] <- r_recsi * trm_shift[i, j] * (1 - p_diag1[i, j]) * PrevE[i]
dim(fn1) <- c(3, 3)


lost_txi[] <- (det0[i] + sum(det1[, i])) * (1 - p_txi[i])
dim(lost_txi) <- 3
lost_txd[] <- r_ltfu[i] * PrevTx[i]
dim(lost_txd) <- 3


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
