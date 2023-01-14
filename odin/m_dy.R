deriv(Asym) <- inc - (r_sym + r_sc + r_death_a + r_acf * sens_acf) * Asym

deriv(Sym) <- r_sym * Asym - (r_csi + r_sc + r_death_s + r_acf * sens_acf) * Sym

deriv(ExCS) <- fn0 - sum(txi1) - (r_sc + r_death_s + r_acf * sens_acf) * ExCS

deriv(TxPubPub) <- txi0[1] + txi1[1] + det_acf - (1 / dur_pub) * TxPubPub
deriv(TxEngPub) <- (txi0[2] + txi1[2]) * (1 - p_pri_on_pub) - (1 / dur_pub) * TxEngPub
deriv(TxEngPri) <- (txi0[2] + txi1[2]) *      p_pri_on_pub  - (1 / dur_pri) * TxEngPri
deriv(TxPriPri) <- txi0[3] + txi1[3] - (1 / dur_pri) * TxPriPri

deriv(EndDead) <- r_death_a * Asym + r_death_s * (Sym + ExCS)
deriv(EndSC) <- r_sc * (Asym + Sym + ExCS)
deriv(EndTxD) <- (p_txd[1] / dur_pub) * TxPubPub + (p_txd[2] / dur_pub) * TxEngPub + (p_txd[2] / dur_pri) * TxEngPri + (p_txd[3] / dur_pri) * TxPriPri
deriv(EndTxL) <- (p_txl[1] / dur_pub) * TxPubPub + (p_txl[2] / dur_pub) * TxEngPub + (p_txl[2] / dur_pri) * TxEngPri + (p_txl[3] / dur_pri) * TxPriPri
deriv(EndTxS) <- (p_txs[1] / dur_pub) * TxPubPub + (p_txs[2] / dur_pub) * TxEngPub + (p_txs[2] / dur_pri) * TxEngPri + (p_txs[3] / dur_pri) * TxPriPri
  

deriv(FpPubPub) <- 0
deriv(FpEngPub) <- 0
deriv(FpEngPri) <- 0
deriv(FpPriPri) <- 0



## Statistics
output(PrevA) <- Asym
output(PrevS) <- Sym
output(PrevC) <- ExCS

PrevUt <- Asym + Sym + ExCS
PrevTx <- TxPubPub + TxEngPub + TxEngPri + TxPriPri
N <- PrevUt + PrevTx

output(PrevUt) <- TRUE
output(PrevTx) <- TRUE
output(N) <- TRUE

output(N_All) <- N + EndDead + EndSC + EndTxD + EndTxS + EndTxL

output(PrA) <- Asym / PrevUt
output(PrS) <- Sym / PrevUt
output(PrC) <- ExCS / PrevUt


## Initial values
initial(Asym) <- prv_a
initial(Sym) <- prv_s
initial(ExCS) <- prv_c

initial(TxPubPub) <- 0
initial(TxEngPub) <- 0
initial(TxEngPri) <- 0
initial(TxPriPri) <- 0

initial(EndDead) <- 0
initial(EndSC) <- 0
initial(EndTxD) <- 0
initial(EndTxL) <- 0
initial(EndTxS) <- 0

initial(FpPubPub) <- 0
initial(FpEngPub) <- 0
initial(FpEngPri) <- 0
initial(FpPriPri) <- 0

prv_a <- user(1)
prv_s <- user(0)
prv_c <- user(0)



## Calculation

det0[] <- r_csi * entry[i] * p_dx0[i] * Sym
dim(det0) <- 3

txi0[] <- det0[i] * p_txi[i]
dim(txi0) <- 3
  

det1[] <- r_recsi * entry[i] * p_dx1[i] * ExCS
dim(det1) <- 3

txi1[] <- det1[i] * p_txi[i]
dim(txi1) <- 3


fn0 <- r_csi * Sym - sum(txi0)


det_acf <- r_acf * (sens_acf * PrevUt)

output(det_acf) <- TRUE
output(det_pub) <- det0[1] + det1[1]
output(det_eng) <- det0[2] + det1[2]
output(det_pri) <- det0[3] + det1[3]
output(det_all) <- sum(det0) + sum(det1) + det_acf


# inc <- r_death_a * PrevA + r_death_s * (PrevS + PrevC) + sum(det0) + sum(det1) - (r_sc + r_acf * sens_acf + adr) * prv
inc <- (r_sym + r_sc + r_death_a + r_acf * sens_acf - adr) * Asym * has_dy



## Load parameters

has_dy <- user(1)


r_sc <- user(0.2)
r_death_a <- user(0)
r_death_s <- user(0.1)

r_sym <- user(1)
r_csi <- user(1)
r_recsi <- user(1)
dur_pub <- user(0.5)
dur_pri <- user(0.5)

r_acf <- user(0.0001)
sens_acf <- user(0.8)
spec_acf <- user(0.995)



p_dx0[1] <- p_dx0_pub
p_dx0[2] <- p_dx0_eng
p_dx0[3] <- p_dx0_pri
dim(p_dx0) <- 3

p_dx1[1] <- p_dx1_pub
p_dx1[2] <- p_dx1_eng
p_dx1[3] <- p_dx1_pri
dim(p_dx1) <- 3

p_dx0_pub <- user(0.5)
p_dx0_eng <- user(0.5)
p_dx0_pri <- user(0.5)
p_dx1_pub <- user(0.5)
p_dx1_eng <- user(0.5)
p_dx1_pri <- user(0.5)


entry[1] <- p_ent_pub
entry[2] <- p_ent_eng
entry[3] <- p_ent_pri
dim(entry) <- 3

p_ent_pub <- user(0.3)
p_ent_eng <- user(0.3)
p_ent_pri <- user(0.4)


ppv[1] <- ppv_pub
ppv[2] <- ppv_eng
ppv[3] <- ppv_pri
dim(ppv) <- 3

ppv_pub <- user(1)
ppv_eng <- user(1)
ppv_pri <- user(1)


adr <- user(0.01)

p_pri_on_pub <- user(0.3)

p_txi[1] <- p_txi_pub
p_txi[2] <- p_txi_eng
p_txi[3] <- p_txi_pri
dim(p_txi) <- 3

p_txi_pub <- user(0.3)
p_txi_eng <- user(0.3)
p_txi_pri <- user(0.4)


p_txs[1] <- p_txs_pub
p_txs[2] <- p_txs_eng
p_txs[3] <- p_txs_pri
dim(p_txs) <- 3

p_txs_pub <- user(0.3)
p_txs_eng <- user(0.3)
p_txs_pri <- user(0.4)


p_txl[1] <- p_txl_pub
p_txl[2] <- p_txl_eng
p_txl[3] <- p_txl_pri
dim(p_txl) <- 3

p_txl_pub <- user(0.3)
p_txl_eng <- user(0.3)
p_txl_pri <- user(0.4)


p_txd[1] <- p_txd_pub
p_txd[2] <- p_txd_eng
p_txd[3] <- p_txd_pri
dim(p_txd) <- 3

p_txd_pub <- user(0.3)
p_txd_eng <- user(0.3)
p_txd_pri <- user(0.4)

