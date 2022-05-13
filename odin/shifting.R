deriv(y[]) <- fn0[i] + sum(fn[, i]) - (r_recsi + mu) * y[i]
initial(y[]) <- 0
dim(y) = 3


det[, ] <- r_recsi * y[i] * trm[i, j] * p_diag1[i, j]
dim(det) <- c(3, 3)

fn[, ] <- r_recsi * y[i] * trm[i, j] * (1 - p_diag1[i, j])
dim(fn) <- c(3, 3)


output(prev_e) <- sum(y)
output(det_pub) <- sum(det[, 1])
output(det_eng) <- sum(det[, 2])
output(det_pri) <- sum(det[, 3])
output(det_all) <- sum(det)


fn0[] <- user()
dim(fn0) <- 3

trm[, ] <- user()
dim(trm) <- c(3, 3)

p_diag1[, ] <- user()
dim(p_diag1) <- c(3, 3)

mu <- user(0.3)
r_recsi <- user(2)