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
  fn0 <- d$r_cs * entry * (1 - pdi0) * d$Prev_C
  
  det1 <- c(d$DetR_pub, d$DetR_eng, d$DetR_pri) - det0
  
  fn <- function(r_recsi, d, trm, pdi1, fn0) {
    prev_e <- solve(a = (r_recsi * t(trm * (1 - pdi1)) - diag(d$r_mu_sym + r_recsi, 3)), b = -fn0)
    (sum(prev_e) - d$Prev_E) ^ 2
  }
  
  
  r_recsi <- nlminb(1, fn, lower = 0, d = d, trm = trm, pdi1 = pdi1, fn0 = fn0)$par
  r_recsi
}



extract_referrals <- function(trm, loc) {
  require(tidyverse)
  
  mx <- rbind(trm[1:3, ], pub_det = 0, eng_det = 0, pri_det = 0)
  mx[4:6, 4:6] <- diag(1, 3)
  
  
  tr_1 <- mx * trm["UC", ]
  
  
  mx_end <- mx
  for (i in 1:20) {
    mx_end <- mx_end %*% mx
  }
  mx_end[, 1:3] <- 0
  mx_end <- mx_end / rowSums(mx_end)
  tr_end <- mx_end * colSums(tr_1)
  
  
  flows <- bind_rows(
    data.frame(as.table(tr_1)) %>% 
      tibble() %>% 
      rename(Sector0 = Var1, Sector1 = Var2) %>% 
      mutate(Stage0 = "Entry", Stage1 = "Second"),
    data.frame(as.table(tr_end)) %>% 
      tibble() %>% 
      rename(Sector0 = Var1, Sector1 = Var2) %>% 
      mutate(Stage0 = "Second", Stage1 = "Diagnosis"),
  )%>% 
    mutate(
      Sector0 = as.character(Sector0),
      Sector0 = case_when(
        Sector0 == "pub_det" ~ "pub",
        Sector0 == "eng_det" ~ "eng",
        Sector0 == "pri_det" ~ "pri",
        T ~ Sector0
      ),
      Sector1 = as.character(Sector1),
      Sector1 = case_when(
        Sector1 == "pub_det" ~ "pub",
        Sector1 == "eng_det" ~ "eng",
        Sector1 == "pri_det" ~ "pri",
        T ~ Sector1
      )
    ) %>% 
    group_by(Sector0, Sector1, Stage0, Stage1) %>% 
    summarise(Freq = sum(Freq)) %>% 
    ungroup() %>% 
    arrange(Stage0, Sector0) %>% 
    mutate(
      Stage0 = factor(Stage0, c("Entry", "Second", "Diagnosis")),
      Stage1 = factor(Stage1, c("Entry", "Second", "Diagnosis")),
      Sector0 = factor(Sector0, c("pri", "eng", "pub")),
      Sector1 = factor(Sector1, c("pri", "eng", "pub"))
    ) %>% 
    arrange(Stage0, Stage1, Sector0, Sector1)
  
  
  
  stocks <- bind_rows(
    flows %>% 
      select(Sector = Sector0, Stage = Stage0, Freq),
    flows %>% 
      select(Sector = Sector1, Stage = Stage1, Freq) %>% 
      filter(Stage == "Diagnosis")
  ) %>% 
    group_by(Stage, Sector) %>% 
    summarise(Freq = sum(Freq)) %>% 
    ungroup() %>% 
    mutate(
      Stage = factor(Stage, c("Entry", "Second", "Diagnosis")),
      Sector = factor(Sector, c("pri", "eng", "pub")),
    ) %>% 
    arrange(Stage, Sector)
  
  list(
    stocks = stocks,
    flows = flows,
    Location = loc
  )
}



vis_referrals <- function(stocks, flows, bar.width=20, interval=70, n.step=50) {
  require(tidyverse)
  
  sts.n <- 3
  width <- sts.n * bar.width + (sts.n-1) * interval
  
  
  stocks <- stocks %>% 
    mutate(
      x0 = (as.numeric(Stage) - 1) * interval,
      x1 = (as.numeric(Stage) - 1) * interval + bar.width
    ) %>% 
    group_by(Stage) %>% 
    mutate(
      y1 = cumsum(Freq),
      y0 = c(0, y1[-length(y1)])
    )
  
  
  labels_x <- stocks %>% 
    select(Stage, x0, x1) %>% 
    mutate(x = (x0 + x1) / 2) %>% 
    distinct()
  
  
  bands <- flows %>% 
    group_by(Stage0) %>%
    arrange(Stage0, Sector0, Sector1) %>% 
    mutate(
      y1s = cumsum(Freq),
      y0s = c(0, y1s[-length(y1s)])
    ) %>% 
    ungroup() %>% 
    group_by(Stage1) %>% 
    arrange(Stage1, Sector1, Sector0) %>% 
    mutate(
      y1t = cumsum(Freq),
      y0t = c(0, y1t[-length(y1t)])
    ) %>% 
    ungroup() %>% 
    left_join(stocks %>% 
                select(Stage0 = Stage, x0 = x1) %>% 
                distinct()) %>% 
    left_join(stocks %>% 
                select(Stage1 = Stage, x1 = x0) %>% 
                distinct()) %>% 
    ungroup() %>% 
    mutate(Key = 1:n()) %>% 
    merge(tibble(xx = seq(-pi/2, pi/2, length.out = n.step))) %>% 
    arrange(Key) %>% 
    group_by(Key) %>% 
    mutate(
      ys.upper = y0s + (y0t-y0s)/2 * (sin(xx) + 1 ),
      ys.lower = y1s + (y1t-y1s)/2 * (sin(xx) + 1 ),
      xs = seq(x0[1], x1[1], length.out = n.step)
    )
  
  g <- ggplot(stocks) +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = Sector), 
              colour = "grey7") +
    geom_ribbon(data=bands, aes(x=xs, ymin=ys.lower, ymax=ys.upper, group = Key), 
                colour = "grey7", fill="darkgreen", alpha=0.3) +
    scale_x_continuous("Stage", breaks = labels_x$x, 
                       labels = c("Initial visit", "Second visit", "Diagnosis")) +
    scale_fill_discrete("Sector", 
                        labels = c(pri = "Private", eng = "Engaged Private", pub = "Public"))
  
  list(
    stocks = stocks,
    flows = flows,
    bands = bands,
    g = g
  )
}

