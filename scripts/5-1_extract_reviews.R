library(tidyverse)
library(rstan)


dur <- read_csv(here::here("data", "review", "sys_dur.csv")) %>% 
  select(Identifier, N, M, L, U)


pr <- read_csv(here::here("data", "review", "sys_pr.csv")) %>% 
  select(Identifier, CutOff, Low, High)


d_dur <- dur %>% 
  select(N_dur = N, M_dur = M) %>% 
  as.list()

d_dur$nd <- nrow(dur)


d_pr <- pr %>% 
  mutate(N_pr = Low + High, X_pr = Low, X_cutoff = CutOff) %>% 
  select(N_pr, X_pr, X_cutoff) %>% 
  as.list()

d_pr$np <- nrow(pr)



dat <- c(d_dur, d_pr)


dat



mod <- rstan::stan_model(here::here("stan", "meta_coxian.stan"))


n_iter <- 5000
n_warmup <- 4000

post <- sampling(mod, data=dat, iter=n_iter, warmup=n_warmup)















