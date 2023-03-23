library(tidyverse)


## Data loading
p_txi <- read_csv(here::here("docs", "tabs", "post_txi.csv"))
p_txi

load(here::here("data", "dat_cas.rdata"))


dat <- local({
  inp <- p_txi %>% 
    filter(State == "India") %>% 
    select(p_txi_pub = p_txi_pub_m, p_txi_eng = p_txi_eng_m) %>% 
    as.list()
  
  d <- dat_noti %>% 
    filter((Year == 2020) & State == "India") %>%
    select(Pop, starts_with("N_")) %>% 
    as.list()
  
  inp <- c(inp, d)
  
  
  d <- dat_tbps %>% 
    filter(State == "India") %>% 
    mutate(
      n_att = N_OnATT_Pri + N_OnATT_Pub,
      p_pub_m = N_OnATT_Pub / n_att,
      p_pub_l = qbinom(0.025, size = n_att, prob = p_pub_m) / n_att,
      p_pub_u = qbinom(0.975, size = n_att, prob = p_pub_m) / n_att
    )
  
  tar <- d %>% 
    select(starts_with("p_pub_"), starts_with("DrugTime_")) %>% 
    as.list()
  
  lower <- c(tar$p_pub_l, tar$DrugTime_L)
  upper <- c(tar$p_pub_u, tar$DrugTime_U)
  
  inp$Drug <- d$DrugTime_M
  inp$Drug_Std <- (d$DrugTime_U - d$DrugTime_L) / 2 / 1.96
    
  inp$Tx_Pub <- round(d$N_OnATT_Pub)
  inp$Tx <- round(d$n_att)

  
  list(
    Inputs = inp,
    Targets = tar,
    Lower = lower,
    Upper = upper
  )
})


dat



### Fitting ----
model <- rstan::stan_model(here::here("stan", "ppv_txdur.stan"))

ds <- dat$Inputs
ds$ppv_pub <- 0.85

post <- rstan::sampling(model, data=ds, iter=5e3, warmup=4e3)

plot(post)

tab <- as.data.frame(summary(post)$summary)

tab

write.csv(tab, here::here("docs", "tabs", "ppm.csv"))


xx = rstan::extract(post, pars = c("ppm", "dur_pri", "ppv_pri")) %>% 
  data.frame() %>% 
  as_tibble()


plot(xx)




