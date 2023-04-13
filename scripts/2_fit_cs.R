library(tidyverse)
library(rstan)

options(mc.cores = 4)
rstan_options(auto_write = TRUE)


library(tidyverse)
library(rstan)


options(mc.cores = 4)
rstan_options(auto_write = TRUE)


## Data loading

targets <- read_csv(here::here("data", "targets_india.csv"))

targets %>% 
  filter(Year == 2019) %>% data.frame()


drug <- targets %>% filter(Index == "DrugTime")
tx <- targets %>% filter(Index == "PrTxiPub") %>% 
  mutate(X = round(M * N))

det <- targets %>% 
  filter(Year == 2019) %>% 
  filter(Index == "TxI") %>% 
  mutate(
    N_Txi = round(N * M),
    N_Det = N
  )


ds <- list(
  N_Det_Pub = det$N_Det[1],
  N_Det_Eng = det$N_Det[2],
  N_Txi_Pub = det$N_Txi[1],
  N_Txi_Eng = det$N_Txi[2],
  Pop = 1383112050,
  Tx = tx$N,
  Tx_Pub = tx$X,
  Drug = drug$M,
  Drug_Std = drug$Error,
  ppv_pub = 0.85,
  ent_pub = 0.483,
  pdx_pub = 0.65
)



### Fitting ----

# "cs_00.stan", "cs_10.stan", "cs_01.stan"

for(src_model in c("cs_11.stan")) {
  model <- rstan::stan_model(here::here("stan", src_model))
  
  post <- rstan::sampling(model, data=ds, iter=5e3, warmup=4e3)
  
  tab <- as.data.frame(summary(post)$summary)
  tab$Name <- rownames(tab)
  tab <- tab %>% as_tibble() %>% relocate(Name)
  
  ext <- rstan::extract(post) %>% 
    data.frame() %>% 
    as.tibble() %>% 
    mutate(
      ppv_pub = ds$ppv_pub,
      ent_pub = ds$entpub,
      pdx_pub = ds$pdx_pub
    )
  
  folder <- gsub(".stan", "", src_model)
  dir.create(here::here("out", folder), showWarnings = F)
  write_csv(tab, file = here::here("out", folder, "summary.csv"))
  write_csv(ext, file = here::here("out", folder, "post.csv"))
}





