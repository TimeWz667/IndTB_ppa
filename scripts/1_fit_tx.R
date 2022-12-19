library(tidyverse)
library(rstan)


dir.create("out/sub_tx")


load(here::here("data", "dat_tx.rdata"))


m_single <- rstan::stan_model(here::here("stan", "txi_single.stan"))
m_multiple <- rstan::stan_model(here::here("stan", "txi_multiple.stan"))





sel_india <- dat_tx %>% filter(State == "India")


die_india <- list(
  N = nrow(sel_india),
  D_Pub = sel_india$N_Tx_Ini_Pub,
  D_Eng = sel_india$N_Tx_Ini_Eng,
  T_Pub = sel_india$N_Tx_Die_Pub,
  T_Eng = sel_india$N_Tx_Die_Eng
)


succ_india <- list(
  N = nrow(sel_india),
  D_Pub = sel_india$N_Tx_Ini_Pub - sel_india$N_Tx_Die_Pub,
  D_Eng = sel_india$N_Tx_Ini_Eng - sel_india$N_Tx_Die_Eng,
  T_Pub = sel_india$N_Tx_Succ_Pub,
  T_Eng = sel_india$N_Tx_Succ_Eng
)



sel_state <- dat_tx %>% filter(State != "India") %>%
  arrange(State) %>% 
  mutate(
    State = factor(State),
    Region = factor(Region)
  )



sel_state_pred <- sel_state %>% select(State, Region) %>% distinct()


die_state <- list(
  N = nrow(sel_state),
  D_Pub = sel_state$N_Tx_Ini_Pub,
  D_Eng = sel_state$N_Tx_Ini_Eng,
  T_Pub = sel_state$N_Tx_Die_Pub,
  T_Eng = sel_state$N_Tx_Die_Eng,
  States = as.numeric(sel_state$State),
  Regions = as.numeric(sel_state$Region),
  N_State = length(levels(sel_state$State)),
  N_Region = length(levels(sel_state$Region)),
  
  States_Pred = as.numeric(sel_state_pred$State),
  Regions_Pred = as.numeric(sel_state_pred$Region)
)


succ_state <- list(
  N = nrow(sel_state),
  D_Pub = sel_state$N_Tx_Ini_Pub - sel_state$N_Tx_Die_Pub,
  D_Eng = sel_state$N_Tx_Ini_Eng - sel_state$N_Tx_Die_Eng,
  T_Pub = sel_state$N_Tx_Succ_Pub,
  T_Eng = sel_state$N_Tx_Succ_Eng,
  States = as.numeric(sel_state$State),
  Regions = as.numeric(sel_state$Region),
  N_State = length(levels(sel_state$State)),
  N_Region = length(levels(sel_state$Region)),
  
  States_Pred = as.numeric(sel_state_pred$State),
  Regions_Pred = as.numeric(sel_state_pred$Region)
)


states <- sel_state_pred$State


## Run fitting
post_txd_india <- sampling(m_single, data=die_india, pars=c("p_pub", "p_eng"))
post_txs_india <- sampling(m_single, data=succ_india, pars=c("p_pub", "p_eng"))

post_txd_state <- sampling(m_multiple, data=die_state, pars=c("p_pub", "p_eng"))
post_txs_state <- sampling(m_multiple, data=succ_state, pars=c("p_pub", "p_eng"))


save(post_txd_india, post_txs_india, post_txd_state, post_txs_state,
     file = here::here("out", "sub_tx", "post.rdata"))


## Extract to json
pd_pub <- extract(post_txd_india, "p_pub")[[1]]
pd_eng <- extract(post_txd_india, "p_eng")[[1]]

ps_pub <- extract(post_txs_india, "p_pub")[[1]]
ps_eng <- extract(post_txs_india, "p_eng")[[1]]

pars <- lapply(1:length(pd_pub), function(i) {
  list(
    p_txd_pub = pd_pub[i],
    p_txs_pub = (1 - pd_pub[i]) * ps_pub[i],
    p_txl_pub = (1 - pd_pub[i]) * (1 - ps_pub[i]),
    p_txd_eng = pd_eng[i],
    p_txs_eng = (1 - pd_eng[i]) * ps_eng[i],
    p_txl_eng = (1 - pd_eng[i]) * (1 - ps_eng[i])
  )
})

jsonlite::write_json(pars, here::here("out", "sub_tx", "pars_India.json"), simplifyVector=T, auto_unbox=T)


tab <- bind_rows(pars) %>% 
  mutate(State = "India")



for (sti in 1:length(states)) {
  state <- states[sti]
  
  pd_pub <- extract(post_txd_state, "p_pub")[[1]][, sti]
  pd_eng <- extract(post_txd_state, "p_eng")[[1]][, sti]
  
  ps_pub <- extract(post_txs_state, "p_pub")[[1]][, sti]
  ps_eng <- extract(post_txs_state, "p_eng")[[1]][, sti]
  
  pars <- lapply(1:length(pd_pub), function(i) {
    list(
      p_txd_pub = pd_pub[i],
      p_txs_pub = (1 - pd_pub[i]) * ps_pub[i],
      p_txl_pub = (1 - pd_pub[i]) * (1 - ps_pub[i]),
      p_txd_eng = pd_eng[i],
      p_txs_eng = (1 - pd_eng[i]) * ps_eng[i],
      p_txl_eng = (1 - pd_eng[i]) * (1 - ps_eng[i])
    )
  })
  
  tab <- bind_rows(tab, 
                   bind_rows(pars)%>% mutate(State = state)) 
  
  jsonlite::write_json(pars, here::here("out", "sub_tx", "pars_" + glue::as_glue(state) + ".json"), simplifyVector=T, auto_unbox=T)
}



## Summarise
tab <- tab %>% 
  group_by(State) %>% 
  summarise(across(everything(), list(
    m = mean,
    l = function(x) quantile(x, 0.025),
    u = function(x) quantile(x, 0.975)
  )))


write_csv(tab, here::here("docs", "tabs", "post_tx.csv"))



