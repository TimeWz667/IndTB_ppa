library(tidyverse)
library(rstan)



load(here::here("data", "dat_txi.rdata"))




sel_india <- dat_txi %>% filter(State == "India")


sel_india <- list(
  N = nrow(sel_india),
  D_Pub = sel_india$N_Det_Pub,
  D_Eng = sel_india$N_Det_Eng,
  T_Pub = sel_india$N_Txi_Pub,
  T_Eng = sel_india$N_Txi_Eng
)



m_single <- rstan::stan_model(here::here("stan", "txi_single.stan"))
m_multiple <- rstan::stan_model(here::here("stan", "txi_multiple.stan"))


sel_state <- dat_txi %>% filter(State != "India") %>%
  filter(N_Det_Pub > N_Txi_Pub) %>% 
  filter(N_Det_Eng > N_Txi_Eng) %>% 
  arrange(State) %>% 
  mutate(
    State = factor(State),
    Region = factor(Region)
  )


sel_state_pred <- sel_state %>% select(State, Region) %>% distinct()


sel_state <- list(
  N = nrow(sel_state),
  D_Pub = sel_state$N_Det_Pub,
  D_Eng = sel_state$N_Det_Eng,
  T_Pub = sel_state$N_Txi_Pub,
  T_Eng = sel_state$N_Txi_Eng,
  States = as.numeric(sel_state$State),
  Regions = as.numeric(sel_state$Region),
  N_State = length(levels(sel_state$State)),
  N_Region = length(levels(sel_state$Region)),
  
  States_Pred = as.numeric(sel_state_pred$State),
  Regions_Pred = as.numeric(sel_state_pred$Region)
)


states <- sel_state_pred$State


## Run fitting
post_txi_india <- sampling(m_single, data=sel_india, pars=c("p_pub", "p_eng"))

post_txi_state <- sampling(m_multiple, data=sel_state, pars=c("p_pub", "p_eng"))


save(post_txi_india, post_txi_state, file = here::here("out", "sub_txi", "post.rdata"))


## Extract to json

p_pub <- extract(post_txi_india, "p_pub")[[1]]
p_eng <- extract(post_txi_india, "p_eng")[[1]]

pars <- lapply(1:dim(p_pub)[1], function(i) list(p_txi_pub = p_pub[i], p_txi_eng = p_eng[i]))

jsonlite::write_json(pars, here::here("out", "sub_txi", "pars_India.json"), simplifyVector=T, auto_unbox=T)


for (sti in 1:length(states)) {
  p_pub <- extract(post_txi_state, "p_pub")[[1]][, sti]
  p_eng <- extract(post_txi_state, "p_eng")[[1]][, sti]
  
  pars <- lapply(1:length(p_pub), function(i) list(p_txi_pub = p_pub[i], p_txi_eng = p_eng[i]))
  
  jsonlite::write_json(pars, here::here("out", "sub_txi", "pars_" + glue::as_glue(states[sti]) + ".json"), simplifyVector=T, auto_unbox=T)
}



## Summarise parameter values
temp <- summary(post_txi_state)$summary
tab_state <- tibble(name=rownames(temp), m=temp[, "mean"], l=temp[, "2.5%"], u=temp[, "97.5%"]) %>% 
  filter(name != "lp__") %>% 
  tidyr::extract(name, c('Index', "State"), "(\\S+)\\[(\\d+)\\]", convert = T) %>% 
  mutate(
    State = states[State],
    State = as.character(State),
    Index = ifelse(Index == "p_pub", "p_txi_pub", "p_txi_eng")
  ) %>% 
  pivot_longer(c(m, l, u)) %>% 
  pivot_wider(names_from = Index, values_from = value) %>% 
  pivot_wider(values_from = c(p_txi_pub, p_txi_eng))


temp <- summary(post_txi_india)$summary
tab_all <- tibble(Index=rownames(temp), m=temp[, "mean"], l=temp[, "2.5%"], u=temp[, "97.5%"]) %>% 
  filter(Index != "lp__") %>% 
  mutate(State = "India",
         Index = ifelse(Index == "p_pub", "p_txi_pub", "p_txi_eng")) %>% 
  pivot_longer(c(m, l, u)) %>% 
  pivot_wider(names_from = Index, values_from = value) %>% 
  pivot_wider(values_from = c(p_txi_pub, p_txi_eng))



tab <- bind_rows(tab_all, tab_state)

write_csv(tab, here::here("docs", "tabs", "post_txi.csv"))
