library(tidyverse)


# Data loading
pop <- read_csv(here::here("data", "Population.csv"))

state_map1 <- read_csv(here::here("data", "StateMap1.csv"))
state_map2 <- read_csv(here::here("data", "StateMap2.csv"))
state_map3 <- read_csv(here::here("data", "StateMap3.csv"))


tbps <- read_csv(here::here("data", "TBPS", "TBPS_ASC_State.csv"))
drug <- read_csv(here::here("data", "DrugSale.csv"))

load(here::here("data", "d_itr.rdata"))



# Collect data mapping


state_map <- state_map1 %>% 
  left_join(state_map2) %>% 
  left_join(state_map3) %>% 
  relocate(State, State_Pop, State_Itr, StateGroup, Region)


# Population data

d_pop <- pop %>% 
  filter(Sex == "Total") %>% 
  select(Year, State_Pop = Location, Pop)


# For treatment initiation analysis

dat_txi <- state_map %>% 
  left_join(d_itr_notif %>% rename(State_Itr = State)) %>%
  select(State, Region, Year, N_Det_Pub = N_Noti_Pub, N_Det_Eng = N_Noti_Pri, N_Txi_Pub, N_Txi_Eng = N_Txi_Pri) %>% 
  filter(N_Txi_Pub > 0)

save(dat_txi, file = here::here("data", "dat_txi.rdata"))



# For treatment outcome analysis
dat_tx <- state_map %>% 
  left_join(d_itr_tx %>% rename(State_Itr = State)) %>% 
  select(State, Region, Year, starts_with("N_")) %>% 
  rename(N_Tx_Ini_Eng = N_Tx_Ini_Pri, N_Tx_Succ_Eng = N_Tx_Succ_Pri, N_Tx_Die_Eng = N_Tx_Die_Pri)

save(dat_tx, file = here::here("data", "dat_tx.rdata"))


# For cascade analysis
dat_noti <- state_map %>% 
  left_join(d_itr_notif %>% rename(State_Itr = State)) %>%
  left_join(d_itr_acf %>% rename(State_Itr = State)) %>% 
  left_join(d_pop) %>% 
  mutate(
    N_Det_ACF = N_ACF_Detected,
    N_Det_Pub = N_Noti_Pub - N_Det_ACF
  ) %>% 
  select(State, Region, Year, Pop, N_Det_Pub, N_Det_Eng = N_Noti_Pri, N_Det_ACF, N_Txi_Pub, N_Txi_Eng = N_Txi_Pri)


state_map %>% 
  left_join(d_itr_notif %>% rename(State_Itr = State)) %>%
  left_join(d_itr_acf %>% rename(State_Itr = State)) %>% 
  left_join(d_pop) %>% 
  mutate(
    N_Det_ACF = N_ACF_Detected,
    N_Det_Pub = N_Noti_Pub - N_Det_ACF,
    Prop_Screened = N_ACF_Screened / Pop,
    Prop_Mapped = N_ACF_Mapped / Pop
  ) %>% 
  filter(Year == 2021) %>% 
  ggplot() +
  geom_point(aes(x = Prop_Screened, y = State, colour = "Screened")) +
  geom_point(aes(x = Prop_Mapped, y = State, colour = "Mapped")) + 
  geom_vline(xintercept = 0.3) + 
  scale_x_continuous("/overall population", labels = scales::percent)


dat_tbps <- state_map %>% 
  left_join(tbps %>% 
              select(State, N, N_Asym, N_NotCS, N_NotDet, N_OnATT_Pub, N_OnATT_Pri)) %>% 
  left_join(drug %>% 
              extract(Tx_month_pri, c("TxMo_M", "TxMo_L", "TxMo_U"), "(\\S+) \\((\\S+), (\\S+)\\)", convert = T) %>% 
              select(State, starts_with("TxMo"))) %>% 
  mutate(Year = 2020) %>% 
  select(State, Region, N, starts_with("N_"), starts_with("TxMo_"))


save(dat_noti, dat_tbps, file = here::here("data", "dat_cas.rdata"))



state_map %>% 
  left_join(d_itr_acf %>% rename(State_Itr = State)) %>% 
  left_join(d_pop) %>% 
  select(State, Region, P1, Year, Pop, starts_with("N_ACF")) %>% 
  mutate(
    N_TB = P1 / 1e5 * Pop,
    N_TP = N_TB * 0.7,
    N_FP = (Pop - N_TB) * (1 - 0.99),
    PPV = N_TP / (N_TP + N_FP)
  ) %>% 
  ggplot() +
  geom_line(aes(x = Year, y = PPV)) +
  facet_wrap(.~State)


(x * 0.7 + (1 - x) * 0.01) *  0.5 = 0.7 * x
0.69 * 0.5 * x + 0.01 * 0.5 = 0.7 * x

0.01 * 0.5 / (1 - 0.99 * 0.5)




  





