library(tidyverse)



d_cas <- read_csv(here::here("data", "cascade", "d_cascade_2019.csv"))

d_drug <- read_csv(here::here("data/drug_sale.csv"))


d_cas %>% 
  mutate(
    DurTxPub = PrevTxPub / (CNR_pub * TxI_pub) 
  ) %>% 
  ggplot() +
  geom_point(aes(x = DurTxPub, y = reorder(State, DurTxPub))) +
  scale_y_discrete("") +
  scale_x_continuous("Treatment duration")



d_cas %>% 
  select(State, Region, CNR_pub, CNR_eng, TxI_pub, TxI_eng, starts_with("Prev")) %>% 
  mutate(
    RatioCNR = CNR_pub / CNR_eng,
    RatioPrevTx = PrevTxPub / PrevTxPri
  ) %>% 
  ggplot() +
  geom_point(aes(x = RatioCNR, y = RatioPrevTx))


d_cas %>% 
  mutate(
    DurTxPri = PrevTxPri / (CNR_eng * TxI_eng)
  ) %>% 
  ggplot() +
  geom_point(aes(x = DurTxPri, y = reorder(State, DurTxPri))) +
  scale_y_discrete("") +
  scale_x_continuous("Treatment duration", limits = c(0, 2))


d_cas %>% 
  mutate(
    DurTxPri = PrevTxPri / (CNR_eng * TxI_eng),
    DurTxPub = PrevTxPub / (CNR_pub * TxI_pub)
  ) %>% 
  ggplot() +
  geom_point(aes(x = DurTxPri, y = DurTxPub)) +
  scale_y_continuous(limits = c(0, 2)) +
  scale_x_continuous(limits = c(0, 2)) +
  geom_text(aes(x = DurTxPri, y = DurTxPub, label = State), hjust = 0, vjust = 1)+
  scale_x_continuous("Treatment duration, private", limits = c(0, 2))+
  scale_y_continuous("Treatment duration, public", limits = c(0, 2))




d_cas %>% 
  mutate(
    DurTxPub = PrevTxPub / (CNR_pub * TxI_pub),
    PrevPri = PrevTxPri - (CNR_eng * TxI_eng * DurTxPub)
  ) %>% 
  ggplot() + 
  geom_point(aes(x = PrevPri, y = reorder(State, PrevPri))) +
  geom_vline(xintercept = 0)




d_cas %>% 
  mutate(
    DurTxPub = PrevTxPub / (CNR_pub * TxI_pub),
    PrevPri = PrevTxPri - (CNR_eng * TxI_eng * DurTxPub)
  ) %>% 
  left_join(d_drug %>% 
              extract(Tx_month_pri, c("M", "L", "U"), "(\\S+) \\((\\S+), (\\S+)\\)", convert = T) %>% 
              mutate(M = M / 12, U = U / 12, L = L / 12) %>% 
              select(State, M, L, U)
  ) %>% 
  ggplot() +
  geom_point(aes(x = M, y = PrevPri * 1e5)) + 
  geom_text(aes(x = M, y = PrevPri * 1e5, label = State), hjust = 0, vjust = 1) +
  geom_hline(yintercept = 0) +
  geom_abline(slope = 1) +
  expand_limits(x = 0) + 
  scale_x_continuous("Prev. on private drug from drug-sale data") +
  scale_y_continuous("Prev. on private drug from TBPS data")







