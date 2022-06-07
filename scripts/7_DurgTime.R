library(tidyverse)




load(here::here("data", "cascade", "d_cascade_2019.rdata"))


drug <- read_csv(here::here("data", "cascade0", "TBPS.csv")) %>% 
  select(State = Location, TxPriPub) %>%
  mutate(
    State = gsub("_", " ", State)
  )
  filter(State %in% c("India", "Delhi", "Kerala", "Tamil Nadu", "Gujarat", "Uttar Pradesh"))



d_cas <- d_cascade %>% 
  select(Region, State, PrevTx, starts_with("CNR"), starts_with("Tx")) %>%
  inner_join(drug) %>% 
  crossing(
    PPV_Pri = seq(0.05, 1, 0.05),
    DurTx_pri = seq(2 / 12, 1, 1 / 24)
  ) %>% 
  mutate(
    TxI_pri = 0.8,
    PPV_pub = 0.85,
    PPV_eng = 0.85,
    rr_ltfu_pri = 1.5,
    r_succ_pub = 2,
    r_succ_eng = 2,
    r_succ_pri = 1 / DurTx_pri,
    r_ltfu_pub = r_succ_pub / TxSucc_pub * TxLTFU_pub,
    r_dead_pub = r_succ_pub / TxSucc_pub * TxDead_pub,
    r_ltfu_eng = r_succ_eng / TxSucc_eng * TxLTFU_eng,
    r_dead_eng = r_succ_eng / TxSucc_eng * TxDead_eng,
    r_ltfu_pri = r_succ_pri / TxSucc_eng * TxLTFU_eng * rr_ltfu_pri,
    r_dead_pri = r_succ_pri / TxSucc_eng * TxDead_eng,
    dur_tx_pub = 1 / (r_succ_pub + r_ltfu_pub + r_dead_pub),
    dur_tx_eng = 1 / (r_succ_eng + r_ltfu_eng + r_dead_eng),
    dur_tx_pri = 1 / (r_succ_pri + r_ltfu_pri + r_dead_pri),
    TxR_pub = CNR_pub * TxI_pub,
    TxR_eng = CNR_eng * TxI_eng,
    PrevTx_pub = TxR_pub * dur_tx_pub,
    PrevTx_eng = TxR_eng * dur_tx_eng,
    PrevTx_pri_TBPS = PrevTx * 1e-5 - PrevTx_pub - PrevTx_eng,
    PrevTx_pri_TBPS = pmax(PrevTx_pri_TBPS, 0),
    PrevTx_pri_drug = PrevTx_pub * TxPriPub - PrevTx_eng,
    PrevTx_pri_drug = pmax(PrevTx_pri_drug, 0)
  )

areas <- unique(d_cas$Region)
areas <- setNames(areas, areas)


gs <- lapply(areas, function(area) {
  d_cas %>% 
    filter(Region == area) %>% 
    ggplot() + 
    geom_point(aes(x = PrevTx_pri_TBPS, y = PrevTx_pri_drug, colour = State)) +
    geom_abline(aes(slope = 1, intercept = 0)) + 
    scale_x_continuous("Pr(Tx, private) with TBPS, per 100k", labels = scales::number_format(scale = 1e5)) +
    scale_y_continuous("Pr(Tx, private) with Drug sale data, per 100k", labels = scales::number_format(scale = 1e5)) +
    expand_limits(x = c(0, 0.005), y = c(0, 0.005)) +
    labs(subtitle = area)
})


g <- ggpubr::ggarrange(gs$Central, gs$Eastern, gs$Northern, gs$Southern, gs$Western, nrow = 3, ncol = 2)
g

ggsave(g, filename = here::here("docs", "TBPS_DrugSale.jpg"), height=10, width=10)







