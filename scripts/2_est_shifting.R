library(tidyverse)


load(here::here("data", "shifting", "shifting.rdata"))


tr_mat <- shifting %>% 
  group_by(From, To) %>% 
  summarise(N = sum(N)) %>% 
  group_by(From) %>% 
  mutate(
    Pr = N / sum(N),
    To = ifelse(To == "det", "pub_det", To)
  )

save(tr_mat, file = here::here("data", "shifting", "shifting_mat.rdata"))
