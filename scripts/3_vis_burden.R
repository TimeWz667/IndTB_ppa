library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


post <- bind_rows(lapply(c("tx_00", "tx_01", "tx_11"), function(folder) {
  read_csv(here::here("out", folder, "post.csv")) %>% 
    select(ppm, dur_pri, ppv_pri, p_pri_on_pub, 
           tp_pri_drug, tp_pri_drug_time, tp_pri_txi) %>% 
    mutate(Scenario = folder)
}))


g_attload <- post %>% 
  extract(Scenario, c("TBPS", "DrugSale"), "_(\\d)(\\d)", remove=F) %>% 
  pivot_longer(starts_with("tp_pri")) %>% 
  ggplot() +
  stat_halfeye(aes(x = value, y = name, fill = Scenario), position = position_dodge(-0.2), alpha=0.4) +
  #facet_grid(TBPS~DrugSale, labeller = label_both) +
  scale_x_continuous("Caseloads, millions", labels=scales::number_format(scale = 1e-6), limits = c(0, 5e6)) +
  scale_y_discrete("", labels=c(
    tp_pri_drug = "Txi with private drugs",
    tp_pri_drug_time = "Caseloads with private drugs",
    tp_pri_txi = "Txi in the unengaged private"
  )) +
  scale_fill_discrete("Data sets", labels=c(tx_00="Nikshay only", tx_01="+Drug sales", tx_11="+TBPS"))

ggsave(g_attload, filename=here::here("doc"))



