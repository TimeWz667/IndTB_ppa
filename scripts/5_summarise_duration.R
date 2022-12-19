library(tidyverse)



## Setup

theme_set(theme_bw() + theme(text = element_text(family = "sans")))

ext <- glue::as_glue(".png")

ds <- dir(here::here("out", "sim", "cas"))


lvs <- c(
  "S_Asym", "S_Sym", "S_ExCs", "S_TxPub", "S_TxPri", 
  "S_Cured", "S_SelfCured", "S_LTFU", "S_DeadU", "S_DeadT"
)



# TTE
tt_labels <- c(
  TTS = "Sympton Developed",
  TTC = "Care sought", TTD = "Diagnosed"
)


for (tp in c("onds", "nods")) {
  tp = glue::as_glue(tp)
  
  
  sim <- bind_rows(lapply(ds[startsWith(ds, "pars_" + tp)], 
                          function(file) read.csv(here::here("out", "sim", "cas", file)))) %>% 
    extract(Source, "Location", "pars_" + tp + "_(\\S+).csv")
  
  
  tte <- sim %>% 
    group_by(Location, Key) %>% 
    mutate(dt = diff(Time)[1]) %>% 
    summarise(across(starts_with("S_"), function(x) sum(x * dt))) %>% 
    select(Location, Key, S_Asym, S_Sym, S_ExCs) %>% 
    ungroup() %>% 
    mutate(
      TTS = S_Asym,
      TTC = TTS + S_Sym,
      TTD = TTC + S_ExCs
    ) %>% 
    pivot_longer(starts_with("TT"), names_to = "Stage") %>% 
    group_by(Location, Stage) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.025),
      U = quantile(value, 0.975)
    ) %>% 
    ungroup()
  
  
  tt_rnk <- tte %>% 
    filter(Stage == "TTD") %>% 
    arrange(M) %>% 
    pull(Location) %>% 
    as.character()
  
  
  g_tte <- tte %>% 
    mutate(
      Stage = factor(Stage, rev(names(tt_labels))),
      Location = factor(Location, tt_rnk)
    ) %>% 
    arrange(Location) %>% 
    ggplot() +
    geom_rect(xmin = -Inf, xmax = Inf, 
              ymin = which(tt_rnk == "India") - 0.5, ymax = which(tt_rnk == "India") + 0.5, 
              fill = "grey", alpha = .02) +
    geom_pointrange(aes(x = M, xmin = L, xmax = U, y = Location, colour = Stage), 
                    position = position_dodge(width = .5)) +
    scale_x_continuous("Time arrival, month", 
                       breaks = c(seq(0, 2, 0.5), seq(3, 5, 1)),
                       labels = scales::number_format(scale = 12, accuracy = 1)) +
    scale_color_discrete("Time to", labels = tt_labels) +
    expand_limits(x = 0) +
    theme(legend.position = c(1, 0), legend.just = c(1.2, -0.3)) +
    guides(color = guide_legend(reverse=TRUE))
  
  
  
  ggsave(g_tte, filename = here::here("docs", "figs", "g_tte_" + tp) + ext, width = 7, height = 8)    
  
  
}
