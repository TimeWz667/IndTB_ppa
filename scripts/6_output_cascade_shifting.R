library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


labels <- c(
  Cas0 = "Sympton Developed", Cas1 = "Symptom Noticed", Cas2 = "Consulted", 
  Cas3 = "Diagnosed", Cas4 = "Tx initialised", Cas5 = "Tx Successful"
)


ext <- glue::as_glue(".png")


for (cnr_year in 2019:2021) {
  folder <- glue::as_glue("cas_") + cnr_year
  load(here::here("out", folder, "cascades_shifting.rdata"))
  
  dir.create(here::here("results", folder), showWarnings = F)
  
  loc <- unique(cascades_shifting$Location)
  loc <- c("India", loc[loc != "India"])
  
  
  cas <- cascades_shifting %>% 
    select(Location, starts_with("Cas")) %>% 
    pivot_longer(-Location, names_to = "Stage") %>% 
    group_by(Location, Stage) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.025),
      U = quantile(value, 0.975)
    ) %>% 
    ungroup() %>% 
    mutate(Location = factor(Location, loc))
  
  
  cas %>% 
    mutate(
      MCI = sprintf("%s (%s-%s)", scales::percent(M), scales::percent(L), scales::percent(U))
    ) %>% 
    select(-M, -L, -U) %>% 
    pivot_wider(names_from = Stage, values_from = MCI) %>% 
    arrange(Location) %>% 
    write_csv(here::here("results", folder, "CascadeShifting.csv"))
  
  
  g <- cas %>% 
    filter(Stage %in% names(labels)) %>% 
    ggplot() +
    geom_bar(aes(x = Stage, y = M, fill = Stage), stat = "identity") +
    geom_linerange(aes(x = Stage, ymin = L, ymax = U)) +
    scale_x_discrete(labels = labels) + 
    scale_y_continuous("Arrival, %", labels = scales::percent, limits = 0:1) +
    guides(fill = guide_none(), colour = guide_none()) +
    facet_wrap(.~Location) +
    theme(axis.text.x = element_text(angle = -50, hjust = 0), plot.margin = margin(5, 30, 5, 5)) +
    labs(caption = sprintf("Year: %d", cnr_year))
  
  
  g <- cas %>% 
    filter(Location == "India") %>% 
    filter(Stage %in% names(labels)) %>% 
    ggplot() +
    geom_bar(aes(x = Stage, y = M, fill = Stage), stat = "identity") +
    geom_linerange(aes(x = Stage, ymin = L, ymax = U)) +
    scale_x_discrete(labels = labels) + 
    scale_y_continuous("Arrival, %", labels = scales::percent, limits = 0:1) +
    guides(fill = guide_none(), colour = guide_none()) +
    facet_wrap(.~Location) +
    theme(axis.text.x = element_text(angle = -50, hjust = 0), plot.margin = margin(5, 30, 5, 5)) +
    labs(caption = sprintf("Year: %d", cnr_year))
  
  ggsave(g, filename = here::here("results", folder, "g_cascade_shifting_national" + ext), width = 4, height = 6)  
  
  ggsave(g, filename = here::here("results", folder, "g_cascade_shifting" + ext), width = 8, height = 8)  
  
}

