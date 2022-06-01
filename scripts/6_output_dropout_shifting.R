library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


ext <- glue::as_glue(".png")


dead_labels <- c(
  DurAsym = "Asym.", DurNA = "Sym., unaware", 
  DurNC = "Sym. not consulted", DurCS = "Sym. not diagnosed"
)

drop_labels <- c(
  Drop0 = "Asym.", Drop1 = "Sym., unaware", 
  Drop2 = "Sym. not consulted", Drop3 = "Sym. not diagnosed",
  Drop4 = "Pre Tx", Drop5 = "On Tx", "NA" = "Location"
)



for (cnr_year in 2019:2021) {
  folder <- glue::as_glue("cas_") + cnr_year
  load(here::here("out", folder, "cascades_shifting.rdata"))
  

  dir.create(here::here("results", folder), showWarnings = F)
  
  
  loc <- unique(cascades_shifting$Location)
  loc <- c("India", loc[loc != "India"])
  
  
  drop <- cascades_shifting %>% 
    select(Location, starts_with("Drop")) %>%
    pivot_longer(-Location, names_to = "Stage") %>% 
    group_by(Location, Stage) %>% 
    summarise(
      M = median(value),
      L = quantile(value, 0.025),
      U = quantile(value, 0.975)
    ) %>% 
    ungroup() %>% 
    mutate(Location = factor(Location, loc)) %>% 
    separate(Stage, c("Stage", "Reason"))
  
  
  drop_n <- cascades_shifting %>% 
    mutate(
      across(starts_with("Drop"), function(x) x * 1) # burden
    ) %>% 
    select(Location, starts_with("Drop")) %>%
    pivot_longer(-Location, names_to = "Stage") %>% 
    group_by(Location, Stage) %>% 
    summarise(M = median(value)) %>% 
    ungroup() %>% 
    mutate(Location = factor(Location, loc)) %>% 
    separate(Stage, c("Stage", "Reason"))
  
  
  dead <- drop %>% filter(Reason == "dead")
  dead_n <- drop_n %>% filter(Reason == "dead")
  
  ltfu <- drop %>% filter(Reason != "dead")
  ltfu_n <- drop_n %>% filter(Reason == "dead")
  
  tree_dead <- treemap::treemap(
    dead_n %>% filter(Location != "India"),
    index = c("Location", "Stage"),
    vSize = "M"
  )$tm %>% 
    mutate(x1 = x0 + w, y1 = y0 + h) %>% 
    mutate(x = (x0+x1)/2, y = (y0+y1)/2) %>% 
    mutate(primary_group = ifelse(is.na(Stage), 1.2, .2)) 
  
  g <- tree_dead %>% 
    filter(!is.na(Stage)) %>% 
    ggplot(aes(xmin = x0, ymin = y0, xmax = x1, ymax = y1)) + 
    geom_rect(aes(fill = Stage, size = primary_group),
              show.legend = T, color = "black", alpha = .3) +
    scale_size(range = range(tree_dead$primary_group)) +
    scale_fill_discrete("Stage", labels = drop_labels) +
    geom_rect(data = filter(tree_dead, is.na(Stage)), 
              aes(fill = Stage, size = primary_group),
              show.legend = F, color = "black", alpha = .1) +
    ggfittext::geom_fit_text(data = filter(tree_dead, is.na(Stage)), aes(label = Location), min.size = 1) +
    guides(size = guide_none()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() + 
    theme(legend.position = "bottom") +
    labs(title = "Deaths due to TB by stage",
         caption = sprintf("Year: %d", cnr_year))
  
  ggsave(g, filename = here::here("results", folder, "g_tree_dead_shifting" + ext), width = 8, height = 8)  
  
  
  g <- drop %>% 
    group_by(Location) %>% 
    arrange(Stage, Reason) %>% 
    mutate(
      Reason = factor(Reason, rev(c("sc", "ltfu", "dead"))),
      Stage = factor(Stage),
      y1 = 1 - cumsum(M),
      y0 = c(1, y1[-n()]),
      StageNum = as.numeric(Stage)
    ) %>% 
    ungroup() %>% 
    ggplot(aes(x = Stage)) +
    geom_rect(aes(xmin = StageNum + 0.5, xmax = StageNum - 0.5, ymin = y1, ymax = y0, fill = Reason)) +
    scale_x_discrete("Stage", labels = drop_labels) +
    scale_y_continuous("Dropout, %", labels = scales::percent, limits = c(0, 1)) +
    scale_fill_discrete(labels = c(dead = "Death", ltfu = "Lost to follow-up", sc = "Self-cured")) +
    facet_wrap(.~Location) +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = -50, hjust = 0), 
          plot.margin = margin(5, 30, 5, 5)) +
    labs(caption = sprintf("Year: %d", cnr_year))
  
  
  ggsave(g, filename = here::here("results", folder, "g_dropout_shifting" + ext), width = 8, height = 8)  
}


  


