library(tidyverse)
library(tidybayes)

theme_set(theme_bw() + theme(text = element_text(family = "sans")))


ext <- glue::as_glue(".png")

labels <- c(
  Cas0 = "Sympton Developed", Cas1 = "Symptom Noticed", Cas2 = "Consulted", 
  Cas3 = "Diagnosed", Cas4 = "Tx initialised", Cas5 = "Tx Successful"
)

tt_labels <- c(
  TTSym = "Sympton Developed", TTAware = "Symptom Noticed", 
  TTCare = "Care sought", TTDet = "Diagnosed"
)

dur_labels <- c(
  DurAsym = "Asym.", DurNA = "Sym., unaware", 
  DurNC = "Sym. not consulted", DurCS = "Sym. not diagnosed"
)



for (scenario in c("shared_pr_asym", "shared_r_onset")) {
  for (cnr_year in 2019:2021) {
    folder <- glue::as_glue("cas_") + cnr_year
    file_cascades <- glue::as_glue("cascades_shifting_") + scenario + ".rdata"
    load(here::here("out", folder, file_cascades))
    
    folder <- folder + "_" + scenario
    dir.create(here::here("results", folder), showWarnings = F)
    
    
    loc <- unique(cascades_shifting$Location)
    loc <- c("India", loc[loc != "India"])
    
    
    dur <- cascades_shifting %>% 
      select(Location, starts_with("Dur")) %>%
      pivot_longer(-Location, names_to = "Stage") %>% 
      group_by(Location, Stage) %>% 
      summarise(
        M = median(value),
        L = quantile(value, 0.025),
        U = quantile(value, 0.975)
      ) %>% 
      ungroup() %>% 
      mutate(Location = factor(Location, loc))
    
    
    tte <- cascades_shifting %>% 
      select(Location, starts_with("TT")) %>%
      pivot_longer(-Location, names_to = "Stage") %>% 
      group_by(Location, Stage) %>% 
      summarise(
        M = median(value),
        L = quantile(value, 0.025),
        U = quantile(value, 0.975)
      ) %>% 
      ungroup() %>% 
      mutate(Location = factor(Location, loc))
    
    
    delay <- cascades_shifting %>% 
      select(Location, starts_with("Delay")) %>%
      pivot_longer(-Location, names_to = "Stage") %>% 
      group_by(Location, Stage) %>% 
      summarise(
        M = median(value),
        L = quantile(value, 0.025),
        U = quantile(value, 0.975)
      ) %>% 
      ungroup() %>% 
      mutate(Location = factor(Location, loc))
      
    
    delay %>% 
      mutate(
        MCI = sprintf("%.1f (%.1f-%.1f)", M * 12, L * 12, U * 12),
        Stage = paste0(Stage, "_mo")
      ) %>% 
      select(-M, -L, -U) %>% 
      pivot_wider(names_from = Stage, values_from = MCI) %>% 
      arrange(Location) %>% 
      write_csv(here::here("results", folder, "DelayShifting.csv"))
    
    
    tt_rnk <- tte %>% 
      filter(Stage == "TTDet") %>% 
      arrange(M) %>% 
      pull(Location) %>% 
      as.character()
    
    g <- tte %>% 
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
                         breaks = seq(0, 2, 0.5),
                         labels = scales::number_format(scale = 12, accuracy = 1)) +
      scale_color_discrete("Time to", labels = tt_labels) +
      expand_limits(x = 0) +
      theme(legend.position = c(1, 0), legend.just = c(1.2, -0.3)) +
      guides(color = guide_legend(reverse=TRUE)) +
      labs(caption = sprintf("Year: %d, Assumption: %s", cnr_year, scenario))
     
    ggsave(g, filename = here::here("results", folder, "g_tte_shifting") + ext, width = 7, height = 8)  
    
    
    g <- dur %>% 
      mutate(
        Stage = factor(Stage, rev(names(dur_labels))),
        Location = factor(Location, rev(loc))
      ) %>% 
      ggplot() +
      geom_bar(aes(x = M, y = Location, fill = Stage), stat = "identity", position = "fill") +
      scale_fill_discrete("Stage", labels = dur_labels) +
      scale_x_continuous("Precentage, %", labels = scales::percent) +
      guides(fill = guide_legend(reverse = TRUE)) + 
      theme(legend.position = "bottom") +
      labs(caption = sprintf("Year: %d, Assumption: %s", cnr_year, scenario))
    
    
    ggsave(g, filename = here::here("results", folder, "g_spent_shifting" + ext), width = 7, height = 8)  
    
  }
}


