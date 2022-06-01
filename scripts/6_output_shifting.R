library(tidyverse)
library(odin)


theme_set(theme_bw() +
          theme(legend.position = "bottom", 
                axis.text.y = element_blank(), axis.ticks.y = element_blank()))

source(here::here("R", "calc_shifting.R"))



for (cnr_year in 2019:2021) {
  folder <- glue::as_glue("cas_") + cnr_year
  load(here::here("out", folder, "shifting.rdata"))
  
  folder <- glue::as_glue("shift_") + cnr_year
  dir.create(here::here("results", folder), showWarnings = F)
  
  locs <- names(sim_shifting)
  
  
  refs <- lapply(locs, function(loc) {
    tr_mat <- sim_shifting[[loc]]$Sim$tr_mat
    
    ref <- extract_referrals(tr_mat, loc)
    ref <- vis_referrals(ref$stocks, ref$flows)
    ref$Location <- loc
    ref
  })
  
  
  stocks <- bind_rows(lapply(refs, function(ref) {
    ref$stocks %>% mutate(Location = ref$Location)
  }))
  
  bands <- bind_rows(lapply(refs, function(ref) {
    ref$bands %>% mutate(Location = ref$Location)
  }))
  
  labx <- stocks %>% 
    mutate(xs = (x0 + x1) / 2) %>% 
    select(xs, Stage) %>% distinct()
  
  g <- ggplot(stocks %>% filter(Location != "India")) +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = Sector), 
              colour = "grey7") +
    geom_ribbon(data=bands %>% filter(Location != "India"), aes(x=xs, ymin=ys.lower, ymax=ys.upper, group = Key), 
                colour = "grey7", fill="darkgreen", alpha=0.3) +
    scale_x_continuous("Stage", breaks = labx$xs, 
                       labels = c("Initial visit", "Second visit", "Diagnosis")) +
    scale_fill_discrete("Sector", 
                        labels = c(pri = "Private", eng = "Engaged Private", pub = "Public")) +
    facet_wrap(.~Location)
  
  ggsave(g, filename = here::here("results", folder, "g_shift.png"), width = 12, height = 12)
  
  for(i in 1:length(refs)) {
    ref <- refs[[i]]
    loc <- glue::as_glue(ref$Location)
    g <- ref$g +
      labs(caption = loc)
    
    
    ggsave(g, filename = here::here("results", folder, "g_shift_" + loc + ".png"), width = 5, height = 4)
  }
}



