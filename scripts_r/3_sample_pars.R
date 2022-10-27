library(tidyverse)



load(here::here("data", "dat_cas.rdata"))

locations <- unique(dat_noti$State)
loc_maps <- setNames(gsub(" ", "_", locations), locations)
loc_maps["Chhatisgarh"] <- "Chhattisgarh"

loc_maps



## Setup
set.seed(1166)

n_samples <- 300

dir.create(here::here("docs", "pars"), showWarnings = F)

## Collect pars
for (loc in names(loc_maps)) {
  loc <- glue::as_glue(loc)
  loc1 <- glue::as_glue(loc_maps[loc])
  
  
  if (loc == "Arunachal Pradesh") {
    txi <- jsonlite::read_json(here::here("out", "sub_txi", "pars_India.json"))
    txi <- txi[sample(1:length(txi), n_samples, rep =T)]
  } else {
    txi <- jsonlite::read_json(here::here("out", "sub_txi", "pars_" + loc + ".json"))
    txi <- txi[sample(1:length(txi), n_samples, rep =T)]
  }

  
  txo <- jsonlite::read_json(here::here("out", "sub_tx", "pars_" + loc + ".json"))
  txo <- txo[sample(1:length(txo), n_samples, rep =T)]
  
  cs <- jsonlite::read_json(here::here("out", "sub_cs", "pars_nods_" + loc + ".json"))
  cs <- cs[sample(1:length(cs), n_samples, rep =T)]
  
  
  pars <- lapply(1:n_samples, function(i) c(cs[[i]], txo[[i]], txi[[i]]))
  
  jsonlite::write_json(pars, here::here("docs", "pars", "pars_nods_" + loc1 + ".json"), auto_unbox=T, digits = 10)
  
  
  cs <- jsonlite::read_json(here::here("out", "sub_cs", "pars_onds_" + loc + ".json"))
  cs <- cs[sample(1:length(cs), n_samples, rep =T)]
  
  
  pars <- lapply(1:n_samples, function(i) c(cs[[i]], txo[[i]], txi[[i]]))
  
  jsonlite::write_json(pars, here::here("docs", "pars", "pars_onds_" + loc1 + ".json"), auto_unbox=T, digits = 10)
  
}





