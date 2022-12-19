library(tidyverse)



load(here::here("data", "dat_cas.rdata"))

locations <- unique(dat_noti$State)


## Setup
set.seed(1166)

n_samples <- 300

dir.create(here::here("docs", "pars"), showWarnings = F)

## Collect pars
for (loc in names(loc_maps)) {
  loc <- glue::as_glue(loc)
  
  
  if (loc == "Arunachal_Pradesh") {
    txi <- jsonlite::read_json(here::here("out", "sub_txi", "pars_India.json"))
    txi <- txi[sample(1:length(txi), n_samples, rep =T)]
  } else {
    txi <- jsonlite::read_json(here::here("out", "sub_txi", "pars_" + loc + ".json"))
    txi <- txi[sample(1:length(txi), n_samples, rep =T)]
  }

  
  txo <- jsonlite::read_json(here::here("out", "sub_tx", "pars_" + loc + ".json"))
  txo <- txo[sample(1:length(txo), n_samples, rep =T)]
  

  cs <- jsonlite::read_json(here::here("out", "sub_cs", "post_cs_s1_" + loc + ".json"))
  cs <- cs[sample(1:length(cs), n_samples, rep =T)]
  
  pars <- lapply(1:n_samples, function(i) c(cs[[i]], txo[[i]], txi[[i]]))
  
  jsonlite::write_json(pars, here::here("docs", "pars", "pars_s1_" + loc + ".json"), auto_unbox=T, digits = 10)
  
  
  cs <- jsonlite::read_json(here::here("out", "sub_cs", "post_cs_s2_" + loc + ".json"))
  cs <- cs[sample(1:length(cs), n_samples, rep =T)]
  
  pars <- lapply(1:n_samples, function(i) c(cs[[i]], txo[[i]], txi[[i]]))
  
  jsonlite::write_json(pars, here::here("docs", "pars", "pars_s2_" + loc + ".json"), auto_unbox=T, digits = 10)
  
  
  cs <- jsonlite::read_json(here::here("out", "sub_cs", "post_cs_s3_" + loc + ".json"))
  cs <- cs[sample(1:length(cs), n_samples, rep =T)]
  
  pars <- lapply(1:n_samples, function(i) c(cs[[i]], txo[[i]], txi[[i]]))
  
  jsonlite::write_json(pars, here::here("docs", "pars", "pars_s3_" + loc + ".json"), auto_unbox=T, digits = 10)
}

