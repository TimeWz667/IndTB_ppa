library(tidyverse)




loc <- "Delhi"



n_samples <- 1500


post_cs <- read_csv(here::here("out", "sub_cs", "post_cs_onds_" + glue::as_glue(loc) + ".csv"))


jsonlite::toJSON(post_cs)

post_txi <- jsonlite::read_json(here::here("out", "sub_txi", "pars_" + glue::as_glue(loc) + ".json"))

post_tx <- jsonlite::read_json(here::here("out", "sub_tx", "pars_" + glue::as_glue(loc) + ".json"))


post_tx[[1]]

