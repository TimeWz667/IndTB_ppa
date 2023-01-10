library(tidyverse)



ds <- dir(here::here("out", "sub_cs"))
ds <- ds[startsWith(ds, "post_cs_s")]
ds <- ds[endsWith(ds, ".rdata")]


res <- bind_rows(lapply(ds, function(d) {
  load(file = here::here("out", "sub_cs", d))
  
  x <- rstan::summary(post)$summary
  as_tibble(x) %>% 
    mutate(name = rownames(x), file=d) %>% 
    relocate(name, file)
}))


res %>% 
  extract(file, c("Scenario", "State"), "post_cs_(s1|s2|s3)_(\\S+).rdata") %>% 
  filter(
    startsWith(name, "nr_") | startsWith(name, "prv") | startsWith(name, "pr_") | startsWith(name, "tbps")
  ) %>% 
  write_csv(here::here("docs", "tabs", "to_fit.csv"))

