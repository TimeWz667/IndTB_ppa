library(tidyverse)
library(tidybayes)

theme_set(theme_bw())



locations <- read_csv(here::here("data", "locations.csv"))

locations



inc <- bind_rows(lapply(locations$Location, function(loc) {
  read_csv(here::here("data", loc, "Post.csv")) %>% mutate(Location = loc)
}))



ggplot(inc) +
  stat_halfeye(aes(x = r_csi, y = Location))

