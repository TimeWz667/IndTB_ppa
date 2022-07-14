library(tidyverse)
library(tidybayes)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))


tars <- read_csv(here::here("docs", "targets.csv"))
region <- read_csv(here::here("data", "StateMap.csv")) %>% select(Location = State, Region) %>% 
  mutate(Location = ifelse(Location == "Chhatisgarh", "Chhattisgarh", Location))


inc <- tars %>%
  left_join(region) 



ord <- inc %>% 
  group_by(Location) %>% 
  summarise(ir = mean(IncR)) %>% 
  arrange(ir) %>% pull(Location)



g_inc <- inc %>% 
  mutate(
    Location = factor(Location, ord),
    Region = factor(Region, region %>% pull(Region) %>% unique())
  ) %>% 
  ggplot() +
  geom_rect(xmin = -Inf, xmax = Inf, 
            ymin = which(ord == "India") - 0.5, ymax = which(ord == "India") + 0.5, 
            fill = "grey", alpha = .02) +
  stat_halfeye(aes(x = IncR, y = Location, fill = Region)) +
  scale_x_continuous("Incidence estimates, per 100 000", 
                     breaks = c(0, 100, 200, 300, 500, 1000, 1500) * 1e-5,
                     labels = scales::number_format(scale = 1e5), limits = c(0, 0.015))


g_inc

ggsave(g_inc, filename = here::here("docs", "inc.png"), height = 10, width = 9)
  
