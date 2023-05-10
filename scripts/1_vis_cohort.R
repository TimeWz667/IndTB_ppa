library(tidyverse)


theme_set(theme_bw() + theme(text = element_text(family = "sans")))


states <- read_csv(here::here("data", "state_space.csv"))

state_colour <- states %>% select(Colour, StageLab) %>% distinct()
state_colour <- set_names(state_colour$Colour, state_colour$StageLab)


runs <- read_csv(here::here("docs", "run_cohort.csv")) %>% 
  group_by(Time) %>% 
  summarise(across(starts_with(c("S_", "E_")), mean)) %>% 
  pivot_longer(-Time, names_to = "Stage") %>%
  left_join(states %>% select(Stage, StageLab)) %>% 
  mutate(
    StageLab = factor(StageLab, rev(names(state_colour)))
  )


head(runs)


cohort <- runs %>% filter(Time != max(Time))
endpt <- runs %>% 
  filter(Time == max(Time)) %>% 
  filter(startsWith(Stage, "E_"))


cohort %>% filter(Time <= 1.5) %>% 
  ggplot() +
  geom_bar(aes(x = Time, y = value, fill = StageLab), width=0.1, stat = "identity") +
  geom_bar(data = endpt, aes(x = 1.7, y = value, fill = StageLab), width=0.1, stat = "identity") +
  scale_fill_manual("Stage", values = state_colour) +
  scale_y_continuous("Proportion, %", labels = scales::percent) +
  scale_x_continuous("Time since incidence, months", labels = scales::number_format(scale=12))


