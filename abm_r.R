library(tidyverse)



n_pop <- 1000


pars <- list(
  beta = 1.5,
  gamma = 0.5
)

summ <- function(df, ti) {
  list(
    Time = ti,
    S = mean(df$State == "Sus"),
    I = mean(df$State == "Inf"),
    R = mean(df$State == "Rec")
  )
}


prog <- function(df, pars, ti, dt) {
  foi <- pars$beta * mean(df$State == "Inf")
  
  df %>% 
    mutate(
      p_infection = pexp(1, rate = foi * dt),
      p_rec = pexp(1, rate = pars$gamma * dt),
      State = case_when(
        State == "Sus" ~ ifelse(runif(n()) < p_infection, "Inf", "Sus"),
        State == "Inf" ~ ifelse(runif(n()) < p_rec, "Rec", "Inf"),
        T ~ State
      )
    ) %>% 
    select(-p_infection, -p_rec)
}

pop <- tibble(ID = 1:n_pop, 
              State = rep(c("Inf", "Sus"), c(100, n_pop - 100)))


ti <- 0
ss <- list(summ(pop, ti))


for (ti in seq(0.2, 10, 0.2)) {
  pop <- prog(pop, pars, ti, 0.2)
  
  ss[[length(ss) + 1]] <- summ(pop, ti)
}

ss <- bind_rows(ss)

ss %>% 
  ggplot() +
  geom_line(aes(x = Time, y = S, colour = "S")) +
  geom_line(aes(x = Time, y = I, colour = "I")) +
  geom_line(aes(x = Time, y = R, colour = "R"))


