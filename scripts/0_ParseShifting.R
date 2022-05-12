library(tidyverse)


txt = read_file(here::here("data", "shifting"))


entries <- strsplit(txt, "\r\n")[[1]]

pid <- 0

shifting <- bind_rows(lapply(entries, function(ent) {
  ent <- strsplit(ent, " ")[[1]]
  
  pid <<- pid + 1
  len <- length(ent) - 2
  
  path <- ent[1:len]
  path <- c("UC", path, "det")
  
  tibble(
    PID = pid,
    From = path[-length(path)],
    To = path[-1],
    N = as.numeric(ent[len + 1]),
    length = len
  )
})) %>% 
  mutate(
    From = case_when(
      From == "Public" ~ "pub",
      From == "Private" ~ "pri",
      T ~ From
    ),
    To = case_when(
      To == "Public" ~ "pub",
      To == "Private" ~ "pri",
      T ~ To
    )
  )

save(shifting, file = here::here("data", "shifting.rdata"))



shifting %>% 
  select(PID, N, length) %>% 
  distinct() %>% 
  group_by(length) %>% 
  summarise(
    M = mean(N)
  ) %>% 
  ungroup() %>% 
  summarise(
    med = sum(M * length) / sum(M)
  )


