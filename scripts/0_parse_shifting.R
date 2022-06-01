library(tidyverse)


txt = read_file(here::here("data", "shifting", "shifting.txt"))


entries <- strsplit(txt, "\r\n")[[1]]

pid <- 0

shifting <- bind_rows(lapply(entries, function(ent) {
  ent <- strsplit(ent, " ")[[1]]
  
  pid <<- pid + 1
  len <- length(ent) - 2
  
  path <- ent[1:len]
  path <- ifelse(path == "Public", "pub", "pri")
  path[len] <- "pub_det"
  path <- c("UC", path)
  
  tibble(
    PID = pid,
    src = paste(path, collapse = ":"),
    From = path[-length(path)],
    To = path[-1],
    N = as.numeric(ent[len + 1]),
    length = len
  )
}))

save(shifting, file = here::here("data", "shifting", "shifting.rdata"))
write_csv(shifting, here::here("data", "shifting", "shifting.csv"))


shifting %>% 
  select(PID, N, length) %>% 
  distinct() %>% 
  summarise(
    med = sum(N * length) / sum(N)
  )




shifting %>% 
  group_by(PID) %>% 
  mutate(
    without_pri = 1 - any(To == "pri")
  ) %>% 
  ungroup() %>% 
  select(PID, N, without_pri, length) %>% 
  distinct() %>% 
  summarise(
    med = sum(N * length * without_pri) / sum(N * without_pri)
  )
