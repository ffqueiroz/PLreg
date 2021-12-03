## code to prepare `PeruVotes` dataset goes here

load(file="data-raw/Peru-realdata.RData")

PeruVotes <- Peru.voters
colnames(PeruVotes) <- c("votes", "HDI")
usethis::use_data(PeruVotes, overwrite = TRUE)
