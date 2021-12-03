## code to prepare `Firm` dataset goes here

load(file="data-raw/Firm.RData")
Firm
colnames(Firm) <- c("firmcost", "assume", "cap", "sizelog", "indcost", "central",
                    "soph")
usethis::use_data(Firm, overwrite = TRUE)
