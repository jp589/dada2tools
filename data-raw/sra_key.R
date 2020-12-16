## code to prepare `sra_key` dataset goes here
sra_key <- read.csv("data-raw/sample_key.csv", stringsAsFactors = FALSE)
sra_key <- sra_key[,2:5]
usethis::use_data(sra_key, overwrite = TRUE)
