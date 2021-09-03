#code to prepare Theis_META
Theis_META <- read.csv("data-raw/Theis_SraRunTable.csv", stringsAsFactors = FALSE)
usethis::use_data(Theis_META, overwrite = TRUE)
