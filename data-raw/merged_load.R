merged <- read.csv("data-raw/merged.csv", stringsAsFactors = FALSE)
#need to edit merged a bit more before it is ready to be used in examples.
usethis::use_data(merged, overwrite = TRUE)
