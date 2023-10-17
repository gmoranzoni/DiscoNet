# wrangle inweb database
# Load necessary library
library(dplyr)
library(stringr)
# load the 2016 version of the inweb database
inweb <- read.delim("inst/extdata/inweb.txt", header = F, stringsAsFactors = F)
col <- c(1, 2, 15) # select the relavnt columns - description of the columns is in the README_confidence_score_inweb2016.txt file
inweb <- inweb[, col]


inweb <- inweb %>%
  # Rename columns
  rename(
    node = V1,
    interactor = V2,
    confidence_score = V15
  ) %>%
  # Remove 'uniprotkb:' from 'node' and 'interactor'
  mutate(
    node = str_replace(node, "uniprotkb:", ""),
    interactor = str_replace(interactor, "uniprotkb:", "")
  ) %>%
  # Keep the score before the '|' separator in 'confidence_score'
  mutate(
    confidence_score = as.numeric(str_extract(confidence_score, "^[^|]+"))
  )


usethis::use_data(inweb, overwrite = TRUE)
