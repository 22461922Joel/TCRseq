library(tidyverse)

path <- file.path("D:/data/experiments/RNAseq/1239Shp15b_Joost_final/PID summary/")

files <- list.files(path = path, pattern = "*summary\ with\ PID.csv")

read_count_table <- function(path, filename, mat) {
  tmp <- read_csv(paste0(path, "/", filename), col_names = T) %>%
    dplyr::group_by(`Subject Id`, `CDR3`) %>%
    summarise(count = n_distinct(PID))
  mat <- bind_rows(mat, tmp)
  rm(tmp)
  return(mat)
}

mat <- tibble(`Subject Id` = character(), count = integer(), CDR3 = character())

for (file in files) {
  mat = read_count_table(path, file, mat)
}

library(data.table)

colnames(mat) = c("ID", "Count", "Sequence")

f = dcast(data = mat, formula = Sequence~ID, fun.aggregate = sum, value.var = "Count")

f_filtered = dplyr::filter(f, !grepl("\\*|_", Sequence))

rm(f)

write_csv(f_filtered, "D:/data/experiments/RNAseq/1239Shp15b_Joost_final/PID_matrix.csv", na = "NA", append = F, quote_escape = "double")
