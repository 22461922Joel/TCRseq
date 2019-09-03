path <- file.path(getwd(), "experiments/RNAseq/1239Shp19 Final MixCR analysis/PID analysis")

clean_data(path)

data <- read_csv(file.path(path, "../cleaned_CDR3s.csv")) %>%
  filter(!str_detect(exp, "neg"))

data$exp <- data$exp %>%
  str_replace_all("-", "_")

stri_sub(data$exp, 6, 5) <- "_"

write_csv(data, file.path(path, "../cleaned_CDR3s.csv"))
