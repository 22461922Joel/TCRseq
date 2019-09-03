clean_data(file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/PID summary"))

data <- read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/cleaned_CDR3s.csv"))

data$exp <- str_remove(data$exp, "mTCR ") %>%
   str_replace(" ", "-") %>%
   str_replace_all("627-RS", "627-RS-4") %>%
   str_replace_all("630-NR", "630-NR-6") %>%
   str_replace_all("633-RS", "633-RS-7") %>%
   str_replace_all("640-RS", "640-RS-8") %>%
   str_replace_all("-", "_") %>%
   str_remove("RENCA_")

unique(data$exp)

data$model <- "RENCA"

stri_sub(data$exp, 2, 1) <- "_"

data <- data %>%
  unite("exp", model, exp, sep = "_")

write_csv(data, file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/cleaned_CDR3s.csv"))

rm(data)

# setwd("D://data//experiments//RNAseq//data_individual")
# 
# for (i in 1:length(data_individual)) {
#   write.csv(data_individual[[i]], paste(unique(data_individual[[i]]$exp), ".csv", sep = ""))
# }



neg_RENCA <- PID_control("D://data//experiments//RNAseq//1239Shp15b_Joost_final")

PID <- PID_control("D://data//experiments//RNAseq//1239Shp15b_Joost_final//PID\ summary")