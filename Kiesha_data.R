clean_data("D:/data/experiments/Kiesha_data/1239Shp17 MixCR analysis all samples/PID summary")

setwd("D:/data/experiments/Kiesha_data")

data <- read.csv("cleaned_CDR3s.csv", stringsAsFactors = F) %>%
  filter(str_detect(exp, "KR"))

data$exp[str_detect(data$exp, "-11")] <- "KR002_2.3_6"

data$exp[str_detect(data$exp, "-10")] <- "KR002_2.0_6"

data$exp[str_detect(data$exp, "-8")] <- "KR002_2.3_4"

data$exp[str_detect(data$exp, "-7")] <- "KR002_2.0_4"

data$exp[str_detect(data$exp, "-2")] <- "KR002_2.3_0"

data$exp[str_detect(data$exp, "-1")] <- "KR002_2.0_0"

neg_KR <- PID_control("D:/data/experiments/Kiesha_data/1239Shp17 MixCR analysis all samples")

factor_extractor <- function(data) {
  data %>% separate(exp, into = c("KRexp", 
                                  "mouse", 
                                  "timepoint"), sep = "_", remove = F)
}
  
data <- factor_extractor(data = data)

summary_TCRseq(data)

summary <- read.csv("summary_stats.csv") %>% factor_extractor()

ggplot(summary, aes(timepoint, diversity, colour = mouse, group = mouse)) +
  geom_point() +
  geom_line()

ggplot(summary, aes(evenness, richness, group = mouse, colour = mouse)) +
  geom_point() +
  facet_grid(. ~ timepoint)

ggplot(summary, aes(timepoint, richness, colour = mouse, group = mouse)) +
  geom_point() +
  geom_line()

aa_data <- data_aa(data)

mouse_CDR3_network(aa_list[[1]])
