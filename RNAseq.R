clean_data("D://data//experiments//RNAseq//1239Shp15b_Joost_final//PID\ summary")

setwd("D://data//experiments//RNAseq")

data <- read.csv("cleaned_CDR3s.csv", stringsAsFactors = F) %>%
  select(-X)

data$exp <- str_remove(data$exp, "mTCR ") %>%
  str_replace(" ", "-") %>%
  str_replace_all("627-RS", "627-RS-4") %>%
  str_replace_all("630-NR", "630-NR-6") %>%
  str_replace_all("633-RS", "633-RS-7") %>%
  str_replace_all("640-RS", "640-RS-8")

factor_separator <- function(data) {
  data %>% separate(exp, into = c("tpmouse", "response", "primer"), sep = "-", remove = F) %>%
    separate(tpmouse, into = c("timepoint", "mouse"), sep = 1, remove = F)
}

data <- data %>%
  factor_separator()
# summary plots
#####

summary_TCRseq(data = data)

summary <- read.csv("summary_stats.csv", stringsAsFactors = F) %>%
  select(-X) %>%
  factor_separator()

summary$timepoint <- if_else(summary$timepoint == "0", "Pre CPB",
                             if_else(summary$timepoint == "2", "day 2",
                                     if_else(summary$timepoint == "4", "day 4", "day 6")))

summary$timepoint <- factor(summary$timepoint, levels = c("Pre CPB", "day 2", "day 4", "day 6", ordered = T))

summary$response <- factor(summary$response, levels = c("NR", "RS"), ordered = T)

response_fill_scale <- scale_fill_manual(values = c("NR" = "firebrick1",
                               "RS" = "deepskyblue1"))


response_colour_scale <- scale_colour_manual(values = c("NR" = "firebrick1",
                                                    "RS" = "deepskyblue1"))

ggplot(summary, aes(timepoint, diversity, fill = response)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(group = response), label = "p.signif") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  response_fill_scale

ggplot(summary, aes(evenness, richness, colour = response)) +
  geom_point() +
  scale_y_log10(breaks = seq(from = 1000, to = max(summary$richness), by = 1000)) +
  scale_x_log10(breaks = seq(from = 0, to = max(summary$evenness), by = 25000)) +
  facet_grid(. ~ timepoint) +
  labs(y = "unique clones", x = "total clones") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  response_colour_scale


#####