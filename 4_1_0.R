# sort_function("D:/data/experiments/4_1_0/sort_PDFs")

setwd("D:/data/experiments/4_1_0")

# sort <- read.csv("sort_data.csv") %>%
#   select(-X, -population, -exp, "exp" = "storage")
# 
# write.csv(sort, "sort_data.csv")
# 
# sort <- read.csv("sort_data.csv")

tumour_growth("D:/data/experiments/4_1_0/tumour_growth")

setwd("D:/data/experiments/4_1_0")

sort <- read.csv("sort_data.csv")
# 
# sort$storage <- sort$storage %>% 
#   str_remove("REAL_") %>%
#   str_replace_all("[^a-zA-Z0-9.]", "_") %>%
#   str_replace("LL", "L_L") %>%
#   str_replace("LR", "L_R") %>%
#   str_replace("TL", "T_L") %>%
#   str_replace("TR", "T_R") %>%
#   str_replace_all("__", "_")
# 
# write.csv(sort, "sort_data.csv")

# clean_data("D:/data/experiments/4_1_0/1239Shp14_Final MiXCR results/PID summary")
# 
# data <- read.csv("cleaned_CDR3s.csv")
# 
# data$exp <- data$exp %>%
#   str_replace("_", ".") %>%
#   str_replace("-", ".") %>%
#   str_replace("_PCR2_1_", "_1.") %>%
#   str_replace("-", "_") %>%
#   str_replace("-", ".") %>%
#   str_replace("LT", "T_L") %>%
#   str_replace("RT", "T_R") %>%
#   str_replace("L4", "T_L_CD4") %>%
#   str_replace("L8", "T_L_CD8") %>%
#   str_replace("R4", "T_R_CD4") %>%
#   str_replace("R8", "T_R_CD8") %>%
#   str_replace_all("-", "_")
# 
# write.csv(data, "cleaned_CDR3s.csv")

neg_4_1_0 <- PID_control("D://data//experiments//4_1_0//1239Shp14_Final MiXCR results")

data <- read.csv("cleaned_CDR3s.csv")

summary_TCRseq(data)

summary <- read.csv("summary_stats.csv")

jk41 <- data_aa(data) %>% factor_extractor()

morisita_network(jk41, "CD4")

morisita_network(jk41, "CD8")

jk41_overlaps <- all_intersects(jk41)

jk41_overlaps$same_mouse <- mice_matching(jk41_overlaps)

ggplot(jk41_overlaps, aes(factor(n_mice, levels = 2:12), n_clones, group = n_mice, colour = same_mouse)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(. ~ pop.x) +
  labs(y = "matching clones", x = "n mice") +
  theme(legend.position = "none")

jk41_pairs <- jk41_overlaps %>%
  filter(n_mice == 2) %>%
  separate(mice, 
           into = c(letters[1:2], rep(NA, times = (str_count(.$mice[1], "-") - 1))), 
           sep = "-") %>%
  mutate(same_mouse = stringdist(.$a, .$b))

ggplot(jk41_pairs, aes(as.character(same_mouse), n_clones)) +
  geom_jitter(width = 0.1) +
  facet_grid(. ~ pop.x)

all_full_joins_graph(jk41, "CD8")

all_full_joins_graph(jk41, "CD4")

hd_dual_tumour_abundance(jk41, "CD8")

#####
# PID sharing between each sample
#####

setwd("D://data//experiments//4_1_0//1239Shp14_Final MiXCR results")

setwd(paste(getwd(), "//PID summary", sep = ""))

names_exp <- dir()[str_detect(dir(), " miXCR clones read names with PIDs.csv")]

comp_list <- rep(list(1), times = length(names_exp))

comp_temp1 <- list()

for (j in 1:24) {
  
  comp_temp1[[j]] <- read_csv(names_exp[j]) %>%
    select(-`Well`, -`Clone count`, -`read name`, -`Clone Id`) %>%
    filter(!str_detect(CDR3, "_")) %>%
    filter(!str_detect(CDR3, "\\*")) %>%
    filter(!str_length(CDR3) > 20) %>%
    filter(!str_length(CDR3) < 8) %>%
    distinct()
    
}

str_length(comp_temp1[[1]]$PID) %>% unique()

map_list <- list()

for (i in 1:length(comp_temp1)) {
  map_list[[i]] <- comp_temp1[1:length(comp_temp1) != i] %>%
      reduce(anti_join, by = c("CDR3", "PID"), .init = comp_temp1[[i]])

  map_list[[i]] <- map_list[[i]] %>%
    group_by(CDR3, `Subject Id`) %>%
    summarise(PID.count = n_distinct(PID))
}

map_list <- map_list %>%
  bind_rows()

 map_list$`Subject Id` <- map_list$`Subject Id` %>%
   str_replace("_", ".") %>%
   str_replace("-", ".") %>%
   str_replace("_PCR2_1_", "_1.") %>%
   str_replace("PCR2_1_", "_1.") %>%
   str_replace("-", "_") %>%
   str_replace("-", ".") %>%
   str_replace("LT", "T_L") %>%
   str_replace("RT", "T_R") %>%
   str_replace("L4", "T_L_CD4") %>%
   str_replace("L8", "T_L_CD8") %>%
   str_replace("R4", "T_R_CD4") %>%
   str_replace("R8", "T_R_CD8") %>%
   str_replace("0T", "0_T") %>%
   str_replace("2T", "2_T") %>%
   str_replace("3T", "3_T") %>%
   str_replace_all("-", "_") %>%
   str_remove("mTCR ")

   
jk41_t <- jk41 %>%
  full_join(map_list, by = c("aaSeqCDR3" = "CDR3", "exp" = "Subject Id")) %>%
  filter(PID.count > 2) %>%
  filter(!is.na(PID.count)) %>%
  na.omit() %>%
  group_by(exp) %>%
  mutate(PID.fraction_sum_t = PID.count / sum(PID.count),
         PID.count_difference = PID.count_sum - PID.count,
         PID.count_difference = scale(PID.count_difference)) %>%
  ungroup() %>%
  mutate(total_transcripts = sum(PID.count_sum)) %>%
  group_by(aaSeqCDR3) %>%
  mutate(overlap_prob = sum(PID.count_sum) / (total_transcripts * 4**9))

top_discord <- jk41_t %>%
  filter(PID.count_difference > 2)

top_agreement <- jk41_t %>%
  filter(PID.count == PID.count_sum)

ggplot(jk41_t, aes(PID.count_difference)) +
  geom_histogram(binwidth = 0.1)
  

ggplot(jk41_t, aes(PID.fraction_sum, PID.fraction_sum_t)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0)

ggplot(jk41_t, aes(PID.count_sum, PID.count, colour = cut(PID.count_difference, c(-Inf, -1, 2, Inf)))) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0) +
  labs(colour = "Z score")

#####
