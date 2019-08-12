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
# 
# clean_data("D:/data/experiments/4_1_0/1239Shp14_Final MiXCR results/PID summary")
# 
# data <- read.csv("cleaned_CDR3s.csv")
# 
# data$exp <- data$exp %>%
#    str_replace("_", ".") %>%
#    str_replace("-", ".") %>%
#    str_replace("_PCR2_1_", "_1.") %>%
#    str_replace("-", "_") %>%
#    str_replace("-", ".") %>%
#    str_replace("LT", "T_L") %>%
#    str_replace("RT", "T_R") %>%
#    str_replace("L4", "T_L_CD4") %>%
#    str_replace("L8", "T_L_CD8") %>%
#    str_replace("R4", "T_R_CD4") %>%
#    str_replace("R8", "T_R_CD8") %>%
#    str_replace_all("-", "_")
# 
# write.csv(data, "cleaned_CDR3s.csv", row.names = F)

neg_4_1_0 <- PID_control("D://data//experiments//4_1_0//1239Shp14_Final MiXCR results")

data <- read.csv("cleaned_CDR3s.csv", stringsAsFactors = F)

data$exp %>% unique()

factors <- c("jkexp", "mouse", "tissue", "flank", "population")

summary_TCRseq(data)

summary <- read.csv("summary_stats.csv")

jk41 <- data_aa(data)

heatmap_dendrogram(jk41, "flank", "L")

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
  
  comp_temp1[[j]]$`Subject Id` <- comp_temp1[[j]]$`Subject Id` %>%
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
    
}

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
 
 map_list_x <- list()
 
for (i in 1:length(comp_temp1)) {
   map_list_x[[i]] <- comp_temp1[1:length(comp_temp1) != i] %>%
     reduce(left_join, by = c("CDR3", "PID"), .init = comp_temp1[[i]])
   # 
   # map_list_x[[i]] <- map_list_x[[i]] %>%
   #   group_by(CDR3, `Subject Id`) %>%
   #   summarise(PID.count = n_distinct(PID))
}
 
map_list_x <- map_list_x %>%
   bind_rows()%>%
   unite("mice", starts_with("Subject"), sep = "-", remove = T) %>%
   mutate(n_mice = str_count(mice, "_[\\d].[\\d]_"))
 
unique(map_list_x$mice)

map_temp <- map_list_x %>%
  select(-mice) %>%
  distinct()
   
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
  group_by(aaSeqCDR3) %>%
  mutate(n_identical = sum(PID.count_sum),
         n_identical_corrected = sum(PID.count),
         overlap_prob = gam_prob(n_identical),
         overlap_prob_corrected = gam_prob(n_identical_corrected))

gam_prob <- function(x) {
  1 - exp(lgamma(4**9 + 1) - lgamma(4**9 - x + 1) - x * log(4**9))
}

poisson_prob <- function(x) {
  1 - exp(-CombN(x = x, m = 100)/(4**9)**2)
}

1-exp(-CombN(1:30, m = 3)/365**2)

poisson_prob(map_temp$PID[map_temp$CDR3 == "CASGETGTNERLFF"])

CombN(map_temp$PID[map_temp$CDR3 == "CASGETGTNERLFF"], m = 2)

jk41_temp <- jk41_t %>%
  select(aaSeqCDR3, overlap_prob) %>%
  distinct()

jk41_gamma <- jk41_temp %>%
  full_join(map_temp, by = c("aaSeqCDR3" = "CDR3")) %>%
  group_by(aaSeqCDR3, overlap_prob) %>%
  summarise(PID.count = n_distinct(PID))

gamma_df <- data.frame(PID.count = seq(from = min(jk41_gamma$PID.count), to = 4500, by = 100)) %>%
  mutate(overlap_prob = gam_prob(PID.count),
         aaSeqCDR3 = "model") %>%
  bind_rows(jk41_gamma)

gamma_df
  

ggplot(gamma_df, aes(y = overlap_prob, x = PID.count, colour = aaSeqCDR3)) +
  geom_point() +
  scale_x_log10(limits = c(1, 2000)) +
  theme(legend.position = "none")

ggplot(jk41_t, aes(n_identical_corrected, overlap_prob_corrected)) +
  geom_point() +
  scale_x_log10()

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
  labs(colour = "Z score", y = "PID counts with no identical PID/RNA combinations", x = "PID counts")

#####
