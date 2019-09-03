dual_tumour_growth("D:/data/experiments/4_2_0/tumour_growth")

# sort_function("D:/data/experiments/4_2_0/sort_PDFs")

working_path <- "experiments/4_2_0"

# sort <- read.csv("sort_data.csv") %>%
#   select(-population, -exp, "exp" = "storage")
# 
# write.csv(sort, "sort_data.csv")

sort <- read.csv("sort_data.csv")
# 
# sort$storage[sort$storage == "JK4.2_2._2_L_R_CD4"] <- "JK4.2_2.2_L_R_CD4"
# 
# sort$storage[sort$storage == "JK4.2_2._2_L_R_CD8"] <- "JK4.2_2.2_L_R_CD8"
# 
# sort$storage[sort$storage == "JK4.2_2._3_L_L_CD4"] <- "JK4.2_2.3_L_L_CD4"
# 
# sort$storage[sort$storage == "JK4.2_2._3_L_L_CD4"] <- "JK4.2_2.3_L_L_CD4"
# 
# sort$storage[sort$storage == "JK4.2_2._3_L_L_CD8"] <- "JK4.2_2.3_L_L_CD8"
# 
# write.csv(sort, "D:/data/experiments/4_2_0/sort_data.csv", row.names = F)


clean_data("D:/data/experiments/4_2_0/Amended and updated  1239Shp17 MixCR analysis all samples/PID summary")

tumour_growth <- read.csv("tumour_growth.csv")

neg_4_2_0 <- PID_control("D://data//experiments//4_2_0//1239Shp17 MixCR analysis all samples")

data <- read_csv(file.path(getwd(), working_path, "Amended and updated  1239Shp17 MixCR analysis all samples/cleaned_CDR3s.csv"))

summary_TCRseq(data)

summary <- read_csv(file.path(getwd(), working_path, "summary_stats.csv"))

factors <- c("jkexp", "mouse", "tissue", "flank", "pop")

jk42 <- data_aa(data)

jk42$exp <- jk42$exp %>%
  str_replace("_8", "_CD8") %>%
  str_replace("_4", "_CD4")

morisita_network(jk42, "4")

morisita_network(jk42, "8")

jk42_overlaps <- all_intersects(jk42)

jk42_overlaps$same_mouse <- mice_matching(jk42_overlaps)

jk42_overlaps <- jk42_overlaps %>%
  mutate(n_mice = factor(n_mice, levels = 2:10))

ggplot(jk42_overlaps, aes(n_mice, n_clones, group = n_mice, colour = same_mouse)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(. ~ pop.x) +
  labs(y = "matching clones", x = "n mice") +
  theme(legend.position = "none")

jk42_pairs <- jk42_overlaps %>%
  filter(n_mice == 2) %>%
  separate(mice, 
           into = c(letters[1:2], rep(NA, times = (str_count(.$mice[1], "-") - 1))), 
           sep = "-") %>%
  mutate(same_mouse = stringdist(.$a, .$b))

ggplot(jk42_pairs, aes(as.character(same_mouse), n_clones)) +
  geom_jitter(width = 0.1) +
  facet_grid(. ~ pop.x)

all_full_joins_graph(jk42, "CD4")

all_full_joins_graph(jk42, "CD8")

jk42_8fj <- all_full_joins_df(jk42, "CD8")

hd_dual_tumour_abundance(jk42, "CD8")
