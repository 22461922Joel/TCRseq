# sort_function("D:/data/experiments/4_1_0/sort_PDFs")

setwd("D:/data/experiments/4_1_0")

# sort <- read.csv("sort_data.csv") %>%
#   select(-X, -population, -exp, "exp" = "storage")
# 
# write.csv(sort, "sort_data.csv")

sort <- read.csv("sort_data.csv")

tumour_growth("D:/data/experiments/4_1_0/tumour_growth")

setwd("D:/data/experiments/4_1_0")

sort <- read.csv("sort_data.csv")

sort$storage <- sort$storage %>% 
  str_remove("REAL_") %>%
  str_replace_all("[^a-zA-Z0-9.]", "_") %>%
  str_replace("LL", "L_L") %>%
  str_replace("LR", "L_R") %>%
  str_replace("TL", "T_L") %>%
  str_replace("TR", "T_R") %>%
  str_replace_all("__", "_")

write.csv(sort, "sort_data.csv")

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

data <- read.csv("cleaned_CDR3s.csv")

summary_TCRseq(data)

summary <- read.csv("summary_stats.csv")

jk41 <- data_aa(data) %>% factor_extractor()

morisita_network(jk41, "CD4")

morisita_network(jk41, "CD8")

jk41_overlaps <- all_intersects(jk41)

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
