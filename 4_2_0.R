tumour_growth("D:/data/experiments/4_2_0/tumour_growth")

# sort_function("D:/data/experiments/4_2_0/sort_PDFs")

setwd("D:/data/experiments/4_2_0")

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

# 
# clean_data("D:/data/experiments/4_2_0/1239Shp17 MixCR analysis all samples/PID summary")

data <- read.csv("cleaned_CDR3s.csv")

summary_TCRseq(data)

jk42 <- data_aa(data) %>% factor_extractor()

morisita_network(jk42, "CD4")

morisita_network(jk42, "CD8")

jk42_overlaps <- all_intersects(jk42)

jk42_pairs <- jk42_overlaps %>%
  filter(n_mice == 2) %>%
  separate(mice, 
           into = c(letters[1:2], rep(NA, times = (str_count(.$mice[1], "-") - 1))), 
           sep = "-") %>%
  mutate(same_mouse = stringdist(.$a, .$b))

ggplot(jk42_pairs, aes(as.character(same_mouse), n_clones)) +
  geom_jitter(width = 0.1) +
  facet_grid(. ~ pop.x)

all_full_joins(jk42, "CD4")

all_full_joins(jk42, "CD8")
