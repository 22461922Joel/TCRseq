library(vegan)

library(MASS)

library(tidyverse)

counts = read_csv("D:/data/experiments/RENCA_HA_AB1_HA_exp12_exp13/combined_data/count_matrix.csv")

colnames(counts) <- colnames(counts) %>%
  str_remove("mTCR") %>%
  str_trim(side = "both") %>%
  str_replace(" ", "_") %>%
  str_replace_all("-", "_") %>%
  str_replace("R_[1-9]?[0-9]", "R") %>%
  str_replace("S_[1-9]?[0-9]", "S")

counts <- column_to_rownames(counts, var = "Sequence")

counts <- t(counts)

counts_dis <- vegdist(counts)

counts_mds0 <- isoMDS(counts_dis)

factors <- c("jcexp", "tissue", "mouse", "population")

counts_mds0_df <- counts_mds0 %>%
  as.data.frame() %>%
  rownames_to_column("exp")

counts_mds0_df$exp <- counts_mds0_df$exp %>%
  str_replace("_1_PCR2", "\\.1") %>%
  str_replace("_1PCR2", "\\.1") %>%
  str_remove("mTCR ") %>%
  str_replace("4-1", "4\\.1")
  str_remove("_12") %>%
  str_remove("T") %>%
  str_replace("11__", "11_T_") %>%
  str_replace("15.2_", "15.2_T_") %>%
  str_replace("__", "_")
  

counts_mds0_df$exp[str_detect(counts_mds0_df$exp, "^[0-9][0-9][0-9]")] <- str_replace(counts_mds0_df$exp[str_detect(counts_mds0_df$exp, "^[0-9][0-9][0-9]")], "^", "RENCA_")

counts_mds0_df <- counts_mds0_df %>%
  factor_extractor() %>%
  filter(tissue == "T", population == "CD8")

ggplot(counts_mds0_df, aes(points.1, points.2, colour = jcexp)) +
  geom_point(size = 3) +
  labs(y = "dim2", x = "dim1")

stressplot(counts_mds0, counts_dis)

ordiplot(counts_mds0, type = "t")

counts_mds <- metaMDS(counts)

counts_mds_df <- counts_mds$points %>%
  as.data.frame() %>%
  rownames_to_column("exp")

counts_mds_df$exp[str_detect(counts_mds_df$exp, "^[0-9][0-9][0-9]")] <- str_replace(counts_mds_df$exp[str_detect(counts_mds_df$exp, "^[0-9][0-9][0-9]")], "^", "RENCA_")

counts_mds_df <- counts_mds_df %>%
  factor_extractor()

ggplot(counts_mds_df, aes(MDS1, MDS2, colour = exp)) +
  geom_point()

counts_pca <- rda(counts)

raremax <- min(rowSums(counts))

rarecurve(counts[1:2,], sample = raremax, step = 1000)

counts_list <- counts %>%
  as.list()

typeof(counts_list[1])

sum(counts_list[[1]])

size_list <- map(counts_list, sum)

sample_list <- map(size_list, seq, from = 1, by = 1000)

counts_rarefy_list <- map2(counts_list, sample_list, rarefy)

counts_rarefy_df <- counts_rarefy_list %>%
  map(t) %>%
  map(as.data.frame) %>%
  map(rownames_to_column) %>%
  bind_rows(.id = "exp")

counts_rarefy_df <- counts_rarefy_df %>%
  factor_extractor() %>%
  separate(mouse, c("timepoint", "mouse"), sep = 1) %>%
  mutate(rowname = str_remove(rowname, "N") %>% as.numeric())

ggplot(counts_rarefy_df, aes(rowname, V1, group = exp, colour = response)) +
  geom_line() +
  theme(legend.position = "none") +
  labs(x = "sample size", y = "species") +
  facet_wrap(model ~ timepoint, scales = "free", ncol = 4)
