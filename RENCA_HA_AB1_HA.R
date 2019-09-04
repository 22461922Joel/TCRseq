library(vegan)

library(MASS)

library(devtools)

library(ggbiplot)

library(ggalt)

working_path <- "experiments/RENCA_HA_AB1_HA_exp12_exp13"

clean_data(file.path(getwd(), working_path, "1239Shp12_20190128_Mixcr analysis results/PID summary"))

clean_data(file.path(getwd(), working_path, "1239Shp13 Mixcr results/PID summary"))

exp11 <- read_csv(file.path(getwd(), working_path, "1239Shp12_20190128_Mixcr analysis results/cleaned_CDR3s.csv")) %>%
  data_aa()

exp15 <- read_csv(file.path(getwd(), working_path, "1239Shp13 Mixcr results/cleaned_CDR3s.csv")) %>%
  data_aa()

exp11$exp <- exp11$exp %>%
  str_remove("_12")

exp11 <- exp11 %>%
  separate(exp, into = c("exp_num", "tissue", "mouse", "population"), sep = "_") %>%
  mutate(flank = "O") %>%
  unite("exp", exp_num, mouse, tissue, flank, population, sep = "_")

exp15 <- exp15 %>%
  mutate(flank = "O") %>%
  separate(exp, into = c("exp_num", "mouse", "tissue", "population"), sep = "_") %>%
  unite("exp", exp_num, mouse, tissue, flank, population, sep = "_")

RENCA_AB1 <- bind_rows(exp11, exp15, jk41)

factors <- c("exp_num", "mouse", "tissue", "flank", "population")

PID_matrix <- RENCA_AB1 %>%
  factor_extractor() %>%
  filter(tissue == "T", population == "CD8") %>%
  dplyr::select(exp, PID.count_sum, aaSeqCDR3) %>%
  spread(exp, PID.count_sum) %>%
  column_to_rownames(var = "aaSeqCDR3")

PID_matrix[is.na(PID_matrix)] <- 0

PID_matrix_t <- t(PID_matrix)

dim(PID_matrix_t)

RENCA_AB1_dist <- vegdist(PID_matrix_t)

RENCA_AB1_mds0 <- isoMDS(RENCA_AB1_dist)


RENCA_AB1_mds0_df <- RENCA_AB1_mds0[[1]] %>%
  as.data.frame() %>%
  rownames_to_column("exp") %>%
  factor_extractor() %>%
  filter(tissue == "T", population == "CD8")

ggplot(RENCA_AB1_mds0_df, aes(V1, V2, colour = exp_num)) +
  geom_point() +
  geom_encircle(data = RENCA_AB1_mds0_df %>%
                  filter(exp_num == "JCCL4Exp11"), aes(V1, V2)) +
  geom_encircle(data = RENCA_AB1_mds0_df %>%
                  filter(exp_num == "JCCL4Exp15.2"), aes(V1, V2)) +
  geom_encircle(data = RENCA_AB1_mds0_df %>%
                  filter(exp_num == "JK4.1"), aes(V1, V2))

PID_matrix_t_reduced <- PID_matrix_t[,colSums(!is.na(PID_matrix_t)) > 3]

dg <- estim_ncpPCA(X = PID_matrix_t_reduced, ncp.min = 2)

dg$ncp

PID_matrix_t_reduced_c <- imputePCA(PID_matrix_t_reduced, ncp = dg$ncp)$completeObs

PID_matrix_t_reduced_c[PID_matrix_t_reduced_c < 0] <- 0

dat_dist <- vegdist(PID_matrix_t_reduced_c)

dat_pca <- prcomp(PID_matrix_t_reduced_c)

dat_pca_group <- dat_pca$x %>%
  rownames() %>%
  str_remove("_.*")

ggbiplot(dat_pca, var.axes = F, groups = dat_pca_group)

dat_mds0 <- isoMDS(dat_dist)

dat_mds0_df <- dat_mds0[[1]] %>%
  as.data.frame() %>%
  rownames_to_column("exp") %>%
  factor_extractor() %>%
  filter(tissue == "T", population == "CD8")

ggplot(dat_mds0_df, aes(V1, V2, colour = exp_num)) +
  geom_point() +
  geom_encircle(data = dat_mds0_df %>%
                  filter(exp_num == "JCCL4Exp11"), aes(V1, V2)) +
  geom_encircle(data = dat_mds0_df %>%
                  filter(exp_num == "JCCL4Exp15.2"), aes(V1, V2)) +
  geom_encircle(data = dat_mds0_df %>%
                  filter(exp_num == "JK4.1"), aes(V1, V2))
