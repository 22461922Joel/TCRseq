
factors <- c("model", "timepoint", "mouse", "response", "primer")
working_path <- file.path(getwd(), "experiments/RNAseq")

data <- bind_rows(read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/cleaned_CDR3s.csv")),
                  read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp19 Final MixCR analysis//cleaned_CDR3s.csv"))) %>%
  factor_extractor() %>%
  filter(!str_detect(v_gene, "TRA")) %>%
  data_aa() %>%
  factor_extractor() %>%
  group_by(exp) %>%
  mutate(PID.fraction = PID.count_sum/sum(PID.count_sum))

#####
# working on V and J gene usage between responders and non responders
#####

library(missMDA)

data_vj <- bind_rows(read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/cleaned_CDR3s.csv")),
                  read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp19 Final MixCR analysis//cleaned_CDR3s.csv"))) %>%
  factor_extractor() %>%
  filter(!str_detect(v_gene, "TRA"), model == "AB1", timepoint == "0") %>%
  group_by(exp, v_gene, j_gene) %>%
  summarise(PID.count = sum(PID.count)) %>%
  filter(v_gene != "TRBV9", 
         v_gene != "TRBV8", 
         v_gene != "TRBV21", 
         v_gene != "TRBV23", 
         v_gene != "TRBV10", 
         j_gene != "TRBJ1-7", 
         j_gene != "TRBJ2-6") %>%
  unite("v_j", v_gene, j_gene, sep = "_") %>%
  group_by(exp) %>%
  mutate(PID.fraction = PID.count/sum(PID.count)) %>%
  select(-PID.count) %>%
  spread(v_j, PID.fraction) %>%
  column_to_rownames("exp")

library(ggcorrplot)

vj_corr <- cor(data_vj)

ggcorrplot(vj_corr, type = "upper") +
  theme(axis.text = element_blank())

corr_graph <- vj_corr[lower.tri(vj_corr)] %>%
  as.data.frame() %>%
  mutate(from = CombSet(colnames(vj_corr), m = 2)[,1],
         to = CombSet(colnames(vj_corr), m = 2)[,2]) %>%
  na.omit()

names(corr_graph) <- c("weight", "from", "to")

corr_graph <- corr_graph %>%
  select(from, to, weight) %>%
  filter(weight > 0.5) %>%
  graph_from_data_frame()

ggnet2(corr_graph, size = "degree")

V(corr_graph)[degree(corr_graph) == max(degree(corr_graph))]

nb <- estim_ncpPCA(data_vj)

data_vj.comp <- imputePCA(as.data.frame(data_vj), ncp = if_else(nb$ncp == 0, 2, as.double(nb$ncp)))$completeObs %>%
  scale() %>%
  as.data.frame() %>%
  rownames_to_column("exp") %>%
  factor_extractor()

data_vj_pca <- prcomp(data_vj.comp[-1:-6])

library(ggfortify)

autoplot(data_vj_pca, loadings = F, loadings.label = F, data = data_vj.comp, shape = "timepoint", colour = "response", size = 3) +
  theme_bw()

# heirarchical clustering for VJ genes and response

library(BBmisc)

model_tp_filter <- function(data) {
  data %>%
  filter(model == "RENCA", timepoint == "6")
}

data_vj <- bind_rows(read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/cleaned_CDR3s.csv")),
                     read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp19 Final MixCR analysis//cleaned_CDR3s.csv"))) %>%
  factor_extractor() %>%
  filter(!str_detect(v_gene, "TRA")) %>%
  model_tp_filter %>%
  group_by(exp, v_gene, j_gene) %>%
  summarise(PID.count = sum(PID.count)) %>%
  filter(v_gene != "TRBV9", 
         v_gene != "TRBV8", 
         v_gene != "TRBV21", 
         v_gene != "TRBV23", 
         v_gene != "TRBV10", 
         j_gene != "TRBJ1-7", 
         j_gene != "TRBJ2-6") %>%
  unite("v_j", v_gene, j_gene, sep = "_") %>%
  group_by(exp) %>%
  spread(v_j, PID.count) %>%
  BBmisc::normalize()

data_vj[is.na(data_vj)] <- 0

x1 <- data_vj %>%
  factor_extractor() %>%
  ungroup() %>%
  dplyr::select(-exp, -timepoint, -primer, -model, -response) %>%
  column_to_rownames("mouse") %>%
  dist() %>%
  hclust(method = "complete") %>%
  as.dendrogram()

y1 <- data_vj %>%
  gather("v_j_gene", "PID.count", -exp) %>%
  spread(exp, PID.count) %>%
  column_to_rownames("v_j_gene") %>%
  BBmisc::normalize() %>%
  dist() %>%
  hclust(method = "complete") %>%
  as.dendrogram()

x_dend <- ggplotGrob(x1 %>%
                       ggdendrogram() +
                       theme_void())

y1_placement <- dendro_data(y1)$labels$label

x1_placement <- dendro_data(x1)$labels$label

x_df_hm <- data_vj %>%
  gather("v_j_gene", "PID.count", -exp) %>%
  factor_extractor() %>%
  mutate(mouse = factor(mouse, levels = x1_placement),
         v_j_gene = factor(v_j_gene, levels = y1_placement),
         dummy = 1) %>%
  distinct()

x_rs <- ggplotGrob(x_df_hm %>%
                     dplyr::select(exp, response, dummy) %>%
                     distinct() %>%
                     ggplot(aes(exp, dummy)) +
                     geom_tile(aes(fill = response)) +
                     theme_void() +
                     scale_fill_manual(guide = F, values = c("NR" = "firebrick1",
                                                             "RS" = "deepskyblue1"),
                                       labels = c("NR" = "non-responder",
                                                  "RS" = "responder")))

x_gg <- ggplotGrob(x_df_hm %>%
                     ggplot(aes(mouse, v_j_gene)) +
                     geom_tile(aes(fill = PID.count)) +
                     theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                     scale_fill_viridis_c(name = "Z score") +
                     labs(y = "V-J pair"))

panel_id <- x_gg$layout[x_gg$layout$name == "panel", c("t", "l", "b")]

x_gg <- gtable_add_rows(x_gg, heights = unit(c(1, 0.3), "in"), 0)

x_gg <- gtable_add_grob(x_gg, grobs = list(x_dend, x_rs), t = c(1, 2), l = c(panel_id$l, panel_id$l), b = c(1, panel_id$b))

grid.newpage()

grid.draw(x_gg)

remove(x_gg, panel_id, x_rs, x_df_hm, x_dend, y1, x1)

#####
# get the threshold 0 and 1 clusters ready for joining to the entire dataset and making the final product, data_hmm
#####

T0_hmm <- read_csv(file.path(working_path, "T0_hmmsearch_left_join.csv")) %>%
  distinct() %>%
  group_by(aaSeqCDR3) %>%
  mutate(match_filter = best_domain_evalue == min(best_domain_evalue)) %>%
  filter(match_filter) %>%
  mutate(n_seq = n(),
         lowest_hmm = hmm_cluster == min(hmm_cluster)) %>%
  filter(lowest_hmm) %>%
  select(-lowest_hmm, -match_filter, -n_seq)

T1_hmm <- read_csv(file.path(getwd(), "experiments/RNAseq/hmmsearch_left_join.csv")) %>%
  distinct() %>%
  group_by(aaSeqCDR3) %>%
  mutate(match_filter = best_domain_evalue == min(best_domain_evalue)) %>%
  filter(match_filter) %>%
  mutate(n_seq = n(),
         lowest_hmm = hmm_cluster == min(hmm_cluster)) %>%
  filter(lowest_hmm) %>%
  select(-lowest_hmm, -match_filter, -n_seq)

names(T0_hmm) <- c("aaSeqCDR3", "T0_query_name", "T0_best_domain_evalue", "T0_hmm_cluster")

names(T1_hmm) <- c("aaSeqCDR3", "T1_query_name", "T1_best_domain_evalue", "T1_hmm_cluster")

data_hmm <- data %>%
  left_join(T1_hmm) %>%
  left_join(T0_hmm) %>%
  left_join(read_csv(file.path(working_path, "ed1_greedy_clusters.csv"))) %>%
  left_join(read_csv(file.path(working_path, "ed0_greedy_clusters.csv")))
#####
# directing a network of a significant cluster to see how the CDR3s move through time as well as by edit distance of 1 or less
#####

clust_4 <- data_hmm %>%
  filter(T1_hmm_cluster == 4) %>%
  ungroup() %>%
  select(aaSeqCDR3, timepoint) %>%
  group_by(timepoint) %>%
  group_split() %>%
  map(distinct)

clust_4_dist <- stringdist::stringdistmatrix(clust_4[[3]]$aaSeqCDR3, clust_4[[4]]$aaSeqCDR3, method = "lv")

colnames(clust_4_dist) <- clust_4[[4]]$aaSeqCDR3

row.names(clust_4_dist) <- clust_4[[3]]$aaSeqCDR3

clust_4_dist_df_3 <- clust_4_dist %>%
  as.data.frame() %>%
  rownames_to_column("from") %>%
  gather("to", "edit_distance", -from) %>%
  filter(edit_distance <= 1) %>%
  mutate(timepoint = 3)

clust_4_dist_df_2 <- clust_4_dist %>%
  as.data.frame() %>%
  rownames_to_column("from") %>%
  gather("to", "edit_distance", -from) %>%
  filter(edit_distance <= 1) %>%
  mutate(timepoint = 2)

clust_4_dist_df <- clust_4_dist %>%
  as.data.frame() %>%
  rownames_to_column("from") %>%
  gather("to", "edit_distance", -from) %>%
  filter(edit_distance <= 1) %>%
  mutate(timepoint = 1)

clust_4_graph <- reduce(list(clust_4_dist_df, clust_4_dist_df_2, clust_4_dist_df_3), bind_rows) %>%
  graph_from_data_frame()

l <- layout_in_circle(clust_4_graph)

plot(clust_4_graph, 
     vertex.label = NA, 
     edge.arrow.size = .2, 
     vertex.size = 5, 
     edge.color = E(clust_4_graph)$timepoint,
     layout = l)

#####
#####
# checking heirarchical clustering of clusters
#####

library(BBmisc)

x_df <- data_hmm %>%
  filter(timepoint == "6" & model == "RENCA") %>%
  group_by(exp, T0_hmm_cluster) %>%
  summarise(count = sum(PID.count_sum)) %>%
  na.omit() %>%
  mutate(fraction = count/sum(count)) %>%
  select(-fraction) %>%
  group_by(T0_hmm_cluster) %>%
  filter(n() == 14) %>%
  spread(T0_hmm_cluster, count) %>%
  BBmisc::normalize() %>%
  gather("T0_hmm_cluster", "fraction", -exp)

x1 <- x_df %>%
  spread(T0_hmm_cluster, fraction) %>%
  factor_extractor() %>%
  dplyr::select(-exp, -timepoint, -primer, -model, -response) %>%
  column_to_rownames("mouse") %>%
  dist() %>%
  hclust(method = "complete") %>%
  as.dendrogram()

x_df_na <- x_df

x_df_na$fraction[is.na(x_df_na$fraction)] <- 0

y1 <- x_df_na %>%
  spread(exp, fraction) %>%
  column_to_rownames("T0_hmm_cluster") %>%
  BBmisc::normalize() %>%
  dist() %>%
  hclust(method = "complete") %>%
  as.dendrogram()

x_dend <- ggplotGrob(x1 %>%
                       ggdendrogram() +
                       theme_void())

y1_placement <- dendro_data(y1)$labels$label

x1_placement <- dendro_data(x1)$labels$label

x_df_hm <- x_df %>%
  factor_extractor() %>%
  mutate(mouse = factor(mouse, levels = x1_placement),
         T0_hmm_cluster = factor(T0_hmm_cluster, levels = y1_placement),
         dummy = 1) %>%
  distinct()

x_rs <- ggplotGrob(x_df_hm %>%
                     dplyr::select(exp, response, dummy) %>%
                     distinct() %>%
                     ggplot(aes(exp, dummy)) +
                     geom_tile(aes(fill = response)) +
                     theme_void() +
                     scale_fill_manual(guide = F, values = c("NR" = "firebrick1",
                                                             "RS" = "deepskyblue1"),
                                       labels = c("NR" = "non-responder",
                                                  "RS" = "responder")))

x_gg <- ggplotGrob(x_df_hm %>%
                     ggplot(aes(mouse, T0_hmm_cluster)) +
                     geom_tile(aes(fill = fraction)) +
                     theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                     scale_fill_viridis_c(name = "Z score") +
                     labs(y = "TCRb CDR3 sequence"))

panel_id <- x_gg$layout[x_gg$layout$name == "panel", c("t", "l", "b")]

x_gg <- gtable_add_rows(x_gg, heights = unit(c(1, 0.3), "in"), 0)

x_gg <- gtable_add_grob(x_gg, grobs = list(x_dend, x_rs), t = c(1, 2), l = c(panel_id$l, panel_id$l), b = c(1, panel_id$b))

grid.newpage()

grid.draw(x_gg)

remove(x_gg, panel_id, x_rs, x_df_hm, x_dend, y1, x1, x_df_na, x_df)

#####
# figuring out what clusters only belong to responders
#####

clust_res_freqtable <- table(data_hmm$response, data_hmm$T0_hmm_cluster)

clust_res_df <- as.data.frame.table(clust_res_freqtable) %>%
  spread(Var1, Freq) %>%
  mutate(dif = RS - NR) %>%
  filter(dif > 100 | dif < 0)

ggplot(clust_res_df, aes(Var2, dif)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))

#####
# for 1 cluster find the distribution of clones that make up that cluster over time
#####

cluster5 <- data_hmm %>%
  filter(T1_hmm_cluster == 5) %>%
  factor_extractor()

cluster5_full <- cluster5 %>%
  dplyr::select(aaSeqCDR3, PID.count_sum, PID_cluster, exp) %>%
  spread(aaSeqCDR3, PID.count_sum)

cluster5_full[is.na(cluster5_full)] <- 0

cluster5_full <- cluster5_full %>%
  gather("aaSeqCDR3", "PID.count", -exp, -PID_cluster) %>%
  factor_extractor() %>%
  filter(model == "RENCA")

ggplot(cluster5_full, aes(reorder(exp, PID.count), PID.count, group = interaction(aaSeqCDR3, as.character(PID_cluster)))) +
  geom_area(aes(fill = aaSeqCDR3), position = "stack") +
  scale_color_brewer() +
  geom_line(position = "stack", colour = "black") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_blank()) +
  facet_wrap(~ response + timepoint, scales = "free_x", ncol = 4) +
  labs(x = "mouse", y = "RNA transcripts from cluster 5")

#####
# trying to get all the top clusters mapped by the top clone of every mouse
#####

test_dat <- data_hmm %>%
  filter(model == "RENCA") %>%
  group_by(exp) %>%
  filter(PID.count_sum == max(PID.count_sum)) %>%
  ungroup() %>%
  select(T0_hmm_cluster) %>%
  unique()

top_dat_AB1 <- data_hmm %>%
  inner_join(test_dat) %>%
  filter(!is.na(T0_hmm_cluster), model == "RENCA") %>%
  select(PID.count_sum, T0_hmm_cluster, exp, aaSeqCDR3) %>%
  spread(T0_hmm_cluster, PID.count_sum)

top_dat_AB1[is.na(top_dat_AB1)] <- 0

top_dat_AB1 <- top_dat_AB1 %>%
  gather("T0_hmm_cluster", "PID.count_sum", -exp, -aaSeqCDR3) %>%
  group_by(exp, T0_hmm_cluster) %>%
  summarise(PID.count_sum = sum(PID.count_sum)) %>%
  factor_extractor()

ggplot(top_dat_AB1, aes(reorder(exp, PID.count_sum), PID.count_sum, group = T0_hmm_cluster)) +
  geom_area(aes(fill = T0_hmm_cluster), position = "stack") +
  scale_fill_viridis_d() +
  geom_line(position = "stack", colour = "black") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  facet_wrap(~ response + timepoint, scales = "free_x", ncol = 4)

#####
#####
# joining data to the McPAS database, code will need to be updated to work with newer raw data set
#####

data_McPAS <- data_hmm %>%
  left_join(hmmsearch_filtered) %>%
  filter(hmm_cluster == 4 | hmm_cluster == 5) %>%
  filter(!is.na(exp)) %>%
  group_by(aaSeqCDR3, pathology, ID, timepoint, model, response) %>%
  summarise(n_mice = n_distinct(exp),
            PID.count = sum(PID.count))

ggplot(McPAS_filtered %>% 
         mutate(pathology = str_remove(pathology, "Human") %>% 
                  str_remove("Experimental_autoimmune_encephalomyelitis") %>%
                  str_remove("_immunodeficiency_virus") %>%
                  str_remove("Lymphocytic_choriomeningitis_virus")), aes(x = pathology)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot(data_McPAS %>% 
         ungroup() %>%
         mutate(pathology = str_remove(pathology, "Human") %>% 
                  str_remove("Experimental_autoimmune_encephalomyelitis") %>%
                  str_remove("_immunodeficiency_virus") %>%
                  str_remove("Lymphocytic_choriomeningitis_virus")), aes(pathology)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggplot(data_McPAS, aes(n_mice, fill = response)) +
  geom_histogram(binwidth = 1) +
  labs(title = "check for groups of mice preferentially selecting Ag specific TCRs-experiment level") +
  facet_grid(model ~ timepoint)

#####
# make a plot looking at the top clusters used in both models with the number of mice involved 
#####

data_clustered_reduced <- data_hmm %>%
  group_by(T0_hmm_cluster) %>%
  mutate(total_PID.count = sum(PID.count_sum)) %>%
  group_by(T0_hmm_cluster, model, timepoint, response) %>%
  summarise(PID.count = sum(PID.count_sum),
            n_mice = n_distinct(mouse),
            total_PID.count = mean(total_PID.count))

ggplot(data_clustered_reduced %>% 
         na.omit() %>% 
         filter(total_PID.count > 50000), aes(timepoint, 
                                              PID.count, 
                                              group = interaction(T0_hmm_cluster, response), 
                                              shape = response)) +
  geom_point(aes(colour = n_mice), size = 5) +
  theme(panel.background = element_rect(fill = "grey")) +
  geom_line() +
  scale_colour_viridis_c() +
  facet_grid(model ~ T0_hmm_cluster) +
  theme_bw()

#####
# makes the images of the top culsters for each model. needs work because it hides how much any given
# mouse has an effect on the total
#####
library(vegan)

clust_ab1 <- data_hmm %>%
  filter(model == "AB1") %>%
  group_by(T0_hmm_cluster, response, timepoint, model) %>%
  summarise(mean_count = mean(PID.count_sum),
            sd_count = sd(PID.count_sum),
            PID.count = sum(PID.count_sum),
            n_mice = n_distinct(mouse)) %>%
  na.omit() %>%
  group_by(T0_hmm_cluster) %>%
  mutate(filter_count = sum(PID.count),
         filter_mice = min(n_mice),
         w_sd = wisconsin(sd_count)) %>%
  filter(filter_count > (mean(.$filter_count) + (4 * sd(.$filter_count))), filter_mice >= 4)

length(unique(clust_ab1$hmm_cluster))

ggplot(clust_ab1, aes(timepoint, PID.count, colour = n_mice, shape = response, group = interaction(T0_hmm_cluster, response))) +
  geom_point(aes(colour = n_mice), size = 4) +
  geom_line(colour = "black") +
  scale_colour_viridis_c() +
  theme_bw() +
  theme(panel.background = element_rect(colour = "grey51")) +
  facet_wrap(~ T0_hmm_cluster) +
  labs(title = "AB1")

clust_renca <- data_hmm %>%
  filter(model == "RENCA") %>%
  group_by(T1_hmm_cluster, response, timepoint, model) %>%
  summarise(mean_count = mean(PID.count_sum),
            sd_count = sd(PID.count_sum),
            PID.count = sum(PID.count_sum),
            n_mice = n_distinct(mouse)) %>%
  na.omit() %>%
  group_by(T1_hmm_cluster) %>%
  mutate(filter_count = sum(PID.count),
         filter_mice = min(n_mice)) %>%
  filter(filter_count > (mean(.$filter_count) + (2 * sd(.$filter_count))), filter_mice >= 4)

length(unique(clust_renca$PID_cluster))

ggplot(clust_renca, aes(timepoint, PID.count, colour = n_mice, shape = response, group = interaction(T1_hmm_cluster, response))) +
  geom_point(aes(colour = n_mice), size = 4) +
  geom_line(colour = "black") +
  scale_colour_viridis_c() +
  theme_bw() +
  theme(panel.background = element_rect(colour = "grey51")) +
  facet_wrap(~ T1_hmm_cluster) +
  labs(title = "RENCA")

#####
# makes the dimensionality reduction for the clusters 
#####
clust_pca_df <- data_hmm %>%
  group_by(T0_hmm_cluster, exp) %>%
  summarise(PID.count = sum(PID.count_sum)) %>%
  na.omit() %>%
  spread(T0_hmm_cluster, PID.count, fill = 0) %>%
  mutate_if(is.numeric, scale)

clust_pca <- prcomp(clust_pca_df[,-1])

summary(clust_pca)

clust_pca_pcs <- clust_pca$x %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  bind_cols(clust_pca_df[,1]) %>%
  factor_extractor()

ggplot(clust_pca_pcs, aes(PC1, PC2, colour = model)) +
  geom_point() +
  theme_bw()

clones_isoMDS_df <- data_hmm %>%
  dplyr::select(exp, PID.count, aaSeqCDR3) %>%
  group_by(exp, aaSeqCDR3) %>%
  summarise(PID.count = sum(PID.count)) %>%
  spread(aaSeqCDR3, PID.count, fill = 0) %>%
  column_to_rownames("exp") %>%
  wisconsin()

clones_isoMDS <- isoMDS(clones_isoMDS_df)

clones_isoMDS_comps <- clones_isoMDS$points %>%
  as.data.frame() %>%
  rownames_to_column("exp") %>%
  factor_extractor()

ggplot(clones_isoMDS_comps, aes(V1, V2, colour = model)) +
  geom_point() +
  theme_bw()

fraction_isoMDS_df <- data_clustered %>%
  dplyr::select(exp, PID.fraction, aaSeqCDR3) %>%
  group_by(exp, aaSeqCDR3) %>%
  summarise(PID.fraction = sum(PID.fraction)) %>%
  spread(aaSeqCDR3, PID.fraction, fill = 0) %>%
  column_to_rownames("exp") %>%
  vegdist()

fraction_isoMDS <- isoMDS(fraction_isoMDS_df)

fraction_isoMDS_comps <- fraction_isoMDS$points %>%
  as.data.frame() %>%
  rownames_to_column("exp") %>%
  factor_extractor()

ggplot(fraction_isoMDS_comps, aes(V1, V2, colour = model)) +
  geom_point() +
  theme_bw()

clust_isoMDS_df <- data_hmm %>%
  group_by(hmm_cluster, exp) %>%
  summarise(PID.count = sum(PID.count)) %>%
  na.omit() %>%
  spread(hmm_cluster, PID.count, fill = 0) %>%
  column_to_rownames("exp") %>%
  decostand(method = "total") %>%
  vegdist()

clust_isoMDS <- isoMDS(clust_isoMDS_df)

clust_isoMDS_comps <- clust_isoMDS$points %>%
  as.data.frame() %>%
  rownames_to_column("exp") %>%
  factor_extractor()

tumour_colour_scale <- scale_colour_manual(values = c("AB1" = "violetred",
                                                        "RENCA" = "seagreen"),
                                             labels = c("AB1" = "AB1",
                                                        "RENCA" = "RENCA"))

ggplot(clust_isoMDS_comps, aes(V1, V2, colour = model)) +
  geom_point() +
  theme_bw() +
  labs(y = "MDS2", x = "MDS1", title = "A") +
  tumour_colour_scale

#####
# produce network of greedy clusters with how the hmm fit into them
#####
cluster_4 <- data_hmm %>%
  ungroup() %>%
  filter(PID_cluster == 4, !is.na(T0_hmm_cluster), T0_hmm_cluster < 200) %>%
  dplyr::select(aaSeqCDR3) %>%
  distinct() %>%
  mutate(residue_group = str_replace_all(aaSeqCDR3, c("[IV]" = "J", 
                                                      "[FY]" = "O", 
                                                      "[QE]" = "U", 
                                                      "[RK]" = "B",
                                                      "[AST]" = "X",
                                                      "[LM]" = "Y",
                                                      "[ND]" = "Z")))

cluster_4_dist_matrix <- data.frame(CombSet(cluster_4$aaSeqCDR3, m = 2), 
                                    lv = stringdistmatrix(cluster_4$residue_group) %>%
                                      as.vector()) %>%
  filter(lv < 2)

table(cluster_4_dist_matrix$lv)

cluster_4_vertices <- data_hmm %>%
  ungroup() %>%
  filter(PID_cluster == 4, !is.na(T0_hmm_cluster), T0_hmm_cluster < 200) %>%
  select(aaSeqCDR3, T0_hmm_cluster, T0_cluster) %>%
  distinct()

cluster_4_vertices$T0_hmm_cluster[is.na(cluster_4_vertices$T0_hmm_cluster)] <- "not assigned"

cluster_4_vertices$T0_cluster[is.na(cluster_4_vertices$T0_cluster)] <- "not assigned"

cluster_4_graph <- graph_from_data_frame(cluster_4_dist_matrix, directed = F, vertices = cluster_4_vertices)

E(cluster_4_graph)$color <- E(cluster_4_graph)$lv

V(cluster_4_graph)

vertex.attributes(cluster_4_graph)

E(cluster_4_graph)$color[E(cluster_4_graph)$color == 0] <- 3

ggnet2(cluster_4_graph, edge.color = "color", node.size = 3, node.color = V(cluster_4_graph)$T0_hmm_cluster, edge.alpha = 0.5) +
  scale_color_viridis_d()

#####
# produces volcano plot
#####
data_meta <- data_clustered %>%
  dplyr::select(model, timepoint, response, PID_cluster, PID.count, mouse) %>%
  group_by(model, timepoint, response, PID_cluster, mouse) %>%
  summarise(cluster_count = sum(PID.count)) %>%
  group_by(PID_cluster, model, timepoint, response) %>%
  mutate(n_mice = n_distinct(mouse),
         mean_count = mean(cluster_count)) %>%
  group_by(PID_cluster, model, timepoint) %>%
  mutate(n_groups = n_distinct(response),
         m_mouse = if_else(min(n_mice) == 1 | min(n_groups) == 1, 1, as.double(min(n_mice))),
         mean_dif = max(mean_count) - min(mean_count)) %>%
  filter(m_mouse != 1, mean_dif > 10) %>%
  na.omit() %>%
  dplyr::select(-mouse, -n_groups, -n_mice, -m_mouse, -mean_dif, -mean_count) %>%
  group_by(timepoint, PID_cluster, model, response) %>%
  nest() %>%
  spread(key = response, value = data) %>%
  mutate(t_test = map2(RS, NR, ~{t.test(.x$cluster_count, .y$cluster_count) %>% tidy()}),
         NR = map(NR, nrow),
         RS = map(RS, nrow)) %>%
  unnest()

data_log2_foldchange <- data_clustered %>%
  dplyr::select(model, timepoint, response, PID_cluster, PID.count, mouse) %>%
  group_by(model, timepoint, response, PID_cluster, mouse) %>%
  summarise(cluster_count = sum(PID.count)) %>%
  group_by(PID_cluster, model, timepoint, response) %>%
  mutate(n_mice = n_distinct(mouse)) %>%
  group_by(PID_cluster, model, timepoint) %>%
  mutate(n_groups = n_distinct(response),
         m_mouse = if_else(min(n_mice) == 1 | min(n_groups) == 1, 1, as.double(min(n_mice)))) %>%
  filter(m_mouse != 1) %>%
  group_by(PID_cluster, model, timepoint, response) %>%
  summarise(mean_cluster_count = mean(cluster_count)) %>%
  spread(response, mean_cluster_count) %>%
  mutate(log2_foldchange = log2(RS) - log2(NR)) %>%
  dplyr::select(-RS, -NR) %>%
  left_join(data_meta) %>%
  na.omit()

ggplot(data_log2_foldchange, aes(log2_foldchange, -log10(p.value), colour = model)) +
  geom_hline(yintercept = -log10(0.05), colour = "red", linetype = "dashed") +
  geom_point() +
  facet_grid(model ~ timepoint)

ggplot(data_log2_foldchange, aes(NR, RS)) +
  geom_jitter(alpha = 0.1) +
  labs(x = "number of non-responding mice", y = "number of responding mice") +
  geom_density_2d(h = c(3, 3), size = 1, colour = "darkorange") +
  facet_grid(model ~ timepoint) +
  theme_bw() +
  scale_x_continuous(breaks = 2:8) +
  scale_y_continuous(breaks = 2:8)

#####

heatmap_dendrogram(aa_data, "timepoint", "0", "response")
# summary plots
#####

summary_TCRseq(data = data, file.path(getwd(), "experiments/RNAseq"))

summary <- read_csv(file.path(getwd(), "experiments/RNAseq/summary_stats.csv")) %>%
  factor_extractor()

summary$timepoint <- if_else(summary$timepoint == "0", "Pre CPB",
                             if_else(summary$timepoint == "2", "day 2",
                                     if_else(summary$timepoint == "4", "day 4", "day 6")))

summary$timepoint <- factor(summary$timepoint, levels = c("Pre CPB", "day 2", "day 4", "day 6", ordered = T))

summary$response <- factor(summary$response, levels = c("NR", "RS"), ordered = T)

response_fill_scale <- scale_fill_manual(values = c("NR" = "firebrick1",
                               "RS" = "deepskyblue1"),
                               labels = c("NR" = "non-responder",
                                          "RS" = "responder"))


response_colour_scale <- scale_colour_manual(values = c("NR" = "firebrick1",
                                                    "RS" = "deepskyblue1"))

table(summary$timepoint, summary$model, summary$response)

summary1 <- ggplot(summary, aes(timepoint, diversity, fill = response)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(group = response), label = "p.signif") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  theme_bw() +
  theme(axis.line = element_line(), legend.position = "bottom", legend.box = "vertical", legend.title = element_blank()) +
  response_fill_scale +
  facet_grid(model ~ .) +
  labs(title = "B")

summary2 <- ggplot(summary, aes(evenness, richness, colour = response)) +
  geom_point() +
  scale_y_log10(breaks = c(500, 2000, 8000, 32000, 128000)) +
  scale_x_log10(breaks = c(190000, 19000, 1900)) +
  facet_grid(model ~ timepoint) +
  scale_colour_manual(labels = c("NR" = "non-responder",
                                 "RS" = "responder"),
                      values = c("NR" = "firebrick1",
                                          "RS" = "deepskyblue1")) +
  labs(y = "unique clones", x = "total clones") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.title = element_blank())


#####
# vegan
#####

counts_list <- data %>%
  dplyr::select(exp, PID.count_sum, aaSeqCDR3) %>%
  spread(exp, PID.count_sum, fill = 0) %>%
  column_to_rownames("aaSeqCDR3") %>%
  as.list()

typeof(counts_list[1])

sum(counts_list[[1]])

size_list <- purrr::map(counts_list, sum)

sample_list <- purrr::map(size_list, seq, from = 1, by = 1000)

plan(multiprocess)

counts_rarefy_list <- future_map2(counts_list, sample_list, rarefy)

counts_rarefy_df <- counts_rarefy_list %>%
  purrr::map(t) %>%
  purrr::map(as.data.frame) %>%
  purrr::map(rownames_to_column) %>%
  bind_rows(.id = "exp") %>%
  factor_extractor() %>%
  mutate(rowname = str_remove(rowname, "N") %>% as.numeric())

counts_rarefy_df$timepoint <- if_else(counts_rarefy_df$timepoint == "0", "Pre CPB",
                             if_else(counts_rarefy_df$timepoint == "2", "day 2",
                                     if_else(counts_rarefy_df$timepoint == "4", "day 4", "day 6")))

counts_rarefy_df$timepoint <- factor(counts_rarefy_df$timepoint, levels = c("Pre CPB", "day 2", "day 4", "day 6", ordered = T))

counts_rarefy_df$response <- factor(counts_rarefy_df$response, levels = c("NR", "RS"), ordered = T)

library(scales)

summary3 <- ggplot(counts_rarefy_df, aes(rowname, V1, group = exp, colour = response)) +
  geom_line() +
  labs(x = "total clones", y = "unique clones", title = "B") +
  facet_grid(model ~ timepoint) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  coord_cartesian(xlim = c(0, max(counts_rarefy_df$rowname)), ylim = c(0, 15000)) +
  response_colour_scale

rm(counts, counts_list, counts_rarefy_df, counts_rarefy_list, raremax)

grd_layout <- rbind(c(1, 1, 2),
                    c(1, 1, 2))

grid.arrange(summary2, summary1, layout_matrix = grd_layout)

#####
#network statistics
#####

CDR3_network(aa_data, "0", "RS")

CDR3_network(aa_data, "0", "NR")

aa_graphs <- list()

for(i in 1:length(aa_list)) {
  aa_graphs[[i]] <- mouse_CDR3_graph(aa_list[[i]])
}

graph_stats_degree <- vector()

graph_stats_exp <- vector()

m <- 1

for (g in 1:length(aa_graphs)) {
  graph_stats_degree[[m]] <- sum(degree(aa_graphs[[g]]) > 0)
  # graph_stats_degree[m] <- igraph::V(aa_graphs[[g]])[
  #     !is.na(vertex_attr(aa_graphs[[g]],
  #                        name = tpmouse[h]))][degree(aa_graphs[[g]],
  #                                                    V(aa_graphs[[g]])[
  #                                                      !is.na(vertex_attr(aa_graphs[[g]],  name = tpmouse[h]))]) >= 1] %>%
  #     length()
  graph_stats_exp[[m]] <- unique(V(aa_graphs[[g]])$exp)
    m <- m + 1
}

graph_stats <- data.frame(exp = graph_stats_exp,
                          connected_nodes = graph_stats_degree)

summary <- summary %>%
  left_join(graph_stats)

write.csv(summary, "summary.csv")

ggplot(summary, aes(timepoint, connected_nodes, fill = response)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  response_colour_scale +
  response_fill_scale +
  stat_compare_means(aes(group = response), label = "p.signif") +
  labs(y = "number of clustered nodes", x = "time point")


setwd("D://RNAseq//mouse_networks")

for (i in 1:length(aa_list)) {
  
  g_titles <- aa_list[[i]] %>%
    select(tpmouse, response, timepoint) %>%
    distinct()
  
  png(filename = paste("Mouse ", g_titles$tpmouse, ", ",
                       "Response ", g_titles$response, ", ",
                       "timepoint ", g_titles$timepoint, ".png", sep = ""),
      width = 1527, height = 883)
  
  mouse_CDR3_network(aa_list[[i]])
  
  dev.off()
}

aa_modified <- aa_data %>%
  filter(PID.count_sum > 3)

degree_dfm <- data.frame(tpmouse = vector(length = nrow(summary)),
                        combined_connected = vector(length = nrow(summary)))

combined_graphs <- list(nr0, rs0, nr2, rs2, nr4, rs4, nr6, rs6)

combined_graphsm <- list(nr0m, rs0m, nr2m, rs2m, nr4m, rs4m, nr6m, rs6m)

m <- 1

for (i in 1:length(combined_graphsm)) {
  vert_names <- names(vertex.attributes(combined_graphsm[[i]]))[str_detect(names(vertex.attributes(combined_graphsm[[i]])), "tpmouse")]
  
  for (j in 1:length(vert_names)) {
    degree_dfm[m,1] <- vertex_attr(combined_graphsm[[i]], vert_names[j]) %>% unique() %>% na.omit()
    
    connected_nodes <- degree(combined_graphsm[[i]], V(combined_graphsm[[i]])[!is.na(vertex_attr(combined_graphsm[[i]], name = vert_names[j]))]) > 0
    
    degree_dfm[m,2] <- sum(connected_nodes)
    
    m <- m + 1
  }
}

datam <- data %>%
  filter(PID.count > 3)

summary_TCRseq(datam)

summarym <- read.csv("summary_stats.csv") %>%
  factor_separator()

summarym <- summarym %>%
  left_join(degree_dfm)

summarym$degreePerRichness <- summarym$combined_connected / summarym$richness

summarym$timepoint <- if_else(summarym$timepoint == "0", "Pre CPB",
                             if_else(summarym$timepoint == "2", "day 2",
                                     if_else(summarym$timepoint == "4", "day 4", "day 6")))

summarym$timepoint <- factor(summarym$timepoint, levels = c("Pre CPB", "day 2", "day 4", "day 6", ordered = T))

summarym$response <- factor(summarym$response, levels = c("NR", "RS"), ordered = T)

write.csv(summarym, "summary_less.csv")

ggplot(summarym, aes(timepoint, degreePerRichness * 100, fill = response)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(group = response), label = "p.signif") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  response_fill_scale +
  labs(y = "% connected nodes")

#####
# AB1 and last of RENCA library prep
#####

ab1 <- read.csv("D://data//experiments//RNAseq//output//ab1_RNA_RT.csv")

renca <- read.csv("D://data//experiments//RNAseq//output//renca_RNA_RT.csv")

IIID <- renca %>%
  bind_rows(ab1) %>%
  select(-sampled) %>%
  na.omit()

template <- read.csv("D://Sequencing running spreadsheet//primer_template.csv", stringsAsFactors = F)

names(template) <- c("primer", "primers")

IIID <- IIID %>%
  left_join(template)

IIID$Mouse_ID <- str_pad(IIID$Mouse_ID, 3, pad = "0")

IIID$Tissue <- toupper(IIID$Tissue)

IIID$PCR <- "PCR2"

IIID <- IIID %>%
  unite("sample", Tissue, Mouse_ID, Group, PCR)

date <- data.frame(octet = unique(IIID$octet),
                   date_temp1 = unique(IIID$octet),
                   date_temp2 = "Jun-2019") %>%
  unite("date", date_temp1, date_temp2, sep = "-")

IIID <- IIID %>%
  left_join(date)

write_excel_csv(IIID, "D://data//experiments//RNAseq//output//IIID_ab1.csv")

#####