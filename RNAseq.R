
factors <- c("model", "timepoint", "mouse", "response", "primer")

data <- bind_rows(read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/cleaned_CDR3s.csv")),
                  read_csv(file.path(getwd(), "experiments/RNAseq/1239Shp19 Final MixCR analysis//cleaned_CDR3s.csv"))) %>%
  factor_extractor()

length(unique(data$exp))

data_clustered <- data %>%
  left_join(as.data.frame(pid_clusters)) %>%
  arrange(desc(PID.fraction), desc(PID_cluster))

clust_dif <- data_clustered %>%
  group_by(PID_cluster) %>%
  summarise(n_mice = n_distinct(exp),
            unique_clones = n_distinct(aaSeqCDR3),
            total_clones = sum(PID.count)) %>%
  group_by(PID_cluster) %>%
  mutate(total_fraction = total_clones/sum(total_clones),
         unique_fraction = unique_clones/sum(unique_clones))

clust_NA <- clust_dif %>%
  filter(is.na(PID_cluster))

ggplot(clust_dif %>% na.omit(), aes(n_mice)) +
  geom_density()

data_clustered_reduced <- data_clustered %>%
  group_by(PID_cluster) %>%
  mutate(total_PID.count = sum(PID.count)) %>%
  group_by(PID_cluster, model, timepoint, response) %>%
  summarise(PID.count = sum(PID.count),
            n_mice = n_distinct(mouse),
            total_PID.count = mean(total_PID.count))

test <- data_clustered_reduced %>% na.omit() %>% filter(total_PID.count > 50000)

ggplot(data_clustered_reduced %>% 
         na.omit() %>% 
         filter(total_PID.count > 50000), aes(timepoint, 
                                              PID.count, 
                                              group = interaction(PID_cluster, response), 
                                              shape = response)) +
  geom_point(aes(colour = n_mice), size = 5) +
  theme(panel.background = element_rect(fill = "grey")) +
  geom_line() +
  scale_colour_viridis_c() +
  facet_grid(model ~ PID_cluster)

clust_pca_df <- data_clustered %>%
  group_by(PID_cluster, exp) %>%
  summarise(PID.count = sum(PID.count)) %>%
  na.omit() %>%
  spread(PID_cluster, PID.count, fill = 0) %>%
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

clones_isoMDS_df <- data_clustered %>%
  dplyr::select(exp, PID.count, aaSeqCDR3) %>%
  group_by(exp, aaSeqCDR3) %>%
  summarise(PID.count = sum(PID.count)) %>%
  spread(aaSeqCDR3, PID.count, fill = 0) %>%
  column_to_rownames("exp") %>%
  vegdist()

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

clust_isoMDS_df <- data_clustered %>%
  group_by(PID_cluster, exp) %>%
  summarise(PID.count = sum(PID.count)) %>%
  na.omit() %>%
  spread(PID_cluster, PID.count, fill = 0) %>%
  column_to_rownames("exp") %>%
  vegdist()

clust_isoMDS <- isoMDS(clust_isoMDS_df)

clust_isoMDS_comps <- clust_isoMDS$points %>%
  as.data.frame() %>%
  rownames_to_column("exp") %>%
  factor_extractor()

ggplot(clust_isoMDS_comps, aes(V1, V2, colour = model)) +
  geom_point() +
  theme_bw()

cluster_1 <- all_clusters %>%
  filter(cluster == "cluster10") %>%
  mutate(residue_group = str_replace_all(aaSeqCDR3, c("[IV]" = "J", 
                                                      "[FY]" = "O", 
                                                      "[QE]" = "U", 
                                                      "[RK]" = "B",
                                                      "[AST]" = "X",
                                                      "[LM]" = "Y",
                                                      "[ND]" = "Z")))

cluster_1_dist_matrix <- data.frame(CombSet(cluster_1$aaSeqCDR3, m = 2), 
                                    lv = stringdistmatrix(cluster_1$residue_group) %>%
                                      as.vector()) %>%
  filter(lv < 2)

table(cluster_1_dist_matrix$lv)

cluster_1_graph <- graph_from_data_frame(cluster_1_dist_matrix, directed = F)

E(cluster_1_graph)$color <- E(cluster_1_graph)$lv

V(cluster_1_graph)

E(cluster_1_graph)$color[E(cluster_1_graph)$color == 0] <- 3

ggnet2(cluster_1_graph, edge.color = "color", node.size = 3)

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

table(data_log2_foldchange$NR, data_log2_foldchange$RS)

data_drilldown <- data_log2_foldchange %>%
  filter(p.value < 0.05)

data_required <- data %>%
  left_join(all_clusters) %>%
  factor_extractor() %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric()) %>%
  dplyr::select(model, timepoint, response, cluster, PID.fraction, mouse) %>%
  group_by(model, timepoint, response, cluster) %>%
  summarise(cluster_fraction = sum(PID.fraction)/124,
            n_mice = n_distinct(mouse)) %>%
  na.omit() %>%
  group_by(cluster, model) %>%
  mutate(n_groups = n_distinct(response)) %>%
  filter(n_groups == 1)

ggplot(data_required, aes(timepoint, cluster_fraction, colour = response, group = cluster)) +
  geom_point() +
  geom_line() +
  facet_grid(model ~ response)

data_logistic <- data %>%
  left_join(all_clusters) %>%
  factor_extractor() %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric(),
         response = if_else(response == "RS", 1, 0)) %>%
  dplyr::select(model, timepoint, response, cluster, PID.fraction, mouse) %>%
  group_by(model, timepoint, response, cluster, mouse) %>%
  summarise(cluster_fraction = sum(PID.fraction)) %>%
  group_by(cluster, model, timepoint, response) %>%
  mutate(n_mice = n_distinct(mouse)) %>%
  group_by(cluster, model, timepoint) %>%
  mutate(n_groups = n_distinct(response),
         m_mouse = if_else(min(n_mice) == 1 | min(n_groups) == 1, 1, as.double(min(n_mice)))) %>%
  filter(m_mouse != 1) %>%
  dplyr::select(-n_groups, -n_mice, -m_mouse) %>%
  spread(cluster, cluster_fraction) %>%
  group_by(model, timepoint) %>%
  group_split()

summary(data_logistic[[1]])

library(Amelia)
library(mlbench)
library(VIM)
library(FactoMineR)
library(naniar)

gg_miss_var(data_logistic[[1]])

data_logistic[[1]] <- data_logistic[[1]][,colSums(is.na(data_logistic[[1]]))<nrow(data_logistic[[1]])]

res <- summary(aggr(data_logistic[[1]], sortVar = T))$combinations

head(res[rev(order(res[,2])),])

matrixplot(data_logistic[[1]], sortby = 2)

data_miss <- data.frame(is.na(data_logistic[[1]]))

data_miss <- apply(X = data_miss, FUN = function(x) if(x) "m" else "o", MARGIN = c(1,2))

res_mca <- MCA(data_miss, graph = F)

plot(res_mca, invis = "ind", cex = 0.5)
  
data_l1 <- estim_ncpPCA(data_logistic[[1]][,5:25], method.cv = "Kfold")

data_l1_complete <- imputePCA()

missmap(data_logistic[[1]][,5:length(data_logistic[[1]])], col = c("blue", "red"), legend = F)

model_null <- glm(response ~ 1, data = data_logistic[[1]], family = binomial(link = "logit"))

model_full <- glm(response ~ `3` + `4` + `5` + `9` + `11` + `13` + `15`, data = data_logistic[[1]], family = binomial(link = "logit"))

step(model_null, scope = list(upper = model_full), direction = "both", test = "Chisq", data = data_logistic[[1]])

clusters_of_consequence <- all_clusters %>%
  right_join(data_clustered) %>%
  select(-clones, -repertoire_fraction)

data_coc <- data %>%
  left_join(clusters_of_consequence) %>%
  factor_extractor() %>%
  mutate(timepoint = as.numeric(timepoint)) %>%
  group_by(response, timepoint, cluster, model) %>%
  summarise(PID.fraction_sum = sum(PID.fraction)) %>%
  na.omit()

ggplot(data_meta, aes(clones, repertoire_fraction)) +
  geom_point()

table(data_clustered$timepoint, data_clustered$response, data_clustered$cluster)


data_clustered %>%
  filter(model == "AB1") %>%
  ggplot(aes(timepoint, PID.fraction_sum, group = interaction(cluster, response), colour = as.character(cluster))) +
  geom_line(aes(linetype = response), size = 2) +
  theme() +
  scale_colour_discrete(guide = F) +
  facet_grid(. ~ response)

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
                               "RS" = "deepskyblue1"))


response_colour_scale <- scale_colour_manual(values = c("NR" = "firebrick1",
                                                    "RS" = "deepskyblue1"))

table(summary$timepoint, summary$model, summary$response)

ggplot(summary %>% filter(exp != "AB1_2_39_RS_4"), aes(timepoint, diversity, fill = response)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(group = response), label = "p.signif") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  theme(axis.line = element_line()) +
  response_fill_scale +
  facet_grid(model ~ .)

ggplot(summary, aes(evenness, richness, colour = response)) +
  geom_point() +
  scale_y_log10(breaks = c(1000, 2000, 4000, 8000)) +
  scale_x_log10(breaks = c(190000, 19000, 1900)) +
  facet_grid(model ~ timepoint) +
  labs(y = "unique clones", x = "total clones") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.line = element_line()) +
  response_colour_scale


#####
# vegan
#####

counts_list <- data %>%
  data_aa() %>%
  select(exp, PID.count_sum, aaSeqCDR3) %>%
  spread(exp, PID.count_sum, fill = 0) %>%
  column_to_rownames("aaSeqCDR3") %>%
  as.list()

typeof(counts_list[1])

sum(counts_list[[1]])

size_list <- map(counts_list, sum)

sample_list <- map(size_list, seq, from = 1, by = 1000)

plan(multiprocess)

counts_rarefy_list <- future_map2(counts_list, sample_list, rarefy)

counts_rarefy_df <- counts_rarefy_list %>%
  map(t) %>%
  map(as.data.frame) %>%
  map(rownames_to_column) %>%
  bind_rows(.id = "exp") %>%
  factor_extractor() %>%
  mutate(rowname = str_remove(rowname, "N") %>% as.numeric())

library(scales)

ggplot(counts_rarefy_df %>% filter(exp != "AB1_2_39_RS_4"), aes(rowname, V1, group = exp, colour = response)) +
  geom_line() +
  scale_x_continuous(labels = scientific) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  labs(x = "sample size", y = "species") +
  facet_grid(model ~ timepoint)

rm(counts, counts_list, counts_rarefy_df, counts_rarefy_list, raremax)

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