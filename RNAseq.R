# clean_data("D://data//experiments//RNAseq//1239Shp15b_Joost_final//PID\ summary")

setwd("D://data//experiments//RNAseq")

data <- read.csv("cleaned_CDR3s.csv", stringsAsFactors = F)
# 
#  data$exp <- str_remove(data$exp, "mTCR ") %>%
#    str_replace(" ", "-") %>%
#    str_replace_all("627-RS", "627-RS-4") %>%
#    str_replace_all("630-NR", "630-NR-6") %>%
#    str_replace_all("633-RS", "633-RS-7") %>%
#    str_replace_all("640-RS", "640-RS-8")
#  
#  write.csv(data, "cleaned_CDR3s.csv")

factor_separator <- function(data) {
  data %>% separate(exp, into = c("tpmouse", "response", "primer"), sep = "-", remove = F) %>%
    separate(tpmouse, into = c("timepoint", "mouse"), sep = 1, remove = F)
}

data <- data %>%
  factor_separator()

aa_data <- data_aa(data) %>% factor_separator()
# summary plots
#####

PID <- PID_control("D://data//experiments//RNAseq//1239Shp15b_Joost_final//PID\ summary")
# 
# summary_TCRseq(data = data)

summary <- read.csv("summary.csv", stringsAsFactors = F, row.names = NULL) %>%
  factor_separator() %>%
  select(-X)

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
#network statistics
#####

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

nr0 <- graph_func(aa_data, "0", "NR")
rs0 <- graph_func(aa_data, "0", "RS")
nr2 <- graph_func(aa_data, "2", "NR")
rs2 <- graph_func(aa_data, "2", "RS")
nr6 <- graph_func(aa_data, "6", "NR")
rs6 <- graph_func(aa_data, "6", "RS")
rs4 <- graph_func(aa_data, "4", "RS")

aa_modified <- aa_data %>%
  filter(PID.count_sum > 3)

nr4 <- graph_func(aa_modified, "4", "NR")
nr0m <- graph_func(aa_modified, "0", "NR")
rs0m <- graph_func(aa_modified, "0", "RS")
nr2m <- graph_func(aa_modified, "2", "NR")
rs2m <- graph_func(aa_modified, "2", "RS")
nr6m <- graph_func(aa_modified, "6", "NR")
rs6m <- graph_func(aa_modified, "6", "RS")
rs4m <- graph_func(aa_modified, "4", "RS")
nr4m <- nr4

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