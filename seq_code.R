library(ggpubr)
library(data.table)
# network
library(RColorBrewer)
library(gridExtra)
library(gtable)
library(grid)
library(sna)
library(network)
library(GGally)
library(treemap)
library(igraph)
library(stringdist)
# phylogram
library(ape)
#tsne
library(stringdist)
library(Rtsne)
library(plotly)
#
library(reshape2)
library(reshape)
library(dendextend)
library(data.tree)
# heatmap + dendrgoram
library(grid)
library(gtable)
library(ggdendro)
library(missMDA)
# standard
library(DescTools)
library(divo)
library(readxl)
library(tidyverse)
library(colorRamps)

#####
# data cleaning
#####

# clean_data expects data shaped as it comes from IIID as of June 2019, ie. <experiment>/PID summary folder structure.
# it returns 2 .csv files with clean data and rejected CDR3s for reference

clean_data <- function(directory) {
  
  setwd(directory)
  
  names_exp <- dir()[str_detect(dir(), " miXCR clones with PID read cutoff 1.tsv")] %>%
    str_remove(" miXCR clones with PID read cutoff 1.tsv")
  
  raw_exp <- dir()[str_detect(dir(), " miXCR clones with PID read cutoff 1.tsv")] %>%
    map(read.delim, stringsAsFactors = F)
  
  names(raw_exp) <- names_exp
  
  exp <- bind_rows(raw_exp, .id = "exp") %>%
    select(-ends_with("R1"), 
           -ends_with("R2"), 
           -ends_with("R4"), 
           -refPoints, 
           -ends_with("FR3"),
           -ends_with("ments"),
           -minQualCDR3,
           -clonalSequenceQuality,
           -allCHitsWithScore,
           -clonalSequence, 
           -Well,
           -cloneId,
           -Subject.id) %>%
    mutate(CDR3_length = str_length(aaSeqCDR3)) %>%
    separate(allVHitsWithScore, c("v_gene", "potential_v_gene"), sep = ",") %>%
    separate(v_gene, c("v_gene", NA), sep = "\\(") %>%
    tidyr::extract(v_gene, into = c("removeV", "v_gene"), regex = "([:lower:])([:alnum:]+)") %>%
    separate(allJHitsWithScore, c("j_gene", "potential_j_gene"), sep = ",") %>%
    separate(j_gene, c("j_gene", NA), sep = "\\(") %>%
    tidyr::extract(j_gene, into = c("removeJ", "j_gene"), regex = "([:lower:])([:alnum:]+)") %>%
    separate(allDHitsWithScore, c("d_gene", "potential_d_gene"), sep = ",") %>%
    separate(d_gene, c("d_gene", NA), sep = "\\(") %>%
    tidyr::extract(d_gene, into = c("removeD", "d_gene"), regex = "([:lower:])([:alnum:]+)") %>%
    select(-starts_with("remove"),
           -starts_with("potential"))
  
  reject_vector <- str_detect(exp$aaSeqCDR3, "_") | 
    str_detect(exp$aaSeqCDR3, "\\*") | 
    str_length(exp$aaSeqCDR3) > 20 | 
    str_length(exp$aaSeqCDR3) < 8 |
    exp$cloneCount < 3 |
    exp$PID.count < 3 #|
  #  nSeqCDR3 != "TGTGCCAGCGGTGAGACAGGGACCAACGAAAGATTATTTTTC"
  
  exp_rejects <- exp %>%
    filter(reject_vector)
  
  setwd(paste(getwd(), "/../.."))
  
  write.csv(exp_rejects, "rejected_CDR3s.csv")
  
  remove(exp_rejects)
  
  exp_clean <- exp %>%
    filter(!reject_vector)
  
  write.csv(exp_clean, "cleaned_CDR3s.csv")
}


#####
# richness eveness
#####

# summary_TCRseq returns a csv file with summary information from the clean_data function

summary_TCRseq <- function(data) {
  
  entropy <- function(x) {
    H <- vector()
    for (i in 1:length(x)) {
      H[i] <- -x[i]*log(x[i])
    }
    H_norm <- sum(H)/log(length(x))
    H_norm
  }
  
  simpsons_index <- function(x) {
    lambda <- vector()
    for (i in 1:length(x)) {
      lambda[i] <- x[i]**2
    }
    lambda <- sum(lambda)
    lambda
  }
  
  
  exp_richness_summary <- data %>%
    group_by(exp, aaSeqCDR3) %>%
    summarise(PID.count = sum(PID.count), PID.fraction = sum(PID.fraction)) %>%
    group_by(exp) %>%
    summarise(richness = n(), 
              evenness = sum(PID.count), 
              average_count = evenness/richness,
              diversity = entropy(PID.fraction),
              simpsons = simpsons_index(PID.fraction))
  
  write.csv(exp_richness_summary, "summary_stats.csv")
}


#####
# only single residues
#####

# data_aa reduces the nucleotide sequences that cause more than aa sequence down to one aa sequence for each repertoire
# it also returns a new column that tells how many nucleotide sequences caused that residue sequence

data_aa <- function(data) {
  data %>%
    group_by(aaSeqCDR3, exp) %>%
    summarise(PID.count_sum = sum(PID.count),
              PID.fraction_sum = sum(PID.fraction),
              n_nuc = n_distinct(nSeqCDR3))
}


#####
# phylgenic tree
#####


exp_clean_only_aa <- exp_clean %>%
  group_by(aaSeqCDR3, tpmouse, response, timepoint, mouse) %>%
  summarise(PID.count_sum = sum(PID.count),
            PID.fraction_sum = sum(PID.fraction),
            n_nuc = n_distinct(nSeqCDR3)) %>%
  mutate(CDR3_length = str_length(aaSeqCDR3)) %>%
  arrange(tpmouse)

hamming_tree_Df <- exp_clean_only_aa %>%
  filter(timepoint == "0", response == "RS", CDR3_length == 14) %>%
  group_by(aaSeqCDR3) %>%
  summarise(n_compartments = n_distinct(tpmouse))

hamming_tree <- stringdistmatrix(hamming_tree_Df$aaSeqCDR3, method = "hamming")

hamming_clust <- hclust(hamming_tree)

hamming_phylo <- as.phylo(hamming_clust)

plot(unroot(hamming_phylo), type = "unrooted", no.margin = T, lab4ut = "axial")

hamming_dend <- as.dendrogram(hamming_clust)

plot(hamming_clust)

hamming_gg <- dendro_data(hamming_dend)

ggdendrogram(hamming_clust) +
  theme_void() +
  scale_y_reverse() +
  coord_polar(theta = "x")

#####
# generate heat map and dendrograms with imputed values
#####

hd_timepoint <- function(df, tp) { 
  x_df <- df %>%
    filter(timepoint == tp) %>%
    group_by(aaSeqCDR3) %>%
    filter(n() >= 3) %>% # change value here to include more/less represented genes
    mutate(Z_fraction = scale(PID.fraction_sum)) %>%
    select(-PID.count_sum, -PID.fraction_sum, -n_nuc, -response, -timepoint, -mouse) %>%
    spread(aaSeqCDR3, Z_fraction)
  
  x_df_ncp <- estim_ncpPCA(x_df[,2:length(x_df)], ncp.max = (length(x_df) - 1), ncp.min = 2)
  
  x_df_complete <- imputePCA(X = as.data.frame(x_df[,2:length(x_df)]), ncp = x_df_ncp$ncp, scale = T)[[1]] %>%
    as.data.frame()
  
  x_df_complete$tpmouse <- x_df$tpmouse
  
  x_df_complete <- x_df_complete %>%
    gather("aaSeqCDR3", "Z_fraction", -tpmouse)
  
  x1 <- x_df_complete %>%
    spread(aaSeqCDR3, Z_fraction) %>%
    column_to_rownames("tpmouse") %>%
    dist() %>%
    hclust(method = "complete") %>%
    as.dendrogram()
  
  y1 <- x_df_complete %>%
    spread(tpmouse, Z_fraction) %>%
    column_to_rownames("aaSeqCDR3") %>%
    dist() %>%
    hclust(method = "complete") %>%
    as.dendrogram()
  
  x_dend <- ggplotGrob(x1 %>%
                         ggdendrogram() +
                         theme_void())
  
  x_df_hm <- df %>%
    ungroup() %>%
    select(-aaSeqCDR3, -timepoint:-n_nuc) %>%
    right_join(x_df_complete) %>%
    mutate(tpmouse = factor(tpmouse, levels = dendro_data(x1)$labels$label),
           aaSeqCDR3 = factor(aaSeqCDR3, levels = dendro_data(y1)$labels$label)) %>%
    distinct()
  
  x_df_hm$dummy <- 1
  
  x_rs <- ggplotGrob(x_df_hm %>%
                       select(tpmouse, response, dummy) %>%
                       distinct() %>%
                       ggplot(aes(tpmouse, dummy)) +
                       geom_tile(aes(fill = response)) +
                       theme_void() +
                       scale_fill_discrete(guide = F))
  
  x_gg <- ggplotGrob(x_df_hm %>%
                       ggplot(aes(tpmouse, aaSeqCDR3)) +
                       geom_tile(aes(fill = Z_fraction)) +
                       theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
                       scale_fill_gradient2(name = "Z score"))
  
  panel_id <- x_gg$layout[x_gg$layout$name == "panel", c("t", "l", "b")]
  
  x_gg <- gtable_add_rows(x_gg, heights = unit(c(1, 0.3), "in"), 0)
  
  x_gg <- gtable_add_grob(x_gg, grobs = list(x_dend, x_rs), t = c(1, 2), l = c(panel_id$l, panel_id$l), b = c(1, panel_id$b))
  
  grid.newpage()
  
  grid.draw(x_gg)
  
}

heatmap_dendrogram(exp_clean_only_aa, "6")


x_df <- exp_clean_only_aa %>%
  filter(timepoint == "6") %>%
  group_by(aaSeqCDR3) %>%
  filter(n() >= 3) %>% 
  mutate(Z_fraction = scale(PID.fraction_sum)) %>%
  select(-PID.count_sum, -PID.fraction_sum, -n_nuc, -response, -timepoint, -mouse) %>%
  spread(aaSeqCDR3, Z_fraction) %>%
  gather("aaSeqCDR3", "Z_fraction", -tpmouse)

x_df[is.na(x_df)] <- 0

x_df_ncp <- estim_ncpPCA(x_df[,2:length(x_df)], ncp.max = (length(x_df) - 1), ncp.min = 2)

x_df_complete <- imputePCA(X = as.data.frame(x_df[,2:length(x_df)]), ncp = x_df_ncp$ncp, scale = T)[[1]] %>%
  as.data.frame()

x_df_complete$tpmouse <- x_df$tpmouse

x_df_complete <- x_df_complete %>%
  gather("aaSeqCDR3", "Z_fraction", -tpmouse)

x1 <- x_df %>%
  column_to_rownames("tpmouse") %>%
  dist() %>%
  hclust(method = "complete") %>%
  as.dendrogram()

y1 <- x_df %>%
  gather("aaSeqCDR3", "Z_fraction", -tpmouse) %>%
  spread(tpmouse, Z_fraction) %>%
  column_to_rownames("aaSeqCDR3") %>%
  dist() %>%
  hclust(method = "complete") %>%
  as.dendrogram()

x_dend <- ggplotGrob(x1 %>%
                       ggdendrogram() +
                       theme_void())

x_df_hm <- exp_clean_only_aa %>%
  ungroup() %>%
  select(-aaSeqCDR3, -timepoint:-n_nuc) %>%
  right_join(x_df) %>%
  mutate(tpmouse = factor(tpmouse, levels = dendro_data(x1)$labels$label),
         aaSeqCDR3 = factor(aaSeqCDR3, levels = dendro_data(y1)$labels$label)) %>%
  distinct()

x_df_hm$dummy <- 1

x_rs <- ggplotGrob(x_df_hm %>%
  select(tpmouse, response, dummy) %>%
  distinct() %>%
  ggplot(aes(tpmouse, dummy)) +
  geom_tile(aes(fill = response)) +
  theme_void() +
    scale_fill_discrete(guide = F))

x_gg <- ggplotGrob(x_df_hm %>%
                     ggplot(aes(tpmouse, aaSeqCDR3)) +
                     geom_tile(aes(fill = Z_fraction)) +
                     theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
                     scale_fill_gradient2(name = "Z score"))

panel_id <- x_gg$layout[x_gg$layout$name == "panel", c("t", "l", "b")]

x_gg <- gtable_add_rows(x_gg, heights = unit(c(1, 0.3), "in"), 0)

x_gg <- gtable_add_grob(x_gg, grobs = list(x_dend, x_rs), t = c(1, 2), l = c(panel_id$l, panel_id$l), b = c(1, panel_id$b))

grid.newpage()

grid.draw(x_gg)

remove(x_gg, panel_id, x_rs, x_df_hm, x_dend, y1, x1, x_df_ncp, x_df_complete, x_df)

#####
# samples sent for TCRseq
#####

sent <- read_xlsx("/Users/joelkidman/Documents/PhD_documents/RNAseq/RENCA_RNAseq_sequencing_primer_order.xlsx", sheet = 2) %>%
  gather(key = "sent", value = "tpmouseresponse", `first run`:`fifth run`) %>%
  select(tpmouseresponse) %>%
  separate(tpmouseresponse, into = c("tpmouse", "response"), sep = " ") %>%
  separate(tpmouse, into = c("timepoint", "mouse"), sep = 1)

table(sent$timepoint, sent$response)

#####
# calculate DE50
######
exp_group <- exp_clean_only_aa %>%
  filter(timepoint == "0") %>%
  group_split(tpmouse)

names(exp_group) <- names_exp

DE50 <- list()

for (j in 1:length(exp_group)) {
  i = 1
  while (sum(exp_group[[j]]$PID.count[1:i]) < sum(exp_group[[j]]$PID.count)/2) {
    i = i + 1
  }
  DE50[[j]] <- i * 100 / nrow(exp_group[[j]])
}
#####
# t-SNE
#####

exp_clean_aa_VJ <- exp_clean %>%
  unite("VJ", v_gene, j_gene, sep = "_", remove = F) %>%
  group_by(VJ) %>%
  group_by(aaSeqCDR3, tpmouse, response, timepoint, mouse, VJ, v_gene, j_gene) %>%
  summarise(PID.count_sum = log(sum(PID.count)),
            PID.fraction_sum = log(sum(PID.fraction)),
            n_nuc = n_distinct(nSeqCDR3),
            Z_score = scale(PID.fraction_sum)) %>%
  mutate(CDR3_length = as.character(str_length(aaSeqCDR3))) %>%
  arrange(tpmouse)

write.csv(exp_clean_aa_VJ, "pyEdgelist_import.csv")
setdiff <- stringdistmatrix(exp_clean_aa_VJ$aaSeqCDR3, method = "lv")

net_tsne <- Rtsne(setdiff, is_distance = T)

remove(setdiff)

tsne_plot_0 <- data.frame(tsne_x = net_tsne$Y[,1], 
                        tsne_y = net_tsne$Y[,2], 
                        mouse = exp_clean_aa_VJ$tpmouse,
                        response = exp_clean_aa_VJ$response,
                        CDR3 = exp_clean_aa_VJ$aaSeqCDR3,
                        CDR3_length = exp_clean_aa_VJ$CDR3_length,
                        VJ_gene = exp_clean_aa_VJ$VJ,
                        v_gene = exp_clean_aa_VJ$v_gene,
                        j_gene = exp_clean_aa_VJ$j_gene,
                        fraction = exp_clean_aa_VJ$PID.fraction_sum,
                        count = exp_clean_aa_VJ$PID.count_sum)

ggplot(tsne_plot_0, aes(tsne_x, tsne_y, colour = count)) +
  geom_point() +
  scale_x_continuous(breaks = seq(from = -50, to = 55, by = 5)) +
  scale_y_continuous(breaks = seq(from = -50, to = 50, by = 5)) +
  scale_color_continuous(low = "blue", high = "red") +
  facet_grid(. ~ response) +
  theme(legend.position = "none")

ggplot(tsne_plot_0, aes(tsne_x, tsne_y)) +
  geom_hex(bins = 50) +
  facet_grid(. ~ response)

tsne_plot %>%
  dplyr::filter(tsne_x > -40, tsne_x < -30, tsne_y > 17, tsne_y < 25) %>%
  select(CDR3)
#####
# CDR3 networks
#####

exp_clean_only_aa <- exp_clean %>%
  select(-nSeqCDR3) %>%
  group_by(aaSeqCDR3, tpmouse, response, timepoint, mouse) %>%
  summarise(PID.count_sum = sum(PID.count),
            PID.fraction_sum = sum(PID.fraction)) %>%
  arrange(tpmouse) %>%
  ungroup() %>%
  group_by(aaSeqCDR3) %>%
  mutate(n_mice = n()) %>%
  group_by(tpmouse) %>%
  arrange(desc(PID.fraction_sum)) %>%
  group_split()

big_rep <- vector()

for (i in 1:length(exp_clean_only_aa)) {
  big_rep[i] <- if_else(length(exp_clean_only_aa[[i]]$aaSeqCDR3) >= 1000, 
                                    T,
                                    F)
}

for (i in 1:length(exp_clean_only_aa[big_rep])) {
  exp_clean_only_aa[big_rep][[i]] <- exp_clean_only_aa[big_rep][[i]][1:1000,]
}

mouse_CDR3_network <- function(df) {
  setdiff_temp <- df %>%
    arrange(n_mice) %>%
    mutate(n_mice = as.character(n_mice))
  
  setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
    as_tibble()
  
  colnames(setdiff) <- setdiff_temp$aaSeqCDR3
  
  setdiff$aa.x <- setdiff_temp$aaSeqCDR3
  
  graph <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select(-starts_with("temp"), -n_mice, -lv) %>%
    graph_from_data_frame(directed = F, vertices = setdiff_temp)
  
  cluster_values <- clusters(graph)$csize %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    arrange(desc(`.`)) %>%
    mutate(rowname_char = as.character(rowname)) %>%
    select(rowname, rowname_char, lay_x = '.') %>%
    filter(lay_x > 2)
  
  cluster_values$rowname <- as.numeric(cluster_values$rowname)
  
  tree_lay <- treemap(cluster_values, 
                      index = c("rowname_char"), 
                      vSize = c("lay_x"),
                      algorithm = "squarified",
                      draw = F)$tm %>%
    arrange(desc(vSize))
  
  graph_components <- decompose.graph(graph)
  
  graph_components_short <- graph_components[cluster_values$rowname]
  
  names(graph_components_short) <- cluster_values$rowname_char
  
  warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
  mypalette <- colorRampPalette(warm)(length(unique(setdiff_temp$n_mice)))
  
  names(mypalette) <- unique(setdiff_temp$n_mice)
  
  mypalette[1] <- "black"
  
  my_pal <- list()
  
  for (col in 1:length(cluster_values$rowname)) {
    my_pal[[col]] <- mypalette[V(graph_components_short[[cluster_values$rowname_char[col]]])$n_mice]
  }
  
  sub_graphing <- function(g, pal) {
    ggplotGrob(g %>%
                 ggnet2(node.size = "PID.fraction_sum",
                        color = "n_mice",
                        palette = pal,
                        color.palette = pal) +
                 theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                       legend.position = "none",
                       plot.background = element_rect(colour = "black")))
  }
  
  graph_graphs <- map2(graph_components_short, my_pal, sub_graphing)
  
  names(graph_graphs) <- cluster_values$rowname_char
  
  annotation_list <- list(grob = graph_graphs,
                          xmin = tree_lay$x0,
                          xmax = tree_lay$w + tree_lay$x0,
                          ymin = tree_lay$y0,
                          ymax = tree_lay$h + tree_lay$y0)
  
  annotation_list <- pmap(annotation_list, annotation_custom)
  
  grid.arrange(ggplot() + annotation_list + labs(title = paste(c("mouse ", unique(setdiff_temp$tpmouse), " ", unique(setdiff_temp$response)), collapse = "")))
}


mouse_CDR3_network(exp_clean_only_aa[[30]])

setwd("D://RNAseq//mouse_networks")

for (i in 1:length(exp_clean_only_aa)) {
  
  g_titles <- exp_clean_only_aa[[i]] %>%
    select(tpmouse, response, timepoint) %>%
    distinct()
  
  png(filename = paste("Mouse ", g_titles$tpmouse, ", ",
                       "Response ", g_titles$response, ", ",
                       "timepoint ", g_titles$timepoint, ".png", sep = ""),
      width = 1527, height = 883)
  
  mouse_CDR3_network(exp_clean_only_aa[[i]])
  
  dev.off()
}

ggsave("mouse_420_RS.png")

#####

setdiff_temp <- exp_clean_only_aa[[1]] %>%
  arrange(n_mice) %>%
  mutate(n_mice = as.character(n_mice),
         PID.fraction_sum = exp(PID.fraction_sum))

setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
  as_tibble()

colnames(setdiff) <- setdiff_temp$aaSeqCDR3

setdiff$aa.x <- setdiff_temp$aaSeqCDR3

graph <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select(-starts_with("temp"), -n_mice, -lv) %>%
  graph_from_data_frame(directed = F, vertices = setdiff_temp)

cluster_values <- clusters(graph)$csize %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(desc(`.`)) %>%
  mutate(rowname_char = as.character(rowname)) %>%
  select(rowname, rowname_char, lay_x = '.') %>%
  filter(lay_x > 2)

cluster_values$rowname <- as.numeric(cluster_values$rowname)

tree_lay <- treemap(cluster_values, 
                    index = c("rowname_char"), 
                    vSize = c("lay_x"),
                    algorithm = "squarified",
                    draw = F)$tm %>%
  arrange(desc(vSize))

graph_components <- decompose.graph(graph)

graph_components_short <- graph_components[cluster_values$rowname]

names(graph_components_short) <- cluster_values$rowname_char

warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
mypalette <- colorRampPalette(warm)(length(unique(setdiff_temp$n_mice)))

names(mypalette) <- unique(setdiff_temp$n_mice)

my_pal <- list()

for (col in 1:length(cluster_values$rowname)) {
  my_pal[[col]] <- mypalette[V(graph_components_short[[cluster_values$rowname_char[col]]])$n_mice]
}

sub_graphing <- function(g, pal) {
  ggplotGrob(g %>%
               ggnet2(node.size = "PID.fraction_sum",
                      color = "n_mice",
                      palette = pal,
                      color.palette = pal) +
               theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                     legend.position = "none", 
                     plot.background = element_rect(fill = "black")))
}

graph_components_short[[1]] %>%
  ggnet2(node.size = "PID.fraction_sum",
         color = "n_mice",
         palette = my_pal[[1]],
         color.palette = my_pal[[1]]) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        legend.position = "none", 
        plot.background = element_rect(fill = "black"))

graph_graphs <- map2(graph_components_short, my_pal, sub_graphing)

names(graph_graphs) <- cluster_values$rowname_char

graph_graphs_order <- graph_graphs[cluster_values$rowname_char]

annotation_list <- list(grob = graph_graphs,
                        xmin = tree_lay$x0,
                        xmax = tree_lay$w + tree_lay$x0,
                        ymin = tree_lay$y0,
                        ymax = tree_lay$h + tree_lay$y0)

annotation_list <- pmap(annotation_list, annotation_custom)

grid.arrange(ggplot() + annotation_list + labs(title = paste(c("mouse ", "003", " ", "RS"), collapse = "")))

#####

CDR3_network <- function(df, tp, rs) {
  setdiff_temp <- df %>%
    filter(timepoint == tp, response == rs) %>%
    group_by(aaSeqCDR3) %>%
    summarise(n_compartments = n_distinct(tpmouse))
  
  setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
    as_tibble()
  
  colnames(setdiff) <- setdiff_temp$aaSeqCDR3
  
  setdiff$aa.x <- setdiff_temp$aaSeqCDR3
  
  vert.x <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select("aa" = aa.x, n_compartments)
  
  vert.y <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select("aa" = aa.y, n_compartments)
  
  vert <- bind_rows(vert.x, vert.y) %>%
    group_by(aa) %>%
    summarise(n_compartments = max(n_compartments))
  
  graph <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select(-starts_with("temp"), -n_compartments, -lv) %>%
    graph_from_data_frame(directed = F, vertices = vert)
  
  cluster_values <- clusters(graph)$csize %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    arrange(desc(`.`)) %>%
    mutate(rowname_char = as.character(rowname)) %>%
    select(rowname, rowname_char, lay_x = '.') %>%
    filter(lay_x > 3)
  
  cluster_values$rowname <- as.numeric(cluster_values$rowname)
  
  tree_lay <- treemap(cluster_values, 
                      index = c("rowname_char"), 
                      vSize = c("lay_x"),
                      algorithm = "squarified",
                      draw = F)$tm %>%
    arrange(desc(vSize))
  
  graph_components <- decompose.graph(graph)
  
  graph_components_short <- graph_components[cluster_values$rowname]
  
  names(graph_components_short) <- cluster_values$rowname_char
  
  warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
  mypalette <- colorRampPalette(warm)(length(unique(setdiff_temp$n_compartments)))
  
  names(mypalette) <- unique(setdiff_temp$n_compartments)
  
  my_pal <- list()
  
  for (col in 1:length(cluster_values$rowname)) {
    my_pal[[col]] <- mypalette[V(graph_components_short[[cluster_values$rowname_char[col]]])$n_compartments]
  }
  
  sub_graphing <- function(g, pal) {
    ggplotGrob(g %>%
                 ggnet2(node.size = 2,
                        color = "n_compartments",
                        palette = pal,
                        color.palette = pal) +
                 theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                       legend.position = "none", 
                       plot.background = element_rect(fill = "black")))
  }
  
  graph_graphs <- map2(graph_components_short, my_pal, sub_graphing)
  
  names(graph_graphs) <- 1:length(graph_graphs)
  
  graph_graphs_order <- graph_graphs[cluster_values$rowname_char]
  
  annotation_list <- list(grob = graph_graphs,
                          xmin = tree_lay$x0,
                          xmax = tree_lay$w + tree_lay$x0,
                          ymin = tree_lay$y0,
                          ymax = tree_lay$h + tree_lay$y0)
  
  annotation_list <- pmap(annotation_list, annotation_custom)
  
  grid.arrange(ggplot() + annotation_list + labs(title = paste(c("timepoint ", tp, " ", rs), collapse = "")))
}

CDR3_network(exp_clean_only_aa, "6", "NR") # example input to function

# function breakdown for trouble shooting

setdiff_temp <- exp_clean_only_aa %>%
  filter(timepoint == "0", response == "NR") %>%
  group_by(aaSeqCDR3) %>%
  summarise(n_compartments = n_distinct(tpmouse))

setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
  as_tibble()

colnames(setdiff) <- setdiff_temp$aaSeqCDR3

setdiff$aa.x <- setdiff_temp$aaSeqCDR3

vert.x <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select("aa" = aa.x, n_compartments)

vert.y <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select("aa" = aa.y, n_compartments)

vert <- bind_rows(vert.x, vert.y) %>%
  group_by(aa) %>%
  summarise(n_compartments = max(n_compartments))

graph <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select(-starts_with("temp"), -n_compartments, -lv) %>%
  graph_from_data_frame(directed = F, vertices = vert)

cluster_values <- clusters(graph)$csize %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(desc(`.`)) %>%
  mutate(rowname_char = as.character(rowname)) %>%
  select(rowname, rowname_char, lay_x = '.') %>%
  filter(lay_x > 3)

cluster_values$rowname <- as.numeric(cluster_values$rowname)

tree_lay <- treemap(cluster_values, 
                    index = c("rowname_char"), 
                    vSize = c("lay_x"),
                    algorithm = "squarified",
                    draw = F)$tm %>%
  arrange(desc(vSize))

graph_components <- decompose.graph(graph)

graph_components_short <- graph_components[cluster_values$rowname]

names(graph_components_short) <- cluster_values$rowname_char

warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
mypalette <- colorRampPalette(warm)(length(unique(setdiff_temp$n_compartments)))

names(mypalette) <- unique(setdiff_temp$n_compartments)

my_pal <- list()

for (col in 1:length(cluster_values$rowname)) {
  my_pal[[col]] <- mypalette[V(graph_components_short[[cluster_values$rowname_char[col]]])$n_compartments]
}

sub_graphing <- function(g, pal) {
  ggplotGrob(g %>%
               ggnet2(node.size = 2,
                      color = "n_compartments",
                      palette = pal,
                      color.palette = pal) +
               theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                     legend.position = "none"))
}

graph_graphs <- map2(graph_components_short, my_pal, sub_graphing)

names(graph_graphs) <- 1:length(graph_graphs)


annotation_list <- list(grob = graph_graphs,
                        xmin = tree_lay$x0,
                        xmax = tree_lay$w + tree_lay$x0,
                        ymin = tree_lay$y0,
                        ymax = tree_lay$h + tree_lay$y0)

annotation_list <- pmap(annotation_list, annotation_custom)

grid.arrange(ggplot() + annotation_list + labs(title = paste(c("timepoint ", "0", " ", "NR"), collapse = "")))
#####
# CDR3 network graph statistics
#####

graph_func <- function(df, tp, rs) {
  setdiff_temp <- exp_clean_only_aa %>%
    filter(timepoint == tp, response == rs) %>%
    group_by(aaSeqCDR3) %>%
    summarise(n_compartments = n_distinct(tpmouse))
  
  setdiff_specific <- exp_clean_only_aa %>%
    filter(timepoint == tp, response == rs) %>%
    group_by(tpmouse) %>%
    group_split()
  
  setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
    as_tibble()
  
  colnames(setdiff) <- setdiff_temp$aaSeqCDR3
  
  setdiff$aa.x <- setdiff_temp$aaSeqCDR3
  
  vert.x <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select("aaSeqCDR3" = aa.x, n_compartments)
  
  vert.y <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select("aaSeqCDR3" = aa.y, n_compartments)
  
  setdiff_specific[[length(setdiff_specific) + 1]] <- bind_rows(vert.x, vert.y) %>%
    group_by(aaSeqCDR3) %>%
    summarise(n_compartments = max(n_compartments))
  
  vert <- reduce(setdiff_specific, 
                 left_join, 
                 by = "aaSeqCDR3", 
                 .init = setdiff_specific[[length(setdiff_specific)]])
  
  
  graph <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select(-starts_with("temp"), -n_compartments, -lv) %>%
    graph_from_data_frame(directed = F, vertices = vert)
  
}

summary_func <- function(df, tp, rs) {
  setdiff_temp <- df %>%
    filter(timepoint == tp, response == rs) %>%
    group_by(aaSeqCDR3) %>%
    summarise(n_compartments = n_distinct(tpmouse))
  
  setdiff_specific <- df %>%
    filter(timepoint == tp, response == rs) %>%
    group_by(tpmouse) %>%
    group_split()
  
  setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
    as_tibble()
  
  colnames(setdiff) <- setdiff_temp$aaSeqCDR3
  
  setdiff$aa.x <- setdiff_temp$aaSeqCDR3
  
  vert.x <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select("aaSeqCDR3" = aa.x, n_compartments)
  
  vert.y <- setdiff %>%
    gather("aa.y", "lv", -aa.x) %>%
    filter(lv == 1) %>% # set lv distance
    unite("temp1", aa.x, aa.y, remove = F) %>%
    mutate(temp2 = pmin(aa.x, aa.y),
           temp3 = pmax(aa.x, aa.y)) %>%
    unite("temp2", temp2, temp3, sep = "_") %>%
    filter(temp1 != temp2) %>%
    left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
    select("aaSeqCDR3" = aa.y, n_compartments)
  
  setdiff_specific[[length(setdiff_specific) + 1]] <- bind_rows(vert.x, vert.y) %>%
    group_by(aaSeqCDR3) %>%
    summarise(n_compartments = max(n_compartments))
  
  vert <- reduce(setdiff_specific, left_join, by = "aaSeqCDR3", .init = setdiff_specific[[length(setdiff_specific)]])
  
  
  graph_stats <- vert %>%
    gather("label", "tpmouse", starts_with("tpmouse")) %>%
    dplyr::ungroup() %>%
    group_by(label, tpmouse) %>%
    summarise(nodes = n()) %>%
    filter(!is.na(tpmouse))
  
}

tp0NR <- graph_func(exp_clean_only_aa, "0", "NR")

tp0NR_stats <- summary_func(exp_clean_only_aa, "0", "NR")

V(tp0NR)

tp_response_graphs <- list()

tp_response_stats <- list()

tp_vector <- rep(unique(exp_clean_only_aa$timepoint), each = length(unique(exp_clean_only_aa$response)))

timepoints <- unique(exp_clean_only_aa$timepoint)

responses <- unique(exp_clean_only_aa$response)

rs_vector <- rep(unique(exp_clean_only_aa$response), times = length(unique(exp_clean_only_aa$timepoint)))

graphs_list <- list(df = as.list(rep(list(exp_clean_only_aa), each = length(rs_vector))),
     tp = as.list(tp_vector),
     rs = as.list(rs_vector))

tp_response_graphs <- pmap(exp_clean_only_aa, summary_func)

m <- 1

for (rs in 1:length(responses)) {
  for (tp in 1:length(timepoints)) {
    tp_response_graphs[[m]] <- graph_func(exp_clean_only_aa, timepoints[tp], responses[rs])
    
    tp_response_stats[[m]] <- summary_func(exp_clean_only_aa, timepoints[tp], responses[rs])
    
    m <- m + 1
  }
}

#####

setdiff_temp <- exp_clean_only_aa %>%
  filter(timepoint == "0", response == "NR") %>%
  group_by(aaSeqCDR3) %>%
  summarise(n_compartments = n_distinct(tpmouse))

setdiff_specific <- exp_clean_only_aa %>%
  filter(timepoint == "0", response == "NR") %>%
  group_by(tpmouse) %>%
  group_split()

setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
  as_tibble()

colnames(setdiff) <- setdiff_temp$aaSeqCDR3

setdiff$aa.x <- setdiff_temp$aaSeqCDR3

vert.x <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select("aaSeqCDR3" = aa.x, n_compartments)

vert.y <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select("aaSeqCDR3" = aa.y, n_compartments)

setdiff_specific[[length(setdiff_specific) + 1]] <- bind_rows(vert.x, vert.y) %>%
  group_by(aaSeqCDR3) %>%
  summarise(n_compartments = max(n_compartments))

vert <- reduce(setdiff_specific, left_join, by = "aaSeqCDR3", .init = setdiff_specific[[5]])
  

graph <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select(-starts_with("temp"), -n_compartments, -lv) %>%
  graph_from_data_frame(directed = F, vertices = vert)

#####

#####
setdiff_temp <- exp_clean_only_aa %>%
  filter(timepoint == "0", response == "NR") %>%
  group_by(aaSeqCDR3) %>%
  summarise(n_compartments = n_distinct(tpmouse))

setdiff_specific <- exp_clean_only_aa %>%
  filter(timepoint == "0", response == "NR") %>%
  group_by(tpmouse) %>%
  group_split()

setdiff <- stringdistmatrix(setdiff_temp$aaSeqCDR3, setdiff_temp$aaSeqCDR3) %>%
  as_tibble()

colnames(setdiff) <- setdiff_temp$aaSeqCDR3

setdiff$aa.x <- setdiff_temp$aaSeqCDR3

vert.x <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select("aaSeqCDR3" = aa.x, n_compartments)

vert.y <- setdiff %>%
  gather("aa.y", "lv", -aa.x) %>%
  filter(lv == 1) %>% # set lv distance
  unite("temp1", aa.x, aa.y, remove = F) %>%
  mutate(temp2 = pmin(aa.x, aa.y),
         temp3 = pmax(aa.x, aa.y)) %>%
  unite("temp2", temp2, temp3, sep = "_") %>%
  filter(temp1 != temp2) %>%
  left_join(setdiff_temp, by = c("aa.x" = "aaSeqCDR3")) %>%
  select("aaSeqCDR3" = aa.y, n_compartments)

setdiff_specific[[length(setdiff_specific) + 1]] <- bind_rows(vert.x, vert.y) %>%
  group_by(aaSeqCDR3) %>%
  summarise(n_compartments = max(n_compartments))

vert <- reduce(setdiff_specific, left_join, by = "aaSeqCDR3", .init = setdiff_specific[[length(setdiff_specific)]])


graph_stats <- vert %>%
  gather("label", "tpmouse", starts_with("tpmouse")) %>%
  dplyr::ungroup() %>%
  group_by(label, tpmouse) %>%
  summarise(nodes = n()) %>%
  filter(!is.na(tpmouse))
#####

mice <- tp0NR_stats$tpmouse

tpmouse <- tp0NR_stats$label

graph_stats_degree <- vector()

m <- 1

for (g in 1:length(tp_response_stats)) {
  tpmouse <- tp_response_stats[[g]]$label
  for (h in 1:length(tp_response_stats[[g]]$tpmouse)) {
    graph_stats_degree[m] <- igraph::V(tp_response_graphs[[g]])[!is.na(vertex_attr(tp_response_graphs[[g]], name = tpmouse[h]))][degree(tp_response_graphs[[g]], V(tp_response_graphs[[g]])[!is.na(vertex_attr(tp_response_graphs[[g]], name = tpmouse[h]))]) > 1] %>%
      length()
    m <- m + 1
  }
}

tpmouse_response <- exp_clean_only_aa %>%
  ungroup() %>%
  select(tpmouse, response, timepoint) %>%
  distinct()

exp_stats <- tp_response_stats %>%
  bind_rows() %>%
  left_join(tpmouse_response)

exp_stats$degree <- graph_stats_degree

exp_stats <- exp_stats %>%
  mutate(degree_norm = degree/nodes)

ggplot(exp_stats, aes(factor(timepoint), degree_norm, fill = response)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(aes(group = response), label = "p.signif") +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  labs(y = "connected nodes", x = "timepoint") +
  ylim(c(min(exp_stats$degree_norm), 0.75))

ggplot(exp_stats, aes(degree_norm, nodes, colour = response)) +
  geom_point() +
  facet_grid(. ~ timepoint)
#####
#####
# tissue network
#####


tissue_network <- function(df, pop) {
  
  exp_intersect <- list()
  
  for (i in 1:length(df)) {
    exp_intersect[[i]] <- lapply(df, inner_join, y = df[[i]], by = "aaSeqCDR3")
  }
  
  names(exp_intersect) <- unique(exp_clean$exp)
  
  for (i in 1:length(exp_intersect)) {
    names(exp_intersect[[i]]) <- unique(exp_clean$exp)
  }
  
  exp_intersect_unique <- exp_intersect %>%
    map(bind_rows) %>%
    bind_rows() %>%
    filter(population.x == population.y, population.x == pop) %>%
    select(flank_mouse_pop.x, flank_mouse_pop.y) %>%
    filter(flank_mouse_pop.x != flank_mouse_pop.y) %>%
    mutate(temp_min = pmin(flank_mouse_pop.x, flank_mouse_pop.y),
           temp_max = pmax(flank_mouse_pop.x, flank_mouse_pop.y)) %>%
    unite("Intersect", temp_min, temp_max) %>%
    select(Intersect) %>%
    unique()
  
  exp_intersect_clean <- exp_intersect %>%
    map(bind_rows) %>%
    bind_rows() %>%
    unite("Intersect", flank_mouse_pop.x, flank_mouse_pop.y, sep = "_", remove = F) %>%
    filter(mouse_flank.x != mouse_flank.y, population.x == population.y) %>%
    arrange(desc(PID.count.x)) %>%
    right_join(exp_intersect_unique) %>%
    mutate(residual = -log(PID.fraction.x) - -log(PID.fraction.y), popul = population.x, norm_residual = scale(residual))
  
  exp_intersect_summary <- exp_intersect_clean %>%
    group_by(Intersect, popul) %>%
    summarise(n_clones = n(), 
              residual = 1/(sum(residual^2)/(n_clones - 1))) %>%
    group_by(popul) %>%
    mutate(norm_n_clones = scale(n_clones),
           norm_residual = scale(residual),
           mm_residual = (((residual - min(residual)) / (max(residual) - min(residual)) + 0.1) * 10),
           log_residual = -log(residual),
           scale_log_residual = scale(log_residual)) %>%
    separate(Intersect, into = c(NA, "mouse_x", NA, NA, "mouse_y", NA), remove = F, sep = "_") %>%
    mutate(mouse_x = as.numeric(mouse_x), mouse_y = as.numeric(mouse_y))
  
  exp_edgelist <- exp_intersect_unique %>%
    left_join(exp_intersect_summary[,c(1,9)]) %>%
    separate(Intersect, into = c("from", "to"), sep = 10) %>%
    mutate(from = str_trunc(from, 9, side = "right", ellipsis = "")) %>%
    filter(mm_residual > 5)
  
  exp_vertexlist <- exp_intersect_unique %>%
    separate(Intersect, into = c("from", "to"), sep = 10) %>%
    mutate(from = str_trunc(from, width = 9, side = "right", ellipsis = "")) %>%
    full_join(., .[,2], by = c("from" = "to")) %>%
    select("vertex" = "from", -to) %>%
    distinct() %>%
    separate(vertex, into = c("flank", "mouse", "pop"), remove = F, sep = "_")
  
  exp_network <- graph_from_data_frame(exp_edgelist, directed = F, vertices = exp_vertexlist)
  
  ggnet2(exp_network, 
         edge.size = E(exp_network)$mm_residual,
         color = V(exp_network)$mouse,
         shape = V(exp_network)$flank,
         color.palette = "Dark2",
         size = 20) +
    labs(title = pop)
}

exp_group <- exp_clean_only_aa %>%
  filter(timepoint == "0") %>%
  group_by(tpmouse) %>%
  group_split()

exp_intersect_summary <- exp_intersect %>%
  group_by(Intersect, popul) %>%
  summarise(n_clones = n(), 
            residual = sum(residual^2)/(n_clones - 1)) %>%
  group_by(popul) %>%
  mutate(norm_n_clones = scale(n_clones),
         norm_residual = scale(residual),
         mm_residual = (((residual - min(residual)) / (max(residual) - min(residual)) + 0.1) * 10),
         log_residual = -log(residual),
         scale_log_residual = scale(log_residual)) %>%
  separate(Intersect, into = c(NA, "mouse_x", NA, NA, "mouse_y", NA), remove = F, sep = "_") %>%
  mutate(mouse_x = as.numeric(mouse_x), mouse_y = as.numeric(mouse_y))

exp_edgelist <- exp_intersect_unique %>%
  left_join(exp_intersect_summary[,c(1,9)]) %>%
  separate(Intersect, into = c("from", "to"), sep = 10) %>%
  mutate(from = str_trunc(from, 9, side = "right", ellipsis = "")) %>%
  filter(mm_residual > 5)

exp_vertexlist <- exp_intersect_unique %>%
  separate(Intersect, into = c("from", "to"), sep = 10) %>%
  mutate(from = str_trunc(from, width = 9, side = "right", ellipsis = "")) %>%
  full_join(., .[,2], by = c("from" = "to")) %>%
  select("vertex" = "from", -to) %>%
  distinct() %>%
  separate(vertex, into = c("flank", "mouse", "pop"), remove = F, sep = "_")

exp_network <- graph_from_data_frame(exp_edgelist, directed = F, vertices = exp_vertexlist)

ggnet2(exp_network, 
       edge.size = E(exp_network)$mm_residual,
       color = V(exp_network)$mouse,
       shape = V(exp_network)$flank,
       color.palette = "Dark2",
       size = 20) +
  labs(title = pop)

tissue_network(exp_group, "CD8")

tissue_network(exp_group, "CD4")
#####
# overlapping clones
#####

exp_group <- exp_clean_only_aa %>%
  filter(timepoint == "0") %>%
  group_by(tpmouse) %>%
  group_split()

# defined vectors and lists

tp_vector <- exp_clean_only_aa$timepoint %>% unique()

rs_vector <- exp_clean_only_aa$response %>% unique()

combinations = 3

exp_calc <- CombSet(exp_group, m = combinations)

map_list <- list()

m <- 1

sum(CombN(1:8, m = 2:8)) + sum(CombN(1:4, m = 2:4)) + sum(CombN(1:4, m = 2:4)) + sum(CombN(1:8, m = 2:8)) + sum(CombN(1:7, m = 2:7)) + sum(CombN(1:8, m = 2:8)) + sum(CombN(1:8, m = 2:8)) + sum(CombN(1:7, m = 2:7)) + sum(CombN(1:5, m = 2:5))

# loops

for (i in 1:length(rs_vector)) {
  for (j in 1:length(tp_vector)) {
    exp_group <- exp_clean_only_aa %>%
      filter(timepoint == tp_vector[j], response == rs_vector[i]) %>%
      group_by(tpmouse) %>%
      group_split()
    for (k in 2:length(exp_group)) {
      exp_calc <- CombSet(exp_group, m = k)
      for (l in 1:(length(exp_calc)/k)) {
        map_list[[m]] <- exp_calc[l, ] %>% reduce(inner_join, by = "aaSeqCDR3")
        m <- m + 1
      }
    }
  }
}

test_rs <- 1:2

test_tp <- 1:4

test_eg <- 2:5

test_ec <- CombN(1:length(test_eg), m = 2:length(test_eg))

for (i in 1:length(rs_vector)) {
  for (j in 1:length(tp_vector)) {
    for (k in 2:length(exp_group)) {
      for (l in 1:(sum(test_ec)/k)) {
        print(paste(i, j, k - 1, l))
      }
    }
  }
}

exp_intersect <- map_list %>%
  bind_rows() %>%
  unite("mice", starts_with("tpmouse"), sep = "_", remove = F) 

exp_int <- exp_intersect %>%
  mutate(clone_sum = rowSums(select(., starts_with("PID.fraction_sum")), na.rm = T)) %>%
  group_by(mice, timepoint.x, response.x) %>%
  summarise(n_clones = n(),
            clonal_proportion_norm = sum(clone_sum)/n_clones,
            clonal_proportion = sum(clone_sum)) %>%
  mutate(n_mice = str_count(mice, "[\\d][\\d][\\d]"))

ggplot(exp_int, aes(x = as.factor(n_mice), y = clonal_proportion, fill = response.x)) +
  geom_boxplot() +
  facet_grid(. ~ timepoint.x, scales = "free_x") +
  labs(x = "n mice", y = "shared clonal proportion per shared clone")

ggplot(exp_int, aes(x = as.factor(n_mice), y = n_clones, fill = response.x)) +
  geom_boxplot() +
  facet_grid(. ~ timepoint.x, scales = "free_x") +
  labs(x = "n mice", y = "shared clones") +
  ylim(c(0, 175))

ggplot(exp_int, aes(x = as.factor(n_mice), y = clonal_proportion, fill = response.x)) +
  geom_boxplot() +
  facet_grid(. ~ timepoint.x, scales = "free_x") +
  labs(x = "n mice", y = "clonal proportion")

#####
# Morisita network
#####
morisita_network <- function(df, pop) {
  
  exp_group_only_aa <- exp_clean_only_aa %>%
    arrange(desc(population)) %>%
    unite("flank_mouse_pop", flank_mouse, population, remove = F) %>%
    group_by(flank_mouse_pop) %>%
    group_split()
  
  exp_intersect <- list()
  
  for (i in 1:length(exp_group_only_aa)) {
    exp_intersect[[i]] <- lapply(exp_group_only_aa, full_join, y = exp_group_only_aa[[i]], by = "aaSeqCDR3")
  }
  
  names(exp_intersect) <- unique(exp_clean$flank_mouse_pop)
  
  for (i in 1:length(exp_intersect)) {
    names(exp_intersect[[i]]) <- unique(exp_clean$flank_mouse_pop)
  }
  
  exp_intersect_unique <- exp_intersect %>%
    map(bind_rows) %>%
    bind_rows() %>%
    filter(population.x == population.y, population.x == pop) %>%
    select(flank_mouse_pop.x, flank_mouse_pop.y) %>%
    filter(flank_mouse_pop.x != flank_mouse_pop.y) %>%
    mutate(temp_min = pmin(flank_mouse_pop.x, flank_mouse_pop.y),
           temp_max = pmax(flank_mouse_pop.x, flank_mouse_pop.y)) %>%
    unite("Intersect", temp_min, temp_max) %>%
    select(Intersect) %>%
    unique()
  
  exp_intersect_clean <- exp_intersect %>%
    map(bind_rows) %>%
    bind_rows()
  
  exp_intersect_clean$PID.count_sum.y[is.na(exp_intersect_clean$PID.count_sum.y)] <- 0
  
  exp_intersect_clean$PID.fraction_sum.y[is.na(exp_intersect_clean$PID.fraction_sum.y)] <- 0
  
  exp_intersect_clean$PID.count_sum.x[is.na(exp_intersect_clean$PID.count_sum.x)] <- 0
  
  exp_intersect_clean$PID.fraction_sum.x[is.na(exp_intersect_clean$PID.fraction_sum.x)] <- 0
  
  exp_intersect_clean <- exp_intersect_clean %>%
    fill(flank_mouse_pop.x:population.x, flank_mouse_pop.y, population.y) %>%
    unite("test", flank_mouse_pop.x, flank_mouse_pop.y, sep = "_", remove = F) %>%
    filter(flank_mouse.x != flank_mouse.y, population.x == population.y, population.x == pop) %>%
    arrange(desc(PID.count_sum.x)) %>%
    mutate(temp_min = pmin(flank_mouse_pop.x, flank_mouse_pop.y),
           temp_max = pmax(flank_mouse_pop.x, flank_mouse_pop.y)) %>%
    unite("Intersect", temp_min, temp_max) %>%
    filter(Intersect == test) %>%
    select(PID.count_sum.x, PID.count_sum.y, aaSeqCDR3, Intersect) %>%
    group_split(Intersect)
  
  names_vector <- vector()
  
  for (j in 1:length(exp_intersect_clean)) {
    names_vector[j] <- unique(exp_intersect_clean[[j]]$Intersect)
    names_vector
  }
  
  names(exp_intersect_clean) <- names_vector
  
  for (i in 1:length(exp_intersect_clean)) {
    
    exp_intersect_clean[[i]] <- exp_intersect_clean[[i]] %>%
      column_to_rownames(var = "aaSeqCDR3") %>%
      select(-Intersect) %>%
      as.matrix() %>%
      mh(PlugIn = T) %>%
      bind_rows() %>%
      as.data.frame()
  }
  
  exp_intersect_clean <- bind_cols(exp_intersect_clean)
  
  colnames(exp_intersect_clean) <- names_vector
  
  exp_edgelist <- exp_intersect_clean[2,] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    separate(rowname, into = c("from", "to"), sep = 10) %>%
    mutate(from = str_trunc(from, 9, side = "right", ellipsis = ""),
           morisita = `2` * 10) %>%
    select(-`2`) %>%
    filter(morisita > 7) # remove irrelevant edges
  
  exp_vertexlist <- exp_intersect_unique %>%
    separate(Intersect, into = c("from", "to"), sep = 10) %>%
    mutate(from = str_trunc(from, width = 9, side = "right", ellipsis = "")) %>%
    full_join(., .[,2], by = c("from" = "to")) %>%
    select("vertex" = "from", -to) %>%
    distinct() %>%
    separate(vertex, into = c("flank", "mouse", "pop"), remove = F, sep = "_")
  
  exp_network <- graph_from_data_frame(exp_edgelist, directed = F, vertices = exp_vertexlist)
  
  ggnet2(exp_network, 
         edge.size = E(exp_network)$morisita,
         color = V(exp_network)$mouse,
         shape = V(exp_network)$flank,
         color.palette = "Dark2",
         size = 20) +
    labs(title = pop)
}

morisita_network(exp_clean_only_aa, "CD8")

morisita_network(exp_clean_only_aa, "CD4")
#####
# intersection plot
#####
write.csv(exp_clone_intersect_summary, "clone_intersections.csv")

exp_clone_intersect <- exp_clone_intersect_summary %>%
  group_by(pop) %>%
  mutate(`scaled clone overlap` = scale(intersect_CDR3)) %>%
  separate(from, into = c("tumour1", NA), sep = 5) %>%
  separate(to, into = c("tumour2", NA), sep = 5)

ggplot(exp_clone_intersect, aes(tumour1, tumour2, fill = `scaled clone overlap`)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(low = "blue",
                      high = "red",
                      mid = "black") +
  facet_wrap(~ pop, scales = "free")
#####
# clones overlap
#####

exp_intersect <- list()

for (i in 1:length(exp_group)) {
  exp_intersect[[i]] <- lapply(exp_group, inner_join, y = exp_group[[i]], by = "aaSeqCDR3")
}

names(exp_intersect) <- unique(exp_clean$exp)

for (i in 1:length(exp_intersect)) {
  names(exp_intersect[[i]]) <- unique(exp_clean$exp)
}

exp_intersect_unique <- exp_intersect %>%
  map(bind_rows) %>%
  bind_rows() %>%
  filter(population.x == population.y) %>%
  select(flank_mouse_pop.x, flank_mouse_pop.y) %>%
  filter(flank_mouse_pop.x != flank_mouse_pop.y) %>%
  mutate(temp_min = pmin(flank_mouse_pop.x, flank_mouse_pop.y),
         temp_max = pmax(flank_mouse_pop.x, flank_mouse_pop.y)) %>%
  unite("Intersect", temp_min, temp_max) %>%
  select(Intersect) %>%
  unique()

exp_intersect_clean <- exp_intersect %>%
  map(bind_rows) %>%
  bind_rows() %>%
  unite("Intersect", flank_mouse_pop.x, flank_mouse_pop.y, sep = "_", remove = F) %>%
  filter(mouse_flank.x != mouse_flank.y, population.x == population.y) %>%
  arrange(desc(PID.count.x)) %>%
  right_join(exp_intersect_unique) %>%
  mutate(residual = PID.fraction.x - PID.fraction.y, popul = population.x, norm_residual = scale(residual))

exp_intersect_clean$mouse_flank.x <- factor(exp_intersect_clean$mouse_flank.x, levels = c("1.0_L", "1.2_L", "1.3_L", "1.0_R", "1.2_R"))
exp_intersect_clean$mouse_flank.y <- factor(exp_intersect_clean$mouse_flank.y, levels = c("1.3_R", "1.2_R", "1.0_R", "1.3_L", "1.2_L"))

gg_plot <- exp_intersect_clean %>%
  filter(popul == "CD4") %>%
  ggplot(aes(PID.fraction.x, PID.fraction.y)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0) +
  facet_grid(mouse_flank.y ~ mouse_flank.x) +
  labs(title = "CD4")

gg_grob <- ggplotGrob(gg_plot)

idx <- which(gg_grob$layout$name %in% c("panel-2-5",
                                        "panel-3-5",
                                        "panel-4-3", 
                                        "panel-5-3", 
                                        "panel-3-4", 
                                        "panel-4-4", 
                                        "panel-4-5", 
                                        "panel-5-2", 
                                        "panel-5-3", 
                                        "panel-5-4", 
                                        "panel-5-5"))

for (i in idx) gg_grob$grobs[[i]] <- nullGrob()

grid.newpage()
grid.draw(gg_grob)

#####
# residue network
#####

residue_network <- function(List) {
  
  residue_profile <- List %>%
    as.data.frame() %>%
    mutate(residue_group = str_replace_all(aaSeqCDR3, c("[AIGLPV]" = "J", 
                                                        "[FWY]" = "O", 
                                                        "[DE]" = "U", 
                                                        "[RHK]" = "B",
                                                        "[ST]" = "X",
                                                        "[CM]" = "Y",
                                                        "[NQ]" = "Z"))) %>%
    separate(aaSeqCDR3, into = c(letters[1:max(.$CDR3_length)]), sep = c(1:max(.$CDR3_length)), remove = F) %>%
    separate(residue_group, into = c(LETTERS[1:max(.$CDR3_length)]), sep = c(1:max(.$CDR3_length)), remove = F)
    
    for (j in match("A", names(residue_profile)):match(LETTERS[max(List$CDR3_length)], names(residue_profile))) {
      residue_profile[,j] <- paste(residue_profile[,j], as.character(j - match("A", names(residue_profile)) + 1), sep = "")
    }
    
    raw_residue_vertex_list <- list()
    
    for (k in 1:max(List$CDR3_length)) {
      raw_residue_vertex_list[[k]] <- residue_profile[,k + match("A", names(residue_profile)) - 1] %>%
        as.data.frame(stringsAsFactors = F)
    }
    
    names(raw_residue_vertex_list) <- colnames(residue_profile)[ match("A", names(residue_profile)):match(LETTERS[max(List$CDR3_length)], names(residue_profile))]
    
    residue_vertex_list <- raw_residue_vertex_list %>%
      bind_rows(.id = "id")
    
    names(residue_vertex_list) <- c("id", "residue")
    
    residue_vertex_list <-  residue_vertex_list %>%
      group_by(residue) %>%
      mutate(counts = n()) %>%
      ungroup() %>%
      distinct() %>%
      separate(residue, into = c("colour", "position"), sep = 1, remove = F)
    
    residue_vertex_list <- residue_vertex_list[str_detect(residue_vertex_list$residue, "[:alpha:]+"),]
    
    residue_vertex_list <-  residue_vertex_list %>%
      group_by(position) %>%
      mutate(counts = (((counts - min(counts)) / (max(counts) - min(counts))) + 0.1) * 5) %>%
      select(-id)
    
    residue_vertex_list$colour <-  if_else(residue_vertex_list$colour == "J", "Aliphatic", 
                                           if_else(residue_vertex_list$colour == "O", "Aromatic",
                                                   if_else(residue_vertex_list$colour == "U", "Acidic", 
                                                           if_else(residue_vertex_list$colour == "B", "Basic",
                                                                   if_else(residue_vertex_list$colour == "X", "Hydroxilic",
                                                                           if_else(residue_vertex_list$colour == "Y", "Sulfuric",
                                                                                   "Amidic"))))))
    
    residue_vertex_list$colour <- factor(residue_vertex_list$colour)
    
    residue_vertex_list$counts[is.nan(residue_vertex_list$counts)] <- 5.5
    
    residue_edge_list <- list()
    
    for (l in 1:(max(List$CDR3_length) - 1)) {
      residue_edge_list[[l]] <- table(residue_profile[, l + match("A", names(residue_profile)) - 1], residue_profile[,l + match("A", names(residue_profile))]) %>%
        as.data.frame(stringsAsFactors = F) %>%
        melt() %>%
        filter(value != 0) %>%
        select(-variable)
    }
    
    residue_edge_list <- residue_edge_list %>%
      bind_rows()
    
    residue_edge_list <- residue_edge_list[str_detect(residue_edge_list$Var1, "[A-Z]+"),]
    
    residue_edge_list <- residue_edge_list[str_detect(residue_edge_list$Var2, "[A-Z]+"),]
    
    residue_edge_list <- residue_edge_list %>%
      separate(Var1, into = c("colour", "position"), sep = 1, remove = F) %>%
      group_by(position) %>%
      mutate(weight = (((value - min(value)) / (max(value) - min(value))) + 0.1) * 5) %>%
      ungroup() %>%
      select(Var1, Var2, weight, position, value)
    
    residue_edge_list$weight[is.nan(residue_edge_list$weight)] <- 5.5
    
    residue_edge_list$edge_colour <- if_else(residue_edge_list$weight == 5.5, residue_edge_list$edge_colour <- 2, residue_edge_list$edge_colour <- 1)
    
    residue_edge_list$edge_size <- if_else(residue_edge_list$weight == 5.5, residue_edge_list$edge_size <- 2, residue_edge_list$edge_size <- 0.25)
    
    residue_edge_list <- residue_edge_list %>%
      filter(weight >= 1)
    
    joined_nodes <- c(residue_edge_list$Var1, residue_edge_list$Var2) %>%
      as.data.frame(stringsAsFactors = F) %>%
      distinct()
    
    residue_vertex_list <- residue_vertex_list %>%
      right_join(joined_nodes, by = c("residue" = "."))
    
    graph <- residue_edge_list %>%
      graph_from_data_frame(vertices = residue_vertex_list)
    
    ggnet2(graph, node.size = V(graph)$counts,
           arrow.size = 5,
           node.color = V(graph)$colour,
           color.palette = "Set3",
           edge.size = E(graph)$edge_size,
           edge.color = E(graph)$edge_colour,
           arrow.gap = 0.02) +
      guides(size = "none") +
      labs(title = List[1,1])

}

residue_networks <- list()

setwd("D://mice_Experiments//4_1_0//residue_networks")

for (i in 1:length(exp_group)) {
  
  gg <- residue_network(exp_group[[i]])
  
  ggsave(paste(exp_group[[i]][1,1], ".png", sep = ""), gg)
  
}

residue_network(exp_group[[12]])
#####
# residue test function
#####


exp_group <- exp_clean_only_aa %>%
  mutate(CDR3_length = str_length(aaSeqCDR3)) %>%
  group_by(timepoint, response) %>%
  group_split()

residue_profile <- exp_group[[8]] %>%
  as.data.frame() %>%
  mutate(residue_group = str_replace_all(aaSeqCDR3, c("[AIGLPV]" = "J", 
                                                      "[FWY]" = "O", 
                                                      "[DE]" = "U", 
                                                      "[RHK]" = "B",
                                                      "[ST]" = "X",
                                                      "[CM]" = "Y",
                                                      "[NQ]" = "Z"))) %>%
  separate(aaSeqCDR3, into = c(letters[1:max(.$CDR3_length)]), sep = c(1:max(.$CDR3_length)), remove = F) %>%
  separate(residue_group, into = c(LETTERS[1:max(.$CDR3_length)]), sep = c(1:max(.$CDR3_length)), remove = F)

for (j in match("A", names(residue_profile)):match(LETTERS[max(exp_group[[8]]$CDR3_length)], names(residue_profile))) {
  residue_profile[,j] <- paste(residue_profile[,j], as.character(j - match("A", names(residue_profile)) + 1), sep = "")
}

raw_residue_vertex_list <- rep(list(NULL), times = max(exp_group[[8]]$CDR3_length))

for (k in 1:max(exp_group[[8]]$CDR3_length)) {
  raw_residue_vertex_list[[k]] <- residue_profile[,k + match("A", names(residue_profile)) - 1] %>%
    as.data.frame(stringsAsFactors = F)
}

names(raw_residue_vertex_list) <- colnames(residue_profile)[ match("A", names(residue_profile)):match(LETTERS[max(exp_group[[8]]$CDR3_length)], names(residue_profile))]

residue_vertex_list <- raw_residue_vertex_list %>%
  bind_rows(.id = "id")

names(residue_vertex_list) <- c("id", "residue")

residue_vertex_list <-  residue_vertex_list %>%
  group_by(residue) %>%
  mutate(counts = n()) %>%
  ungroup() %>%
  distinct() %>%
  separate(residue, into = c("colour", "position"), sep = 1, remove = F)

residue_vertex_list <- residue_vertex_list[str_detect(residue_vertex_list$residue, "[:alpha:]+"),]

residue_vertex_list <-  residue_vertex_list %>%
  group_by(position) %>%
  mutate(counts = (((counts - min(counts)) / (max(counts) - min(counts))) + 0.1) * 5) %>%
  select(-id)

residue_vertex_list$colour <-  if_else(residue_vertex_list$colour == "J", "Aliphatic", 
                                       if_else(residue_vertex_list$colour == "O", "Aromatic",
                                               if_else(residue_vertex_list$colour == "U", "Acidic", 
                                                       if_else(residue_vertex_list$colour == "B", "Basic",
                                                               if_else(residue_vertex_list$colour == "X", "Hydroxilic",
                                                                       if_else(residue_vertex_list$colour == "Y", "Sulfuric",
                                                                               "Amidic"))))))

residue_vertex_list$colour <- factor(residue_vertex_list$colour)

residue_vertex_list$counts[is.nan(residue_vertex_list$counts)] <- 5.5

residue_edge_list <- rep(list(NULL), times = (max(exp_group[[8]]$CDR3_length) - 1))

for (l in 1:(max(exp_group[[8]]$CDR3_length) - 1)) {
  residue_edge_list[[l]] <- table(residue_profile[, l + match("A", names(residue_profile)) - 1], residue_profile[,l + match("A", names(residue_profile))]) %>%
    as.data.frame(stringsAsFactors = F) %>%
    melt() %>%
    filter(value != 0) %>%
    select(-variable)
}

residue_edge_list <- residue_edge_list %>%
  bind_rows()

residue_edge_list <- residue_edge_list[str_detect(residue_edge_list$Var1, "[A-Z]+"),]

residue_edge_list <- residue_edge_list[str_detect(residue_edge_list$Var2, "[A-Z]+"),]

residue_edge_list <- residue_edge_list %>%
  separate(Var1, into = c("colour", "position"), sep = 1, remove = F) %>%
  group_by(position) %>%
  mutate(weight = (((value - min(value)) / (max(value) - min(value))) + 0.1) * 5) %>%
  ungroup() %>%
  select(Var1, Var2, weight, position, value)

residue_edge_list$weight[is.nan(residue_edge_list$weight)] <- 5.5

residue_edge_list$edge_colour <- if_else(residue_edge_list$weight == 5.5, residue_edge_list$edge_colour <- 2, residue_edge_list$edge_colour <- 1)

residue_edge_list$edge_size <- if_else(residue_edge_list$weight == 5.5, residue_edge_list$edge_size <- 2, residue_edge_list$edge_size <- 0.25)

residue_edge_list <- residue_edge_list %>%
  filter(weight >= 1)

joined_nodes <- c(residue_edge_list$Var1, residue_edge_list$Var2) %>%
  as.data.frame(stringsAsFactors = F) %>%
  distinct()

residue_vertex_list <- residue_vertex_list %>%
  right_join(joined_nodes, by = c("residue" = "."))

graph <- residue_edge_list %>%
  graph_from_data_frame(vertices = residue_vertex_list)

ggnet2(graph, node.size = V(graph)$counts,
       arrow.size = 5,
       node.color = V(graph)$colour,
       color.palette = "Set3",
       edge.size = E(graph)$edge_size,
       edge.color = E(graph)$edge_colour,
       arrow.gap = 0.02) +
  guides(size = "none")

#####
# overlapping PIDs and clones
#####

neg_PIDs <- read.csv("D://mice_experiments//neg_control//output//neg_PIDs.csv") %>%
  select(-X)

setwd("D://RNAseq//1239Shp15b_Joost_final//PID summary")

exp_PIDs <- map(dir()[str_detect(dir(), " miXCR clones read names with PIDs.csv")], read.csv, stringsAsFactors = F)

exp_PIDs <- exp_PIDs %>%
  map(select, -`Well`, -`Clone.Id`, -`read.name`) %>%
  map(distinct)

exp_PIDs <- exp_PIDs %>%
  bind_rows()

exp_PIDs <- exp_PIDs %>%
  filter(!str_detect(CDR3, "_") | !str_detect(CDR3, "\\*") | !str_length(CDR3) > 20 | !str_length(CDR3) < 8)

exp_PIDs <- exp_PIDs %>%
  group_by(Subject.Id, CDR3) %>%
  mutate(PID.count = n_distinct(PID)) %>%
  filter(PID.count > 3)

exp_PIDs <- exp_PIDs %>%
  ungroup() %>%
  left_join(neg_PIDs, by = c("CDR3", "PID"))

PID_overlap <- exp_PIDs %>%
  filter(!is.na(exp)) %>%
  distinct()

PID_overlap_summary <- PID_overlap %>%
  group_by(Subject.Id, CDR3) %>%
  summarise(matching_PIDs = n_distinct(PID))

#####
