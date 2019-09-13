library(tidyverse)
library(seqinr)
library(furrr)

path <- file.path(getwd(), "experiments/RNAseq/1239Shp15b_Joost_final/PID_clusters")



# files <- paste0(path, "/", dir(path))[str_detect(paste0(path, "/", dir(path)), ".*cluster[1-9]?[0-9]?[0-9]?_.*")] # for subsets

files <- paste0(path, "/", dir(path))[!str_detect(paste0(path, "/", dir(path)), "nohup")]

plan(multiprocess)

clusters <- future_map(files, read.fasta, seqtype = "AA", as.string = T) %>%
  future_map(getSequence, as.string = T)

cluster_names <- files %>%
  str_remove(".*_clusters/PID_") %>%
  str_remove("_t.*")

names(clusters) <- cluster_names

pid_clusters <- future_map(clusters, unlist) %>%
  future_map(as.data.frame, stringsAsFactors = F) %>%
  bind_rows(.id = "cluster") %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric())

names(pid_clusters) <- c("PID_cluster", "aaSeqCDR3")

rm(clusters, cluster_names, files, path)
