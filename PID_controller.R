library(tidyverse)

#####
# data cleaning
#####

setwd("D://mice_experiments//neg_control//1239Shp16_JK_final//PID summary") # set the path the directory that contains your TCRseq PID files

names_exp <- dir()[str_detect(dir(), " miXCR clones with PID read cutoff 1.tsv")] %>%
  str_remove(" miXCR clones with PID read cutoff 1.tsv")
  
raw_exp <- dir()[str_detect(dir(), " miXCR clones with PID read cutoff 1.tsv")] %>%
  map(read.delim, stringsAsFactors = F)

#setwd("/Users/joelkidman/Documents/PhD_documents/RNAseq/")

setwd("D://mice_experiments//neg_control//output//") # change directory to where you want the pictures to go

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

remove(raw_exp)

#remove stop and out of frame sequences, rejects to separate file

reject_vector <- str_detect(exp$aaSeqCDR3, "_") | 
  str_detect(exp$aaSeqCDR3, "\\*") | 
  str_length(exp$aaSeqCDR3) > 20 | 
  str_length(exp$aaSeqCDR3) < 8 |
  exp$cloneCount < 3 |
  exp$PID.count < 3 #|
#  nSeqCDR3 != "TGTGCCAGCGGTGAGACAGGGACCAACGAAAGATTATTTTTC"

exp_rejects <- exp %>%
  filter(reject_vector)

write.csv(exp_rejects, "rejected_CDR3s.csv")

remove(exp_rejects)

exp_clean <- exp %>%
  filter(!reject_vector)

remove(exp, reject_vector)


#####
# PID matching
#####

setwd("D://mice_experiments//neg_control//1239Shp16_JK_final//PID summary") # set the path the directory that contains your TCRseq PID files

names_exp <- dir()[str_detect(dir(), " miXCR clones read names with PIDs.csv")] %>%
  str_remove(" miXCR clones read names with PIDs.csv")

raw_exp <- dir()[str_detect(dir(), " miXCR clones read names with PIDs.csv")] %>%
  map(read.csv, stringsAsFactors = F)

#setwd("/Users/joelkidman/Documents/PhD_documents/RNAseq/")

setwd("D://mice_experiments//neg_control//output") # change directory to where you want the pictures to go

names(raw_exp) <- names_exp

neg_PIDs <- raw_exp %>%
  as.data.frame()

names(neg_PIDs) <- c("exp", "Well", "clone.Id", "clone.count", "CDR3", "PID", "read.name")

neg_PIDs <- neg_PIDs %>%
  filter(clone.count > 3) %>%
  select(-read.name, -Well, -clone.Id, -clone.count) %>%
  distinct() %>%
  group_by(CDR3) %>%
  mutate(PID.count = n_distinct(PID)) %>%
  ungroup()

reject_vector <- str_detect(neg_PIDs$CDR3, "_") | 
  str_detect(neg_PIDs$CDR3, "\\*") | 
  str_length(neg_PIDs$CDR3) > 20 | 
  str_length(neg_PIDs$CDR3) < 8  |
  neg_PIDs$PID.count < 3

exp_rejects <- neg_PIDs %>%
  filter(reject_vector)

neg_PIDs <- neg_PIDs %>%
  filter(!reject_vector)

remove(reject_vector, exp_rejects)

write.csv(neg_PIDs, "neg_PIDs.csv")
  
#####