# PID matching
#####

setwd("D:/data/experiments/neg_control/1239Shp16_JK_final/PID summary") # set the path the directory that contains your TCRseq PID files

names_exp <- dir()[str_detect(dir(), " miXCR clones read names with PIDs.csv")] %>%
  str_remove(" miXCR clones read names with PIDs.csv")

raw_exp <- dir()[str_detect(dir(), " miXCR clones read names with PIDs.csv")] %>%
  map(read.csv, stringsAsFactors = F)

#setwd("/Users/joelkidman/Documents/PhD_documents/RNAseq/")

setwd("D://data//experiments//neg_control//output") # change directory to where you want the pictures to go

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
# relationship with other samples
#####

neg_ctl <- PID_control("D://data//experiments//neg_control//1239Shp16_JK_final")

neg <- neg_4_1_0 %>%
  bind_rows(neg_4_2_0) %>%
  bind_rows(neg_KR) %>%
  bind_rows(neg_RENCA) %>%
  bind_rows(neg_ctl)

neg_summary <- neg %>%
  group_by(machine_id, run_number) %>%
  summarise(n_PIDs = sum(matching_PIDs))

#####