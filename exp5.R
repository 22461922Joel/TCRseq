
AE140_columns <- c("numeric", 
                   "numeric",
                   rep("skip", times = 2),
                   "text",
                   "text",
                   "numeric",
                   rep("skip", times = 7),
                   "numeric",
                   "numeric",
                   rep("skip", times = 6),
                   "numeric",
                   rep("skip", times = 5),
                   "text",
                   rep("skip", times = 8),
                   "text")

single_tumour_growth <- function(directory, columns) {
  setwd(directory) 
  
  tg <- dir()[1] %>% 
    read_xls(col_types = columns) %>%
    fill(`Item ID`) %>%
    filter(!is.na(`Version No.`)) %>%
    group_by(`Item ID`) %>%
    fill(`Animal ID`, Box, `Expt ID`, `Treatment Type`, .direction = "up") %>%
    fill(`Treatment Type`) %>%
    filter(!is.na(`Expt ID`)) %>%
    unite("mouse", Box, `Animal ID`, sep = ".") %>%
    group_by(mouse) %>%
    fill(`R Tumour L`, `R Tumour W`, `Days Post Inoc`, .direction = "up") %>%
    filter(!is.na(`R Tumour L`)) %>%
    group_by(`Days Post Inoc`, mouse) %>%
    mutate(`R Tumour L` = `R Tumour L`[match(max(`Version No.`), `Version No.`)],
           `R Tumour W` = `R Tumour W`[match(max(`Version No.`), `Version No.`)]) %>%
    select(-`Version No.`) %>%
    distinct() %>%
    mutate(R_Tumour_size = `R Tumour W` * `R Tumour L`)
  
  
  tg$`Days Post Inoc` <- str_remove(tg$`Days Post Inoc`, "float;#") %>% as.numeric()
  
  tg$`Treatment Type` <- if_else(tg$`Treatment Type` == "100ug aPD-L1, 100 ug aCTLA-4 in 200uL IP", "CP", "PBS")
  
  tg
}

tg <- single_tumour_growth("D://data//experiments//5_0_0//tumour_growth", AE140_columns)


setwd("D://data//experiments//5_0_0//tumour_growth") 

tg <- dir()[1] %>% 
  read_xls(col_types = AE140_columns) %>%
  fill(`Item ID`) %>%
  filter(!is.na(`Version No.`)) %>%
  group_by(`Item ID`) %>%
  fill(`Animal ID`, Box, `Expt ID`, Strain, `Treatment Type`, .direction = "up") %>%
  filter(!is.na(`Expt ID`)) %>%
  unite("mouse", Box, `Animal ID`, sep = ".") %>%
  group_by(mouse) %>%
  mutate(`Days Post Inoc` = str_remove(`Days Post Inoc`, "float;#"),
         `Days Post Inoc` = as.numeric(`Days Post Inoc`)) %>%
  fill(`R Tumour L`, `R Tumour W`, `Days Post Inoc`, .direction = "up") %>%
  fill(`Treatment Type`) %>%
  filter(!is.na(`R Tumour L`)) %>%
  group_by(`Days Post Inoc`, mouse, Strain) %>%
  mutate(`R Tumour L` = `R Tumour L`[match(max(`Version No.`), `Version No.`)],
         `R Tumour W` = `R Tumour W`[match(max(`Version No.`), `Version No.`)]) %>%
  select(-`Version No.`) %>%
  distinct() %>%
  mutate(R_Tumour_size = `R Tumour W` * `R Tumour L`) %>%
  group_by(`Days Post Inoc`, Strain)



tg$`Treatment Type` <- if_else(tg$`Treatment Type` == "100ug aPD-L1, 100 ug aCTLA-4 in 200uL IP", "CP", "PBS")

ggplot(tg, aes(`Days Post Inoc`, R_Tumour_size, group = `Treatment Type`, colour = `Treatment Type`)) +
  # geom_line() +
  # geom_point() +
  geom_smooth() +
  facet_grid(. ~ Strain)

ggplot(tg, aes(`Days Post Inoc`, n_alive, group = `Treatment Type`, colour = `Treatment Type`)) +
  geom_step() +
  facet_grid(`Treatment Type` ~ Strain)

tg %>%
  filter(`Days Post Inoc` == 10) %>%
  ggplot(aes(`Treatment Type`, R_Tumour_size, fill = `Treatment Type`)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  facet_grid(. ~ Strain) +
  theme(axis.title.x = element_blank()) +
  labs(y = "tumour size at first treatment")

write.csv(tg, "RENCA_CPB_tumour_growth.csv")

