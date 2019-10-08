
working_path <- file.path(getwd(), "experiments/4_0_0")

#####
# tumour growth for grow out and monitor mice
#####

tg <- file.path(working_path, "AE140_joel.xlsx") %>% 
  read_excel() %>%
  select(ends_with("L", ignore.case = F),
         ends_with("W", ignore.case = F),
         `Item ID`,
         `Version No.`,
         starts_with("Exp"),
         "day" = starts_with("day", ignore.case = T),
         Box,
         starts_with("Animal")) %>%
  fill(`Item ID`) %>%
  filter(!is.na(`Version No.`)) %>%
  group_by(`Item ID`) %>%
  fill(`Animal ID`, Box, `Expt ID`, .direction = "up") %>%
  filter(!is.na(`Expt ID`)) %>%
  unite("mouse", Box, `Animal ID`, sep = "-") %>%
  group_by(mouse) %>%
  fill(`R Tumour L`, `R Tumour W`, `L Tumour L`, `L Tumour W`, day, .direction = "up") %>%
  fill(`R Tumour L`, `R Tumour W`, `L Tumour L`, `L Tumour W`, day) %>%
  filter(!is.na(`R Tumour L`)) %>%
  group_by(day, mouse) %>%
  mutate(`R Tumour L` = `R Tumour L`[match(max(`Version No.`), `Version No.`)],
         `R Tumour W` = `R Tumour W`[match(max(`Version No.`), `Version No.`)],
         `L Tumour L` = `L Tumour L`[match(max(`Version No.`), `Version No.`)],
         `L Tumour W` = `L Tumour W`[match(max(`Version No.`), `Version No.`)]) %>%
  select(-`Version No.`) %>%
  distinct() %>%
  mutate(L = as.numeric(`L Tumour W`) * as.numeric(`L Tumour L`), 
         R = as.numeric(`R Tumour W`) * as.numeric(`R Tumour L`)) %>%
  drop_na() %>%
  gather(key = "flank", value = "tumour_size", L, R) %>%
  ungroup() %>%
  mutate(day = str_remove(.$day, "float;#") %>% as.numeric(),
         mouse = as.character(mouse)) %>%
  filter(day < 100) %>%
  separate(mouse, into = c("box", "mouse"), sep = "-", remove = F)

tg_aes <- aes(day, tumour_size, group = interaction(mouse, box), colour = box)

tg_line <- geom_line(size = 1)

tg_treatment <- geom_vline(xintercept = c(11, 14))

R <- ggplot(tg %>% filter(flank == "R"), tg_aes) +
  tg_line +
  theme_bw() +
  theme(legend.position = c(0.1, 0.55)) +
  scale_y_continuous(position = "right", limits = c(0, max(tg$tumour_size))) +
  tg_treatment

L <- ggplot(tg %>% filter(flank == "L"), tg_aes) +
  tg_line +
  scale_x_reverse() +
  theme_bw() +
  theme(legend.position = "none") +
  tg_treatment +
  ylim(0, max(tg$tumour_size))

grid.arrange(L, R, ncol = 2)  
