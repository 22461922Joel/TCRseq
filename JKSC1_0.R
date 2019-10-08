library(readxl)
library(lubridate)

working_path <- file.path(getwd(), "experiments/JKSC1.0")

tg <- read_excel(file.path(working_path, "AE140.xlsx")) %>%
  select(ends_with("L", ignore.case = F),
         ends_with("W", ignore.case = F),
         `Item ID`,
         `Version No.`,
         starts_with("Exp"),
         "day" = contains("monitor", ignore.case = T),
         Box,
         starts_with("Animal"),
         contains("treatment type")) %>%
  mutate(treatment_group = if_else(str_detect(`Treatment Type`, "PBS"), "PBS", "combo")) %>%
  select(-ends_with("type")) %>%
  fill(`Item ID`) %>%
  filter(!is.na(`Version No.`)) %>%
  group_by(`Item ID`) %>%
  fill(`Animal ID`, Box, `Expt ID`, .direction = "up") %>%
  filter(!is.na(`Expt ID`)) %>%
  unite("mouse", Box, `Animal ID`, sep = "-") %>%
  group_by(mouse) %>%
  fill(`R Tumour L`, `R Tumour W`, `L Tumour L`, `L Tumour W`, day, treatment_group, .direction = "up") %>%
  fill(`R Tumour L`, `R Tumour W`, `L Tumour L`, `L Tumour W`, day, treatment_group) %>%
  filter(!is.na(`R Tumour L`)) %>%
  group_by(day, mouse) %>%
  mutate(`R Tumour L` = `R Tumour L`[match(max(`Version No.`), `Version No.`)],
         `R Tumour W` = `R Tumour W`[match(max(`Version No.`), `Version No.`)],
         `L Tumour L` = `L Tumour L`[match(max(`Version No.`), `Version No.`)],
         `L Tumour W` = `L Tumour W`[match(max(`Version No.`), `Version No.`)],
         treatment_group = treatment_group[match(max(`Version No.`), `Version No.`)]) %>%
  select(-`Version No.`) %>%
  distinct() %>%
  mutate(L = as.numeric(`L Tumour W`) * as.numeric(`L Tumour L`), 
         R = as.numeric(`R Tumour W`) * as.numeric(`R Tumour L`)) %>%
  drop_na() %>%
  gather(key = "flank", value = "tumour_size", L, R) %>%
  ungroup() %>%
  mutate(day = str_remove(.$day, "datetime;#") %>% ymd_hms() %>% round_date(unit = "day"),
         mouse = as.character(mouse)) %>%
  separate(mouse, into = c("box", "mouse"), sep = "-") %>%
  mutate(treatment_group = if_else(box == "1", "PBS", "combo")) %>%
  unite("mouse", box, mouse, sep = "-", remove = F)

tg$day <- as.interval(tg$day, start = read_excel(file.path(working_path, "AE140.xlsx")) %>%
                        select(contains("tumour inoc")) %>%
                        na.omit() %>%
                        unique() %>%
                        as.character() %>%
                        dmy()) %>%
  as.period(unit = "day") %>%
  as.numeric("days") * -1

tg <- tg %>%
  filter(box != "2" | flank != "L" | day <= 13, day != 19) %>%
  filter(mouse != "3-1" | flank == "R" | day <= 13) %>%
  filter(mouse != "3-2" | flank == "R" | day <= 13) %>%
  filter(mouse != "3-3" | flank == "R" | day <= 13)

write_csv(tg, file.path(working_path, "tumour_growth.csv"))

tg <- read_csv(file.path(working_path, "tumour_growth.csv")) %>%
  filter(day <= 16 | day > 17)

tg_line <- geom_line(size = 2)

tg_aes <- aes(day, tumour_size, group = interaction(mouse, flank), colour = interaction(treatment_group, mouse))

y_axis <- ylim(c(0, max(tg$tumour_size)))

treatment <- geom_vline(xintercept = c(7, 9, 11))

surgery <- geom_vline(xintercept = 13, colour = "red")

R <- ggplot(tg %>% filter(flank == "R"), tg_aes) +
  tg_line +
  theme_bw() +
  theme(legend.position = c(0.1, 0.55), legend.title = element_blank()) +
  scale_y_continuous(position = "right", limits = c(0, max(tg$tumour_size))) +
  treatment

L <- ggplot(tg %>% filter(flank == "L"), tg_aes) +
  tg_line +
  scale_x_reverse() +
  theme_bw() +
  theme(legend.position = "none") +
  y_axis +
  treatment +
  surgery

grid.arrange(L, R, ncol = 2)  

