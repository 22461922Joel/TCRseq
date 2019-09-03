library(readxl)
library(grid)
library(gridExtra)
library(grDevices)
library(RColorBrewer)
library(tidyverse)

# tumour growth function expects a directory containing an .xls 
# spreadsheet from the export list version history tab from sharepoint.
# the function will try and guess the right columns to include in the tumoumr growth.
# odd or duplicate column names will throw it off
# it will create a .csv file with clean data from each mouse for further analysis
# as well as a graph to check correctness

setwd("D:/data/experiments/4_1_0/tumour_growth")

tg <- dir()[1] %>% 
  read_xls() %>%
  select(ends_with("L", ignore.case = F),
         ends_with("W", ignore.case = F),
         `Item ID`,
         `Version No.`,
         starts_with("Exp"),
         "day" = starts_with("day", ignore.case = F),
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
  filter(day < 100)

dual_tumour_growth <- function(directory) {
  
  setwd(directory) 
  
  tg <- dir()[1] %>% 
    read_xls() %>%
    select(ends_with("L", ignore.case = F),
           ends_with("W", ignore.case = F),
           `Item ID`,
           `Version No.`,
           starts_with("Exp"),
           "day" = starts_with("day", ignore.case = F),
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
    filter(day < 100)
  
  setwd(paste(getwd(), "/.."))
  
  write.csv(tg, file = "tumour_growth.csv", row.names = F)
  l_text <- "LEFT"
  r_text <- "RIGHT"
  l_grob <- grid.text(l_text, just = c(1,1), x = 1, y = 1, gp = gpar(fontsize = 14))
  r_grob <- grid.text(r_text, just = c(1,0), x = 1, y = 0, gp = gpar(fontsize = 14))
  Line <- geom_line(size = 1, lineend = "round")
  
  colour_scale <- scale_colour_manual(values = brewer.pal(8, "Dark2"))
  
  D <- ggplot(tg, aes(day, diff_size, group = mouse, colour = mouse)) +
    Line + 
    theme(legend.position = "none") +
    labs(y = "proportional tumour size difference",
         x = "days post innoculation") +
    annotation_custom(l_grob) +
    annotation_custom(r_grob) +
    colour_scale +
    ylim(c(-1,1))
  
  y_limits <- c(0, max(tg$tumour_size))
  
  
  L <- ggplot(filter(tg, flank == "L"), aes(day, tumour_size, group = mouse, colour = mouse)) +
    Line +
    scale_x_reverse() +
    theme(legend.position = "none") +
    ylim(y_limits) +
    colour_scale +
    labs(y = expression(paste("tumour size ", mm^2)),
         x = "days post innoculation")
  
  R <- ggplot(filter(tg, flank == "R"), aes(day, tumour_size, group = mouse, colour = mouse)) +
    Line +
    scale_y_continuous(position = "right", limits = y_limits) +
    theme(legend.position = c(0.05, 0.95),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1)) +
    colour_scale +
    labs(y = expression(paste("tumour size ", mm^2)),
         x = "days post innoculation")
  
  lay <- rbind(c(1,2))
  
  grid.newpage()
  
  grid.arrange(arrangeGrob(L, R, layout_matrix = lay))
}

single_tumour_growth <- function(directory, columns) {
  setwd(directory) 
  
  tg <- dir()[1] %>% 
    read_xls(col_types = columns) %>%
    fill(`Item ID`) %>%
    filter(!is.na(`Version No.`)) %>%
    group_by(`Item ID`) %>%
    fill(`Animal ID`, Box, `Expt ID`, .direction = "up") %>%
    filter(!is.na(`Expt ID`)) %>%
    unite("mouse", Box, `Animal ID`, sep = ".") %>%
    group_by(mouse) %>%
    fill(`R Tumour L`, `R Tumour W`, `L Tumour L`, `L Tumour W`, `days (tumour)`, .direction = "up") %>%
    filter(!is.na(`R Tumour L`)) %>%
    group_by(`days (tumour)`, mouse) %>%
    mutate(`R Tumour L` = `R Tumour L`[match(max(`Version No.`), `Version No.`)],
           `R Tumour W` = `R Tumour W`[match(max(`Version No.`), `Version No.`)],
           `L Tumour L` = `L Tumour L`[match(max(`Version No.`), `Version No.`)],
           `L Tumour W` = `L Tumour W`[match(max(`Version No.`), `Version No.`)]) %>%
    select(-`Version No.`) %>%
    distinct() %>%
    mutate(L_Tumour_size = `L Tumour W` * `L Tumour L`, 
           R_Tumour_size = `R Tumour W` * `R Tumour L`) %>%
    drop_na() %>%
    mutate(diff_size = ((L_Tumour_size - R_Tumour_size) / (L_Tumour_size + R_Tumour_size))) %>%
    gather(key = "flank", value = "tumour_size", L_Tumour_size, R_Tumour_size) %>%
    mutate(log_tumour_size = log(tumour_size)) %>%
    unite("mouse_flank", mouse, flank, remove = F)
  
  
  tg$`days (tumour)` <- str_remove(tg$`days (tumour)`, "float;#") %>% as.numeric()
  tg$diff_size[is.na(tg$diff_size)] <- 0
}


setwd("D://mice_experiments//4_1_0//tumour_growth") #`R Tumour L`, `R Tumour W`, `L Tumour L`, `L Tumour W`, 

tg <- dir()[1] %>% 
  read_xls(col_types = columns) %>%
  fill(`Item ID`) %>%
  filter(!is.na(`Version No.`)) %>%
  group_by(`Item ID`) %>%
  fill(`Animal ID`, Box, `Expt ID`, .direction = "up") %>%
  filter(!is.na(`Expt ID`)) %>%
  unite("mouse", Box, `Animal ID`, sep = ".") %>%
  group_by(mouse) %>%
  fill(`R Tumour L`, `R Tumour W`, `L Tumour L`, `L Tumour W`, `days (tumour)`, .direction = "up") %>%
  filter(!is.na(`R Tumour L`)) %>%
  group_by(`days (tumour)`, mouse) %>%
  mutate(`R Tumour L` = `R Tumour L`[match(max(`Version No.`), `Version No.`)],
         `R Tumour W` = `R Tumour W`[match(max(`Version No.`), `Version No.`)],
         `L Tumour L` = `L Tumour L`[match(max(`Version No.`), `Version No.`)],
         `L Tumour W` = `L Tumour W`[match(max(`Version No.`), `Version No.`)]) %>%
  select(-`Version No.`) %>%
  distinct() %>%
  mutate(L_Tumour_size = `L Tumour W` * `L Tumour L`, 
         R_Tumour_size = `R Tumour W` * `R Tumour L`) %>%
  drop_na() %>%
  mutate(diff_size = ((L_Tumour_size - R_Tumour_size) / (L_Tumour_size + R_Tumour_size))) %>%
  gather(key = "flank", value = "tumour_size", L_Tumour_size, R_Tumour_size) %>%
  mutate(log_tumour_size = log(tumour_size)) %>%
  unite("mouse_flank", mouse, flank, remove = F)

tg$`days (tumour)` <- str_remove(tg$`days (tumour)`, "float;#") %>% as.numeric()
tg$diff_size[is.na(tg$diff_size)] <- 0

setwd("D://mice_experiments//4_1_0")

write.csv(tg, file = "tumour_growth.csv")

ggplot(tg, aes(`days (tumour)`, tumour_size, group = mouse_flank)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 160), 
                     breaks = seq(from = 0, to = 160, by = 20), 
                     minor_breaks = seq(from = 10, to = 150, by = 20)) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(from = 0, to = 70, by = 5)) +
  labs(x = "days post inoculation",
       y = expression(paste("tumour size ", mm^2)))

l_text <- "LEFT"
r_text <- "RIGHT"
l_grob <- grid.text(l_text, just = c(1,1), x = 1, y = 1, gp = gpar(fontsize = 14))
r_grob <- grid.text(r_text, just = c(1,0), x = 1, y = 0, gp = gpar(fontsize = 14))
Line <- geom_line(size = 1, lineend = "round")

colour_scale <- scale_colour_manual(values = brewer.pal(8, "Dark2"))

D <- ggplot(tg, aes(`days (tumour)`, diff_size, group = mouse, colour = mouse)) +
  Line + 
  theme(legend.position = "none") +
  labs(y = "proportional tumour size difference",
       x = "days post innoculation") +
  annotation_custom(l_grob) +
  annotation_custom(r_grob) +
  colour_scale +
  ylim(c(-1,1))

y_limits <- c(0, max(tg$tumour_size))


L <- ggplot(filter(tg, flank == "L_Tumour_size"), aes(`days (tumour)`, tumour_size, group = mouse, colour = mouse)) +
  Line +
  scale_x_reverse() +
  theme(legend.position = "none") +
  ylim(y_limits) +
  colour_scale +
  labs(y = expression(paste("tumour size ", mm^2)),
       x = "days post innoculation")

R <- ggplot(filter(tg, flank == "R_Tumour_size"), aes(`days (tumour)`, tumour_size, group = mouse, colour = mouse)) +
  Line +
  scale_y_continuous(position = "right", limits = y_limits) +
  theme(legend.position = c(0.05, 0.95),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1)) +
  colour_scale +
  labs(y = expression(paste("tumour size ", mm^2)),
       x = "days post innoculation")

lay <- rbind(c(1,2))

tumour_growth <- grid.arrange(arrangeGrob(L, R, layout_matrix = lay))

tg_log <- tg %>%
  filter(tumour_size > 0)

y_limits <- c(0, max(tg_log$log_tumour_size))


L <- ggplot(filter(tg_log, flank == "L_Tumour_size"), aes(`days (tumour)`, log_tumour_size, group = mouse, colour = mouse)) +
  Line +
  scale_x_reverse() +
  theme(legend.position = "none") +
  ylim(y_limits) +
  colour_scale +
  labs(y = expression(paste("tumour size ", mm^2)),
       x = "days post innoculation")

R <- ggplot(filter(tg_log, flank == "R_Tumour_size"), aes(`days (tumour)`, log_tumour_size, group = mouse, colour = mouse)) +
  Line +
  scale_y_continuous(position = "right", limits = y_limits) +
  theme(legend.position = c(0.05, 0.95),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1)) +
  colour_scale +
  labs(y = expression(paste("tumour size ", mm^2)),
       x = "days post innoculation")

lay <- rbind(c(1,2))

tumour_growth <- grid.arrange(arrangeGrob(L, R, layout_matrix = lay))

ggsave("tumour_growth.png", plot = tumour_growth)

remove(columns, l_text, r_text, l_grob, r_grob, Line, D, y_limits, L, R, lay, tumour_growth)

tg_14 <- tg %>%
  select(mouse, flank, tumour_size, "days" = `days (tumour)`) %>%
  as_tibble() %>%
  filter(days == 14) %>%
  separate(flank, into = c("flank", NA, NA), sep = "_") %>%
  as_tibble()
