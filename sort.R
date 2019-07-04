library(lubridate)
library(tidyverse)
library(tabulizer)
library(plotly)
library(reshape2)

setwd("D://mice_Experiments//4_1_0//4.1_sort_data")

sort_area <- list(c(260, 40, 288, 265), c(7, 13, 26, 402))
sort_data <- dir() %>% map(extract_tables, area = sort_area, guess = F, pages = c(1,1))
test <- do.call("rbind", sort_data[1:length(sort_data)])
sort_data_1 <- do.call("rbind", test[1:(length(test)/2)])
sort_data_2 <- do.call("rbind", test[(length(test)/2 + 1):length(test)])
sort_data_1 <- cbind(sort_data_1, rep(sort_data_2[,2], each = 2))
sort_data_1 <- cbind(sort_data_1, rep(sort_data_2[,1], each = 2))
sort_data_final <- as.data.frame(sort_data_1, stringsAsFactors = F)
remove(sort_data_1, sort_data_2, sort_data, test, sort_area)


percent <- function(x) {x / 100}

sort_colnames <- c("population", "target", "cells", "rate", "efficiency", "sort_time", "exp_mouse_tissueFlank", "date")

names(sort_data_final) <- sort_colnames

sort_data_final$date <- mdy(sort_data_final$date)

sort_data_final[, 2:4] <- map(sort_data_final[, 2:4], str_remove, pattern = ",")

sort_data_final[, 2:4] <- map(sort_data_final[, 2:4], as.numeric)

sort_data_final$sort_time <- sort_data_final$sort_time %>% 
  ms() %>%
  as.duration() %>%
  as.numeric()

sort_data_final$efficiency <- sort_data_final$efficiency %>% str_remove("%") %>% as.numeric() %>% percent()

sort_data_final$exp_mouse_tissueFlank <- str_replace_all(sort_data_final$exp_mouse_tissueFlank, "Â", "_")

sort_data_final <- sort_data_final %>% 
  separate(exp_mouse_tissueFlank, c("exp", "mouse", "tissue"), "_", remove = F) %>% 
  separate(tissue, c("tissue", "flank"), 2) %>%
  as_tibble()

sort_data_final$mouse <- str_trunc(sort_data_final$mouse, width = 3, side = "left", ellipsis = "")

sort_data_final$tissue <- str_trunc(sort_data_final$tissue, width = 1, side = "left", ellipsis = "")

str_length(sort_data_final$tissue)

sort_data_final <- sort_data_final %>%
  filter(mouse != "1.1")

#####
# IIID lot specific
#####

IIID_subject_ID <- sort_data_final %>%
  filter(str_detect(exp_mouse_tissueFlank, "2."),
         str_detect(tissue, "T")) %>%
  separate(exp_mouse_tissueFlank, c("exp_mouse"), 10, remove = F) %>%
  unite("exp_mouse_tissueFlank_pop", exp_mouse_tissueFlank, population) %>%
  select(exp_mouse_tissueFlank_pop, exp_mouse) %>%
  arrange(exp_mouse_tissueFlank_pop)

setwd("D://mice_Experiments//4_2_0")

write.csv(IIID_subject_ID, "IIID_subject_IDs.csv")

#####
# Sort spreadsheet info
#####

sort_tumours <- sort_data_final %>%
  filter(str_detect(tissue, "T"),
         mouse != "1.1") %>%
  left_join(tg_14) %>%
  mutate(norm_tumour_size = (tumour_size - min(tumour_size)) / (max(tumour_size) - min(tumour_size)),
         norm_cells = ((cells - min(cells)) / (max(cells) - min(cells))))

#####, by = c("mouse" = "mouse", "flank" = "flank")
# PCR page
#####

PCR_export <- sort_tumours %>%
  filter(mouse != "1.1") %>%
  select(-target, -rate, -efficiency, -sort_time, -exp:-date, -days:-norm_cells)

primers <- rep(1:12, length.out = length(PCR_export$exp_mouse_tissueFlank))

tube <- rep(1:8, length.out = length(PCR_export$exp_mouse_tissueFlank))

PCR_export$primers <- primers

PCR_export$tube <- tube

write.csv(PCR_export, "PCR_primers.csv")

#####
# storage locale
#####

setwd("D://mice_Experiments")

storage_perm <- read.csv("storage_locale.csv", as.is = T) %>%
  select(-X) %>%
  mutate(mouse = as.character(mouse),
         date = ymd(date))

storage <- sort_tumours %>%
  select(population, exp:date)

storage$tissue[storage$tissue == "T"] <- "tum"

last_locale <- max(storage_perm$PCR1_location)

new_locale <- vector()

new_locale[1] <- 1

for (i in 1:(2 *length(storage$tissue))) {
  new_locale[i] <- if_else(max(new_locale) < 96, 
                           (last_locale + i), 
                           as.integer(last_locale + i - 96))
}

RNA_box <- vector()

last_box <- max(storage_perm$PCR1_box)

for (i in 1:length(storage$tissue)) {
  RNA_box[i] <- if_else(new_locale[i] == 1,
                         as.integer(last_box + 1),
                         last_box)
  last_box <- if_else(new_locale[i] == 1,
                      as.integer(last_box + 1),
                      last_box)
}

PCR1_box <- vector()

for (i in 1:length(storage$tissue)) {
  PCR1_box[i] <- if_else(new_locale[i + length(storage$tissue)] == 1,
                         as.integer(last_box + 1),
                         last_box)
  last_box <- if_else(new_locale[i + length(storage$tissue)] == 1,
                      as.integer(last_box + 1),
                      last_box)
}

storage$RNA_box <- RNA_box

storage$RNA_location <- new_locale[1:length(storage$population)]

storage$PCR1_box <- PCR1_box

storage$PCR1_location <- new_locale[(1 + length(storage$population)):(2 * length(storage$population))]

storage_perm <- storage_perm %>%
  bind_rows(storage)

write.csv(storage_perm, "storage_locale.csv")

setwd("D://mice_Experiments//4_1_0")

write.csv(storage, "Storage.csv")

#####


ggplot(sort_tumours, aes(norm_tumour_size, norm_cells)) +
  geom_point(size = 3)

ggplot(sort_tumours, aes(population, cells)) +
  geom_boxplot() +
  labs(y = "sorted cells", x = "T cell subset")

ggplot(sort_tumours, aes(flank, evenness, group = mouse_pop)) +
  geom_point() +
  scale_y_log10() +
  geom_line(size = 1) +
  labs(y = "total returned clones")

resolution = 0.1

umi_mod <- lm(norm_evenness ~ norm_cells * norm_tumour_size, data = sort_tumours)

axis_x <- seq(min(sort_tumours$norm_cells), max(sort_tumours$norm_cells), by = resolution)
axis_y <- seq(min(sort_tumours$norm_tumour_size), max(sort_tumours$norm_tumour_size), by = resolution)
surface <- expand.grid(norm_cells = axis_x, norm_tumour_size = axis_y, KEEP.OUT.ATTRS = F)
surface$norm_evenness <- predict.lm(umi_mod, newdata = surface)
surface <- acast(surface, norm_cells ~ norm_tumour_size, value.var = "norm_evenness")

plot_ly(data = sort_tumours, x = ~norm_cells, y = ~norm_tumour_size, z = ~norm_evenness) %>%
  add_markers() %>%
  add_surface(z = surface, x = axis_x, y = axis_y)

sort_tumour_umi_limits <- c(1000, 50000)

sort_tumour_umi <- ggplot(sort_tumours, aes(cells, evenness, size = tumour_size, colour = mouse, shape = flank)) +
  geom_point() +
  scale_y_log10(limits = sort_tumour_umi_limits) +
  scale_x_log10(limits = sort_tumour_umi_limits) +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "total clones")

ggplotly(sort_tumour_umi)

sort_tumour_size <- ggplot(sort_tumours, aes(cells, tumour_size, colour = Animal_ID)) +
  geom_point(size = 3) +
  facet_wrap(~population) +
  labs(y = expression(paste("tumour size ", mm^2)))

ggsave("sort_tumour_size.png", sort_tumour_size)

sort_tumours %>% group_by(population) %>% summarise(mean = mean(cells), max = max(cells), min = min(cells), median = median(cells))

#####
# repertoire yield model
#####

raw_train <- read_xlsx("UMI_yield_train.xlsx")

raw_seq <- read_xlsx("raw_seq.xlsx") %>% gather(exp, reads, - population) %>% drop_na()

yield_train <- raw_train %>% select(-3, -4, -7, -8) %>% gather(exp, UMIs, -cells, -population) %>% drop_na() 

flank_col <- c(rep("L", times = 3), rep("R", times = 2), rep(c("L", "R"), each = 3, times = 3), rep(NA, times = 5))

flank_col2 <- c(rep(c("L", "R"), each = 3, length.out = 23), rep(NA, times = 5))

yield_train <- cbind(yield_train, flank_col2)

yield_train$exp <- str_remove(yield_train$exp, "_UMIs")

raw_seq$exp <- tolower(raw_seq$exp) %>% str_replace("dln", "LN")

raw_seq <- cbind(raw_seq, flank_col)

raw_seq <- arrange(raw_seq, exp)

yield_train <- arrange(yield_train, exp)

yield_train <- cbind(yield_train, raw_seq$reads)

ytn <- c("population", "cells", "exp", "UMIs", "flank", "reads")

names(yield_train) <- ytn

H_read_mod <- lm(UMIs ~ cells, data = yield_train[yield_train$reads > 500000,])

L_read_mod <- lm(UMIs ~ cells, data = yield_train[yield_train$reads < 500000,])

new <- data.frame(cells = sort_tumours$cells)

sort_tumours$pred_UMIs_high <- predict(H_read_mod, newdata = new)

sort_tumours$pred_UMIs_low <- predict(L_read_mod, newdata = new)

sort_tumours <- sort_tumours %>% gather(model, UMIs, pred_UMIs_high, pred_UMIs_low)

yield_train$type <- rep("training", times = 28)

sort_tumours$type <- rep("test", times = 32)

sort_plot <- rbind(sort_tumours[, c(3, 16, 15)], yield_train[, c(2, 5, 8)])

sort_predict <- ggplot(sort_plot, aes(cells, UMIs, colour = type)) + 
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0, show.legend = T) +
  theme(legend.position = c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1))

ggsave("yield_model.png", sort_predict)

yield_train <- yield_train %>% separate(exp, into = c("exp", "tissue"), sep = "_")

exp_yield <- yield_train[yield_train$exp == c("1.1-1", "1.3-4"),]

exp_yield$exp <- factor(exp_yield$exp)

normalise <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

exp_yield <- exp_yield %>% gather(param, normal_value, reads, UMIs, cells)

ggplot(exp_yield, aes(exp, normal_value, colour = param)) + geom_boxplot() + labs(y = "min max scale")

exp_yield[,c(2,5,7)] <- map(exp_yield[,c(2,5,7)], normalise)

model_yield <- rpart(reads ~ exp + UMIs, data = exp_yield)

rpart.plot(model_yield, box.palette = "blue")

reads_split <- yield_train %>% 
  filter(exp == c("1.1-1", "1.3-4")) %>% 
  ggplot(aes(exp, reads)) + 
  geom_point() +
  geom_hline(yintercept = 655810) +
  scale_y_log10()

ggsave("model_split.png", reads_split)

#####
