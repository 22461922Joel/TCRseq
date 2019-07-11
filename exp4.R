jk4 <- list(jk41 = read.csv("D:/data/experiments/4_1_0/cleaned_CDR3s.csv"),
            jk42 = read.csv("D:/data/experiments/4_2_0/cleaned_CDR3s.csv")) %>%
  map(data_aa) %>%
  map(factor_extractor)

jk4_sort <- list(jk41 = read.csv("D:/data/experiments/4_1_0/sort_data.csv"),
                 jk42 = read.csv("D:/data/experiments/4_2_0/sort_data.csv")) %>%
  map(factor_extractor) %>%
  bind_rows(.id = "treated") %>%
  select(-X)

jk4_summary <- list(jk41 = read.csv("D:/data/experiments/4_1_0/summary_stats.csv"),
                    jk42 = read.csv("D:/data/experiments/4_2_0/summary_stats.csv")) %>%
  map(factor_extractor) %>%
  bind_rows(.id = "treated")

jk4_ss <- jk4_sort %>%
  left_join(jk4_summary) %>%
  na.omit()

jk4_ss$treated[jk4_ss$treated == "jk41"] <- "ctl"

jk4_ss$treated[jk4_ss$treated == "jk42"] <- "combo"

ggplot(jk4_ss, aes(pop, cells, fill = treated)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(label = ..p.signif..)) +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  labs(y = "sorted cells")

ggplot(jk4_ss, aes(evenness, cells, colour = treated)) +
  geom_point() +
  facet_grid(. ~ pop) +
  labs(y = "sorted cells", x = "total clones") +
  geom_abline(slope = 1, intercept = 0)

ggplot(jk4_ss, aes(richness, cells, colour = treated)) +
  geom_point() +
  facet_grid(. ~ pop) +
  labs(y = "sorted cells", x = "unique clones") +
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10()

ggplot(jk4_ss, aes(richness, evenness, colour = treated)) +
  geom_point() +
  facet_grid(. ~ pop) +
  labs(y = "total clones", x = "unique clones") +
  geom_abline(slope = 1, intercept = 0) +
  theme(legend.title = element_blank()) +
  scale_x_log10(limits = c(150, max(jk4_ss$richness)))+
  scale_y_log10(limits = c(150, max(jk4_ss$evenness)))

ggplot(jk4_ss, aes(pop, diversity, fill = treated)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(label = ..p.signif..)) +
  theme(legend.title = element_blank(), axis.title.x = element_blank())

jk4_intersects <- map(jk4, all_intersects)

jk4_unions

pairing <- function(df) {
  df %>%
    filter(n_mice == 2) %>%
    separate(mice, 
             into = c(letters[1:2], rep(NA, times = (str_count(.$mice[1], "-") - 1))), 
             sep = "-") %>%
    separate(a, into = c("jkexpa", "mousea", "tissuea", "flanka", "popa"), sep = "_", remove = F) %>%
    separate(b, into = c("jkexpb", "mouseb", "tissueb", "flankb", "popb"), sep = "_", remove = F) %>%
    mutate(same_mouse = mousea == mouseb)
}

jk4_pairs <- map(jk4_intersects, pairing) %>%
  bind_rows(.id = "treated")

jk4_pairs$same_mouse[jk4_pairs$same_mouse == F] <- "between"

jk4_pairs$same_mouse[jk4_pairs$same_mouse == T] <- "within"

jk4_pairs$treated[jk4_pairs$treated == "jk41"] <- "ctl"

jk4_pairs$treated[jk4_pairs$treated == "jk42"] <- "combo"

ggplot(jk4_pairs, aes(same_mouse, n_clones, fill = treated)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  facet_grid(. ~ pop.x) +
  stat_compare_means(aes(label = ..p.signif..)) +
  theme(legend.title = element_blank(), axis.title.x = element_blank())

ggplot(jk4_pairs, aes(treated, clonal_proportion, fill = same_mouse)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(label = ..p.signif..)) +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  facet_grid(. ~ pop.x)



