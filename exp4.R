jk4 <- list(jk41 = read.csv("D:/data/experiments/4_1_0/cleaned_CDR3s.csv"),
            jk42 = read.csv("D:/data/experiments/4_2_0/cleaned_CDR3s.csv")) %>%
  map(data_aa) %>%
  bind_rows()

hd_dual_tumour(jk4)

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

jk4_unions <- jk4 %>%
  bind_rows() %>%
  group_by(exp) %>%
  arrange(desc(PID.fraction_sum), .by_group = T) %>%
  group_split()

for (i in 1:length(jk4_unions)) {
  jk4_unions[[i]]$rank <- 1:length(jk4_unions[[i]]$aaSeqCDR3)
}

jk4_unions_2 <- list()

for (i in 1:length(jk4_unions)) {
  jk4_unions_2[[i]] <- jk4_unions[[i]][1:2,]
}

jk4_unions_2 <- bind_rows(jk4_unions_2) %>%
  factor_extractor() %>%
  group_by(JKexp, mouse, pop) %>%
  unite("exp_m_p", JKexp, mouse, pop, sep = "-", remove = F)

table(jk4_unions_1$aaSeqCDR3, jk4_unions_1$exp_m_p)

ggplot(jk4_unions_1, aes(aaSeqCDR3, fill = mouse)) +
  geom_histogram(stat = "count") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(. ~ JKexp, scales = "free_y") +
  scale_color_brewer(type = "div") +
  geom_hline(yintercept = 2)

jk4_unions_20 <- jk4_unions %>%
  bind_rows() %>%
  factor_extractor() %>%
  filter(rank <= 20, JKexp == "JK4.2", mouse == "1.3", pop == "CD8") %>%
  arrange(desc(PID.fraction_sum)) %>%
  group_by(aaSeqCDR3) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  mutate(aaSeqCDR3 = factor(aaSeqCDR3, ordered = T))

jk4_unions_20 %>%
  ggplot(aes(aaSeqCDR3, rank, fill = flank, group = flank)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_blank(),
        panel.grid = element_line(colour = "grey"),
        legend.position = "none") +
  facet_wrap(~ aaSeqCDR3, scales = "free_x")

jk4_fj <- all_full_joins_df(jk4, "CD8")  
