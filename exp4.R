
# data
#####

factors = c("jkexp",      "mouse",      "tissue",     "flank",      "population")

working_path <- file.path(getwd(), "experiments")

fig1 <- bind_rows(read.csv(file.path(working_path, "4_1_0", "tumour_growth.csv"), stringsAsFactors = F),
                   read.csv(file.path(working_path, "4_2_0", "tumour_growth.csv"), stringsAsFactors = F)) %>%
  unite("exp", Expt.ID, mouse, flank, sep = "_", remove = F)

jk4 <- bind_rows(jk41, jk42) %>%
  factor_extractor()

jk4_sort <- list(jk41 = read.csv("D:/data/experiments/4_1_0/sort_data.csv", stringsAsFactors = F),
                 jk42 = read.csv("D:/data/experiments/4_2_0/sort_data.csv", stringsAsFactors = F)) %>%
  map(factor_extractor) %>%
  bind_rows(.id = "treated") %>%
  select(-X) %>%
  filter(tissue == "T")

jk4_summary <- list(jk41 = read.csv("D:/data/experiments/4_1_0/summary_stats.csv"),
                    jk42 = read.csv("D:/data/experiments/4_2_0/summary_stats.csv")) %>%
  map(factor_extractor) %>%
  bind_rows(.id = "treated") %>%
  mutate(population = if_else(str_detect(population, "4"), "CD4", "CD8"))
#####
# figure 1
#####

sequenced_only <- jk4_summary %>%
  select(jkexp, mouse, flank) %>%
  mutate(mouse = str_replace(.$mouse, fixed("."), "-")) %>%
  unite("exp", sep = "_")

fig1A <- fig1 %>%
  right_join(sequenced_only) %>%
  group_by(day, Expt.ID) %>%
  mutate(mean_tg = mean(tumour_size),
         sd_tg = sd(tumour_size))

table(fig1A$Expt.ID, fig1A$mouse)

fig1A$Expt.ID[fig1A$Expt.ID == "JK4.1"] <- "control"

fig1A$Expt.ID[fig1A$Expt.ID == "JK4.2"] <- "OX-40 + CTLA-4"

ggplot(fig1A, aes(day, tumour_size, group = exp, colour = Expt.ID)) +
  geom_errorbar(aes(ymin = mean_tg - sd_tg, 
                    ymax = mean_tg + sd_tg, 
                    group = interaction(day, Expt.ID))) +
  geom_smooth(aes(group = Expt.ID)) +
  labs(y = "tumour size mm^2", x = "days post inoculation", title = "A") +
  geom_vline(xintercept = 10, linetype = 2) +
  geom_vline(xintercept = 13, linetype = 2) +
  scale_x_continuous(breaks = seq(from = min(fig1A$day), to = max(fig1A$day), by = 2)) +
  scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10),
                     labels = c(0, "", 20, "", 40, "", 60)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.line = element_line())

fig1A %>%
  filter(day == 14) %>%
  ggplot(aes(Expt.ID, tumour_size, fill = Expt.ID)) +
  geom_boxplot() +
  stat_compare_means(aes(label = ..p.signif..), label.x = 1.5) +
  labs(y = "tumour size mm^2 at day 14", title = "B") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank())

fig1B <- jk4_sort %>%
  select(-exp) %>%
  mutate(mouse = str_replace(.$mouse, fixed("."), "-")) %>%
  unite("exp", jkexp, mouse, flank, sep = "_", remove = F) %>%
  right_join(sequenced_only)

#####
# figure 2
#####

jk4_8 <- morisita_df(jk4, "CD8")

jk4_m <- jk4_8 %>% 
  separate(from, into = c("exp1", "mouse1", NA, NA, NA), remove = F, sep = "_") %>%
  separate(to, into = c("exp2", "mouse2", NA, NA, NA), remove = F, sep = "_") %>%
  mutate(same_mouse = if_else(exp1 == exp2 & mouse1 == mouse2, T, F),
         same_exp = if_else(exp1 == exp2, T, F),
         population = if_else(str_detect(.$from, "CD4"), "CD4", "CD8"))

index_lim <- 0.15

jk4.1_4_net <- jk4_m %>%
  filter(exp1 == "JK4.1", population == "CD4", same_exp, index > index_lim) %>%
  select(from, to, index) %>%
  graph_from_data_frame(vertices = jk4_m %>%
                          filter(exp1 == "JK4.1", population == "CD4", same_exp, index > index_lim) %>%
                          gather("tmp", "mouse", from, to) %>%
                          select("exp" = mouse) %>%
                          distinct() %>%
                          factor_extractor())

jk4.1_4_grob <- ggnet2(jk4.1_4_net, 
       color = V(jk4.1_4_net)$mouse, 
       shape = V(jk4.1_4_net)$flank,
       color.palette = "Dark2",
       size = 7) +
  guides(shape = "none") +
  theme(legend.position = "bottom",
        panel.background = element_rect(colour = "black"))

jk4.1_8_net <- jk4_m %>%
  filter(exp1 == "JK4.1", population == "CD8", same_exp, index > index_lim) %>%
  select(from, to, index) %>%
  graph_from_data_frame(vertices = jk4_m %>%
                          filter(exp1 == "JK4.1", population == "CD8", same_exp, index > index_lim) %>%
                          gather("tmp", "mouse", from, to) %>%
                          select("exp" = mouse) %>%
                          distinct() %>%
                          factor_extractor())

jk4.1_8_grob <- ggnet2(jk4.1_8_net, color = V(jk4.1_8_net)$mouse, shape = V(jk4.1_8_net)$flank,
       color.palette = "Dark2",
       size = 7) +
  guides(shape = "none") +
  theme(legend.position = "bottom",
        panel.background = element_rect(colour = "black"))

jk4.2_4_net <- jk4_m %>%
  filter(exp1 == "JK4.2", population == "CD4", same_exp, index > index_lim) %>%
  select(from, to, index) %>%
  graph_from_data_frame(vertices = jk4_m %>%
                          filter(exp1 == "JK4.2", population == "CD4", same_exp, index > index_lim) %>%
                          gather("tmp", "mouse", from, to) %>%
                          select("exp" = mouse) %>%
                          distinct() %>%
                          factor_extractor())

jk4.2_4_grob <- ggnet2(jk4.2_4_net, color = V(jk4.2_4_net)$mouse, shape = V(jk4.2_4_net)$flank,
       color.palette = "Dark2",
       size = 7) +
  guides(shape = "none") +
  theme(legend.position = "top",
        panel.background = element_rect(colour = "black"))

jk4.2_8_net <- jk4_m %>%
  filter(exp1 == "JK4.2", population == "CD8", same_exp, index > index_lim) %>%
  select(from, to, index) %>%
  graph_from_data_frame(vertices = jk4_m %>%
                          filter(exp1 == "JK4.2", population == "CD8", same_exp, index > index_lim) %>%
                          gather("tmp", "mouse", from, to) %>%
                          select("exp" = mouse) %>%
                          distinct() %>%
                          factor_extractor())

jk4.2_8_grob <- ggnet2(jk4.2_8_net, color = V(jk4.2_8_net)$mouse, shape = V(jk4.2_8_net)$flank,
       color.palette = "Dark2",
       size = 7) +
  guides(shape = "none") +
  theme(legend.position = "top",
        panel.background = element_rect(colour = "black"))

l_matrix <- rbind(c(1, 1, 2, 2), 
            c(1, 1, 2, 2))

grid.arrange(fig1_b, 
             grid.arrange(jk4.2_4_grob, 
             jk4.2_8_grob, 
             jk4.1_4_grob, 
             jk4.1_8_grob, 
             top = text_grob("C", x = 0, hjust = 0, face = "bold")),
             layout_matrix = l_matrix,
             ncol = 2)

jk4_m_anti <- jk4_m %>%
  filter(same_mouse == F & same_exp == T) %>%
  anti_join(jk4_m %>% filter(exp1 == "JK4.1", population == "CD4", same_mouse == F, same_exp == T) %>% dplyr::sample_n(54)) %>%
  anti_join(jk4_m %>% filter(exp1 == "JK4.1", population == "CD8", same_mouse == F & same_exp == T) %>% dplyr::sample_n(54)) %>%
  anti_join(jk4_m %>% filter(exp1 == "JK4.2", population == "CD4", same_mouse == F & same_exp == T) %>% dplyr::sample_n(54)) %>%
  anti_join(jk4_m %>% filter(exp1 == "JK4.2", population == "CD8", same_mouse == F & same_exp == T) %>% dplyr::sample_n(54))

jk4_mf <- jk4_m %>%
  filter(same_mouse == T) %>%
  bind_rows(jk4_m_anti) %>%
  select(from, to, index, same_mouse) %>%
  separate(from, into = c("exp", NA, NA, NA, "population"), remove = F, sep = "_")

jk4_mf$exp <- if_else(jk4_mf$exp == "JK4.1", "control", "aCTLA-4 + aOX-40")

fig1_b <- ggplot(jk4_mf, aes(same_mouse, index, fill = same_mouse)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.1) +
  facet_grid(exp ~ population) +
  stat_compare_means(label.x = 1.5, label = "p.signif") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(labels = c("between", "within"),
                      values = c("slateblue", "chartreuse3")) +
  labs(title = "B")

#####




hd_dual_tumour(jk4)




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

jk4_intersects <- all_intersects(jk4)

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

jk4_8fj <- all_full_joins_df(jk4, "CD8") %>%
  factor_extractor_union() %>%
  filter(JKexp.x == JKexp.y)

jk4_8fj$n_mice <- jk4_8fj$mouse.x == jk4_8fj$mouse.y & jk4_8fj$JKexp.x == jk4_8fj$JKexp.y

jk4_8fj$n_mice[jk4_8fj$n_mice == T] <- "within"

jk4_8fj$n_mice[jk4_8fj$n_mice == F] <- "between"


sum(jk4_8fj$n_mice == "between")

sum(jk4_8fj$n_mice == "within")

jk4_8fj_sample <- jk4_8fj %>%
  filter(JKexp.x == "JK4.1") %>%
  ungroup() %>%
  sample_n(size = sum(jk4_8fj$JKexp.x == "JK4.1") - sum(jk4_8fj$JKexp.x == "JK4.2"))

jk4_8fj_sample2 <- jk4_8fj %>%
  anti_join(jk4_8fj_sample)

jk4_8fj_sample3 <- jk4_8fj_sample2 %>%
  filter(n_mice == "between") %>%
  ungroup() %>%
  sample_n(size = sum(.$n_mice == "between") - sum(jk4_8fj_sample2$n_mice == "within"))

jk4_8fj_top <- jk4_8fj_sample2 %>%
  anti_join(jk4_8fj_sample3) %>%
  filter(rank.x < 50 | rank.y < 50) %>%
  na.omit()

jk4_8fj_s <- jk4_8fj_sample2 %>%
  anti_join(jk4_8fj_sample3) %>%
  na.omit()

ggplot(jk4_8fj_top, aes(rank.x, rank.y, colour = n_mice)) +
  geom_jitter(width = 0.3, height = 0.3) +
  scale_y_log10() +
  scale_x_log10() +
  facet_grid(JKexp.x ~ JKexp.y)

ggplot(jk4_8fj_top, aes(rank.x, rank.y, colour = JKexp.x, alpha = 0.3)) +
  geom_jitter(width = 0.2, height = 0.2) +
  scale_y_log10() +
  scale_x_log10() +
  facet_grid(JKexp.x ~ n_mice)

ggplot(jk4_8fj_s, aes(n_mice, y = n_mice, group = JKexp.x, fill = JKexp.x)) +
  geom_bar(position = "dodge", aes(y = stat(count))) +
  labs(y = "count")

ggplot(jk4_8fj_s, aes(x = JKexp.x, y = n_mice)) +
  geom_point(aes(y = stat(count)))

#####
# removing similar clones
#####



#####
