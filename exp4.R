jk4 <- list(jk41 = read.csv("D:/data/experiments/4_1_0/cleaned_CDR3s.csv"),
            jk42 = read.csv("D:/data/experiments/4_2_0/cleaned_CDR3s.csv")) %>%
  map(data_aa) %>%
  map(factor_extractor)

jk4_overlaps <- map(jk4, all_overlaps)

jk4_unions

pairing <- function(df) {
  df %>%
    filter(n_mice == 2) %>%
    separate(mice, 
             into = c(letters[1:2], rep(NA, times = (str_count(.$mice[1], "-") - 1))), 
             sep = "-") %>%
    mutate(same_mouse = stringdist(.$a, .$b))
}

jk4_pairs <- map(jk4_overlaps, pairing) %>%
  bind_rows(.id = "treated")

ggplot(jk4_pairs, aes(as.character(same_mouse), n_clones, fill = treated)) +
  geom_boxplot() +
  labs(x = "single mouse, different mice/same flank, different mice/different flank") +
  facet_grid(. ~ pop.x)

ggplot(jk4_pairs, aes(as.character(same_mouse), clonal_proportion, fill = treated)) +
  geom_boxplot() +
  labs(x = "single mouse, different mice/same flank, different mice/different flank") +
  facet_grid(. ~ pop.x)



