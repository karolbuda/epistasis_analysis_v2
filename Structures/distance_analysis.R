library(tidyverse)

full_set = tibble()

for(dir in 1:length(list.dirs(path = "Structures")[-1])) {
  cur_set = try(read_csv(paste0(list.dirs(path = "Structures")[-1][dir], "/output.txt"), col_names = F), silent = T)
  enzyme = str_remove(list.dirs(path = "Structures")[-1][dir], "Structures/")
  if(class(cur_set) != "try-error") {
    cur_set$enzyme = enzyme
    full_set = bind_rows(full_set, cur_set)
    enzyme = c()
  }
}

#full_set = bind_rows(full_set, cur_set)

# Import 3hbr distances for the 3 trajectories

cur_set = read_csv(paste0("Structures/OXA-48/3hbr_solo", "/trajectory1.txt"), col_names = F)
cur_set$enzyme = "OXA-48/traj1"

full_set = bind_rows(full_set, cur_set)

cur_set = read_csv(paste0("Structures/OXA-48/3hbr_solo", "/trajectory2.txt"), col_names = F)
cur_set$enzyme = "OXA-48/traj2"

full_set = bind_rows(full_set, cur_set)

cur_set = read_csv(paste0("Structures/OXA-48/3hbr_solo", "/trajectory3.txt"), col_names = F)
cur_set$enzyme = "OXA-48/traj3"

full_set = bind_rows(full_set, cur_set)

# Import mira and weinreich TEM

cur_set = read_csv(paste0("Structures/TEM", "/mira.txt"), col_names = F)
cur_set$enzyme = "TEM_MIRA"

full_set = bind_rows(full_set, cur_set)

cur_set = read_csv(paste0("Structures/TEM", "/weinreich.txt"), col_names = F)
cur_set$enzyme = "TEM_WEINREICH"

full_set = bind_rows(full_set, cur_set)


##############

colnames(full_set) = c("res1", "res2", "dist", "enzyme") 

### What does the distance distribution look like

full_set %>%
  mutate(dist = as.numeric(dist)) %>%
  ggplot(aes(x = dist)) +
  geom_histogram(binwidth = 1, fill = "#edae49") +
  facet_wrap(~ enzyme) +
  theme_bw()


##### SNIPPET OF HIGHER ORDER CODE

library(ggExtra)
library(ggrepel)
library(ggpubr)

higher_order_list = list.files(path = "Output", pattern="higher_box", recursive = T)

higher_df = data.frame()
for(higher_file in 1:length(higher_order_list)){
  Condition = str_split(str_split(higher_order_list[higher_file], "_")[[1]][3], "/")[[1]][1]
  Measurement = str_split(higher_order_list[higher_file], "_")[[1]][2]
  Enzyme = str_split(higher_order_list[higher_file], "_")[[1]][1]
  appending_file = read_csv(paste0("Output/", higher_order_list[higher_file], collapse = ""))
  appending_file$Condition = Condition
  appending_file$Measurement = Measurement
  appending_file$Enzyme = Enzyme
  higher_df = bind_rows(higher_df, appending_file)
}


higher_df = higher_df %>%
  unite("unique_id", c(pos, Condition, Measurement, Enzyme), remove = F)

higher_all_sign_check = higher_df %>%
  mutate(mutations = factor(paste("Order", mutations, sep = " "))) %>%
  group_by(unique_id, mutations) %>%
  summarise(min_effect = min(avg),
            max_effect = max(avg))


##########

full_set$combo = paste0("p", full_set$res1, "|", "p", full_set$res2)

full_set = distinct(full_set, combo, .keep_all = T)

# gets average effects of all pairwise combos

higher_df_twos = higher_df %>% 
  filter(mutations == 2) %>%
  group_by(unique_id) %>%
  summarise(avg_effect = mean(avg),
            idio_index = 2*sd(avg),
            pos = pos[1],
            enzyme = Enzyme[1])

distances = c()

for(i in 1:length(higher_df_twos$pos)) {
  if(sum(full_set$combo %in% higher_df_twos$pos[i]) > 0) {
    location = which(full_set$combo %in% higher_df_twos$pos[i])
    distances = c(distances, full_set[location,]$dist)
  } else {
    distances = c(distances, NA)
  }
}

higher_df_twos$dist = distances

## There doesn't necessarily seem to be a big correlation between effect and distance

higher_df_twos %>%
  drop_na() %>%
  ggplot(aes(x = dist, y = avg_effect)) +
  geom_point()

## No strong correlation between idiosyncrasy at a position and distance between pairs

higher_df_twos %>%
  drop_na() %>%
  ggplot(aes(x = dist, y = idio_index)) +
  geom_point()

#### FOR WT no strong correlation between distance of interaction and epistatic contribution

higher_df_twos_wt = higher_df[which(str_count(higher_df$genotype, "1") == 0), ] %>%
  filter(mutations == 2)
  
distances = c()

for(i in 1:length(higher_df_twos_wt$pos)) {
  if(sum(full_set$combo %in% higher_df_twos_wt$pos[i]) > 0) {
    location = which(full_set$combo %in% higher_df_twos_wt$pos[i])
    distances = c(distances, full_set[location,]$dist)
  } else {
    distances = c(distances, NA)
  }
}

higher_df_twos_wt$dist = distances

# Defined adaptive_trajectories from statistical_analysis.Rmd

adaptive_trajectories = c("DHFR_ki_trajg", 
                          "DHFR_ki_trajr",
                          "MPH_catact_ZnPTM",
                          "NfsA_ec50_2039",
                          "NfsA_ec50_3637",
                          "OXA-48_ic50_CAZtraj1",
                          "OXA-48_ic50_CAZtraj2",
                          "OXA-48_ic50_CAZtraj3",
                          "PTE_catact_2NH",
                          "TEM_MIC_weinreich",
                          "DHFR_ic75_palmer")

higher_df_twos_wt %>%
  unite("partial_id", c(Enzyme, Measurement, Condition), sep = "_", remove = F) %>%
  filter(partial_id %in% adaptive_trajectories) %>%
  drop_na() %>%
  ggplot(aes(x = dist, y = abs(avg), color = partial_id)) +
  geom_point() +
  labs(x = "Distance between alpha carbon pairs", 
       y = expression(epsilon*" between pairs")) +
  theme_classic()


###### Playing with specific enzymes


higher_df_twos_plot = higher_df_twos %>%
  drop_na() %>%
  ggplot(aes(x = dist, y = abs(avg_effect))) +
  geom_point() +
  stat_cor(aes(label = paste(..r.label.., sep = "~~")), method = "spearman", cor.coef.name = "rho", label.x.npc = "middle",
           label.y.npc = "top", r.accuracy = 0.01, size = 2.5, color = "black") +
  geom_smooth(method = 'loess', formula = y ~ x) +
  labs(x = "Distance between Positions (\u00c5)",
       y = expression("Pairwise Epistasis Absolute Value (|"*epsilon*"|)")) +
  facet_wrap(~ enzyme, scales = "free") +
  theme_classic() +
  theme(text = element_text(size=9), axis.text = element_text(size=8, color = "black"),
        legend.position = "none", axis.line=element_line())

higher_df_twos_plot

# Uncomment this if you want to save the plot

#ggsave("supp_fig_distance.svg", higher_df_twos_plot, width = 180, height = 247/2, dpi = 300, units = "mm")

