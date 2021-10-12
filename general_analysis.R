setwd("C:/Users/karol/OneDrive - The University Of British Columbia (1)/PhD/Epistasis_Lit/Bulk/Output")

library(tidyverse)
library(ggExtra)
library(ggrepel)

### FUNCTION FOR THE HALF VIOLIN

`%notin%` <- Negate(`%in%`)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x - width / 2,
                     xmax = x)
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, 
                              xmaxv = x,
                              xminv = x + violinwidth * (xmin - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

## Reads all pred_dfs to look at epistasis from WT reference

all_ratios = do.call("rbind", lapply(list.files(pattern="ratio_export.csv", recursive = T), read_csv, skip = 1, col_names = c("id", "effect", "muts", "cv", "colors")))

all_ratios = all_ratios %>% filter(muts > 1)

length(all_ratios$effect) # examined muts

sum(all_ratios$effect > 1.5) # synergystic
paste("Synergystic", sum(all_ratios$effect > 1.5)/length(all_ratios$effect)*100, "%", collapse = "") # synergystic %

sum(all_ratios$effect >= 1/1.5 & all_ratios$effect <= 1.5) # antagonistic
paste("Neutral", sum(all_ratios$effect >= 1/1.5 & all_ratios$effect <= 1.5)/length(all_ratios$effect)*100, "%", collapse = "") # antagonistic %

sum(all_ratios$effect < 1/1.5) # antagonistic
paste("Antagonistic", sum(all_ratios$effect < 1/1.5)/length(all_ratios$effect)*100, "%", collapse = "") # antagonistic %

# Quick histogram

all_ratios %>%
  ggplot(aes(x = log10(effect))) +
  geom_histogram(bins = 100, alpha = 0.8) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  geom_vline(xintercept = log10(1/1.5), lty = 2) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

## 2nd step only

all_ratios_2 = all_ratios %>% filter(muts == 2)

length(all_ratios_2$effect) # examined muts

sum(all_ratios_2$effect > 1.5) # synergystic
paste("Synergystic", sum(all_ratios_2$effect > 1.5)/length(all_ratios_2$effect)*100, "%", collapse = "") # synergystic %

sum(all_ratios_2$effect >= 1/1.5 & all_ratios_2$effect <= 1.5) # antagonistic
paste("Neutral", sum(all_ratios_2$effect >= 1/1.5 & all_ratios_2$effect <= 1.5)/length(all_ratios_2$effect)*100, "%", collapse = "") # antagonistic %

sum(all_ratios_2$effect < 1/1.5) # antagonistic
paste("Antagonistic", sum(all_ratios_2$effect < 1/1.5)/length(all_ratios_2$effect)*100, "%", collapse = "") # antagonistic %


## 3rd step only

all_ratios_3 = all_ratios %>% filter(muts == 3)

length(all_ratios_3$effect) # examined muts

sum(all_ratios_3$effect > 1.5) # synergystic
paste("Synergystic", sum(all_ratios_3$effect > 1.5)/length(all_ratios_3$effect)*100, "%", collapse = "") # synergystic %

sum(all_ratios_3$effect >= 1/1.5 & all_ratios_3$effect <= 1.5) # antagonistic
paste("Neutral", sum(all_ratios_3$effect >= 1/1.5 & all_ratios_3$effect <= 1.5)/length(all_ratios_3$effect)*100, "%", collapse = "") # antagonistic %

sum(all_ratios_3$effect < 1/1.5) # antagonistic
paste("Antagonistic", sum(all_ratios_3$effect < 1/1.5)/length(all_ratios_3$effect)*100, "%", collapse = "") # antagonistic %


## 4th step only

all_ratios_4 = all_ratios %>% filter(muts == 4)

length(all_ratios_4$effect) # examined muts

sum(all_ratios_4$effect > 1.5) # synergystic
paste("Synergystic", sum(all_ratios_4$effect > 1.5)/length(all_ratios_4$effect)*100, "%", collapse = "") # synergystic %

sum(all_ratios_4$effect >= 1/1.5 & all_ratios_4$effect <= 1.5) # antagonistic
paste("Neutral", sum(all_ratios_4$effect >= 1/1.5 & all_ratios_4$effect <= 1.5)/length(all_ratios_4$effect)*100, "%", collapse = "") # antagonistic %

sum(all_ratios_4$effect < 1/1.5) # antagonistic
paste("Antagonistic", sum(all_ratios_4$effect < 1/1.5)/length(all_ratios_4$effect)*100, "%", collapse = "") # antagonistic %

## Reads all csvs for raw fitness landscapes effects

fit_land = do.call("rbind", lapply(list.files(pattern="observed_values.csv", recursive = T), read_csv, skip = 1, col_names = c("id", "effect", "muts", "2", "3", "enz")))

fit_land = fit_land %>% select(-c('2', '3'))

# neg, neutral, positive fold-changes

c(sum(fit_land$effect < log10(1/1.5)), sum(fit_land$effect >= log10(1/1.5) & fit_land$effect <= log10(1.5)), sum(fit_land$effect > log10(1.5))) / dim(fit_land)[1]

fit_land %>%
  ggplot(aes(x = effect, fill = enz)) +
  geom_histogram(bins = 100, alpha = 0.8) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  geom_vline(xintercept = log10(1/1.5), lty = 2) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

fit_land %>%
  ggplot(aes(x = effect)) +
  geom_histogram(bins = 100, alpha = 0.8) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  geom_vline(xintercept = log10(1/1.5), lty = 2) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# density plot

fit_land %>%
  ggplot(aes(x = effect, fill = enz)) +
  geom_density(alpha = 0.8) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  geom_vline(xintercept = log10(1/1.5), lty = 2) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

## Reads all csvs in collated folder and combines them

d = do.call("rbind", lapply(list.files(pattern="*2d_box.csv", recursive = T), read_csv, col_names = c("positions", "identity", "mut", "effects", "enzyme", "type", "cond")))


d1 = d %>% filter_all(any_vars(!is.na(.))) # removes rows with all NA before making unique ID

d1 = d1 %>%
  unite("unique_id", c(positions, enzyme, type, cond), remove = F)


## Histogram like above but for d1

d1 %>%
  ggplot(aes(x = effects, fill = enzyme)) +
  geom_histogram(bins = 100, alpha = 0.8) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  geom_vline(xintercept = log10(1/1.5), lty = 2) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

### General Info

# Sig Neg, Neutral, Sig Pos

general = c(sum(d1$effects < log10(1/1.5)) / length(d1$effects),
  sum(d1$effects <= log10(1.5) & d1$effects >= log10(1/1.5)) / length(d1$effects),
  sum(d1$effects > log10(1.5)) / length(d1$effects)
)

general

sum(general)


##

d1 %>%
  mutate(enzyme = factor(enzyme)) %>%
  ggplot(aes(x = unique_id, y = effects, fill = enzyme)) +
  geom_boxplot(coef = 10) +
  geom_jitter() +
  geom_hline(yintercept = 0, lty = 2,) +
  geom_hline(yintercept = -log(1.5, 10), lty = 2, col = "red") +
  geom_hline(yintercept = log(1.5, 10), lty = 2, col = "red") +
  facet_wrap(~ enzyme, scale = "free_x") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



#####

aic_list = list.files(pattern="aic", recursive = T)

aic_df = data.frame()
for(aic_file in 1:length(aic_list)){
  Condition = str_split(str_split(aic_list[aic_file], "_")[[1]][3], "/")[[1]][1]
  Measurement = str_split(aic_list[aic_file], "_")[[1]][2]
  appending_file = read_csv(aic_list[aic_file])
  appending_file$Condition = Condition
  appending_file$Measurement = Measurement
  aic_df = bind_rows(aic_df, appending_file)
}

aic_df %>%
  mutate(Step = factor(Step, levels = unique(Step))) %>%
  ggplot(aes(x = Step, y = MAE)) +
  geom_point(size = 3.5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  stat_summary(fun=mean, colour="darkred", geom="crossbar", width=0.2) +
  scale_y_continuous(limits  = c(0, 8)) +
  theme_classic() +
  labs(x = "Mutational Trajectory Step",
       y = "Mean Absolute Error of Prediction") +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))


### Exploring idiosyncrasies in MAE

aic_df %>%
  mutate(Step = factor(Step, levels = unique(Step))) %>%
  filter(Enzyme == "AP" | Enzyme == "PTE") %>%
  ggplot(aes(x = Step, y = MAE, color = Enzyme)) +
  geom_point(size = 3.5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  scale_y_continuous(limits  = c(0, 8)) +
  theme_classic() +
  labs(x = "Mutational Trajectory Step",
       y = "Mean Absolute Error of Prediction") +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

aic_df %>%
  mutate(Step = factor(Step, levels = unique(Step))) %>%
  filter(Enzyme == "OXA-48") %>%
  ggplot(aes(x = Step, y = MAE, group = Condition)) +
  geom_point(size = 3.5) +
  stat_summary(fun=sum, geom="line") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  #scale_y_continuous(limits  = c(0, 8)) +
  theme_classic() +
  labs(x = "Mutational Trajectory Step",
       y = "Mean Absolute Error of Prediction") +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))



#####

sign_unique_ids = c()
for(i in 1:length(unique(d1$unique_id))) {
  testing_effects = d1 %>%
    filter(unique_id == unique(d1$unique_id)[i]) %>%
    pull(effects)
  
  # Checks whether any effects are less than lower significance threshold AND upper significance threshold
  if(any(testing_effects > log(1.5, 10)) & any(testing_effects < log(1/1.5, 10))) {
    sign_unique_ids = c(sign_unique_ids, unique(d1$unique_id)[i])
  }
  
}

# Add TRUE star argument for those with sign

d1$sign = d1$unique_id %in% sign_unique_ids

# Add TRUE where its WT

d1$wt = d1$mut == 1

all_sign_check = d1 %>%
  group_by(unique_id) %>%
  summarise(min_effect = min(effects),
            max_effect = max(effects))

### Sign sum
sign_sum = 0

threshold = 1.5

# All
length(all_sign_check$max_effect)

# Only Neutral
sum(all_sign_check$min_effect >= log10(1/threshold) & all_sign_check$max_effect <= log10(threshold))
sum(all_sign_check$min_effect >= log10(1/threshold) & all_sign_check$max_effect <= log10(threshold)) / length(all_sign_check$max_effect)

sign_sum = sign_sum + sum(all_sign_check$min_effect >= log10(1/threshold) & all_sign_check$max_effect <= log10(threshold))

# Only Negative
sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect < log10(1/threshold))
sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect < log10(1/threshold)) / length(all_sign_check$max_effect)

sign_sum = sign_sum + sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect < log10(1/threshold))

# Only Positive
sum(all_sign_check$min_effect > log10(threshold) & all_sign_check$max_effect > log10(threshold))
sum(all_sign_check$min_effect > log10(threshold) & all_sign_check$max_effect > log10(threshold)) / length(all_sign_check$max_effect)

sign_sum = sign_sum + sum(all_sign_check$min_effect > log10(threshold) & all_sign_check$max_effect > log10(threshold))

# Only Neutral and Negative
sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect <= log10(threshold))
sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect <= log10(threshold)) / length(all_sign_check$max_effect)

sign_sum = sign_sum + sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect <= log10(threshold))

# Only Neutral and Positive
sum(all_sign_check$min_effect >= log10(1/threshold) & all_sign_check$max_effect > log10(threshold))
sum(all_sign_check$min_effect >= log10(1/threshold) & all_sign_check$max_effect > log10(threshold)) / length(all_sign_check$max_effect)

sign_sum = sign_sum + sum(all_sign_check$min_effect >= log10(1/threshold) & all_sign_check$max_effect > log10(threshold))

# Both Negative and Positive
sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect > log10(threshold))
sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect > log10(threshold)) / length(all_sign_check$max_effect)

sign_sum = sign_sum + sum(all_sign_check$min_effect < log10(1/threshold) & all_sign_check$max_effect > log10(threshold))

sign_sum

#### Something isn't working so I will put labels

types = c()
for(each in 1:dim(all_sign_check)[1]){
  if(all_sign_check$min_effect[each] < log10(1/threshold)){
    # Negative branch
    if(all_sign_check$max_effect[each] < log10(1/threshold)) {
      # Full negative
      types = c(types, "Negative")
    } else if(all_sign_check$max_effect[each] <= log10(threshold)) {
      # Negative Neutral
      types = c(types, "Neutral Negative")
    } else {
      # Negative Positive
      types = c(types, "Positive Negative")
    }
  } else {
    # Neutral or Positive
    if(all_sign_check$min_effect[each] > log10(threshold)) {
      # Full positive
      types = c(types, "Positive")
    } else if(all_sign_check$max_effect[each] > log10(threshold)) {
      # Neutral Positive
      types = c(types, "Neutral Positive")
    } else {
      # Neutral
      types = c(types, "Neutral")
    }
  }
}

all_sign_check$type = types

all_sign_check %>%
  ggplot(aes(x = unique_id)) +
  geom_errorbar(aes(ymin = min_effect, ymax = max_effect)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(threshold), lty = 2) +
  geom_hline(yintercept = log10(1/threshold), lty = 2) +
  facet_grid(~ type, scale = "free_x") +
  theme_classic() +
  theme(text = element_text(size=18), axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

all_sign_check$type = factor(all_sign_check$type)

summary(all_sign_check$type)
summary(all_sign_check$type) / length(all_sign_check$type) * 100



########################


# get ranges for each enzyme

d1_range_res = d1 %>%
  mutate(enzyme = factor(enzyme)) %>%
  group_by(unique_id, enzyme) %>%
  summarize(range_effects = abs(max(effects)) + abs(min(effects)), 
            wt_effect = effects[1],
            dev_median = abs(effects[1] - median(effects)),
            dev_mean = abs(effects[1] - mean(effects)),
            pos_mean = mean(effects),
            pos_sd = sd(effects)) 

# get wt residuals for each enzyme

d1_range_res$wt_res = (d1 %>% 
                         mutate(enzyme = factor(enzyme)) %>%
                         group_by(unique_id, enzyme) %>% 
                         mutate(effects = effects + (-1*effects[1])) %>% # normalize to WT
                         summarize(wt_res = sum(abs(effects))) %>% pull(wt_res))


## z-score for each

d1_range_res$zscore = (d1_range_res$wt_effect - d1_range_res$pos_mean) / (d1_range_res$pos_sd)

d1_range_res %>%
  ggplot(aes(x = dev_mean, fill = enzyme)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

sum(d1_range_res$dev_mean > log10(1.5)) / length(d1_range_res$dev_mean)

# Plot correlation between wt effect and mean per position, and add histogram

cor = summary(lm(wt_effect ~ pos_mean, d1_range_res))$r.squared

p = ggplot(d1_range_res, aes(x = pos_mean, y = wt_effect)) +
  geom_smooth(formula = y ~ x, method = 'lm', color = "black") +
  geom_point(size = 3) +
  ggtitle(paste0(round(cor * 100, 1), "% R-squared")) +
  theme_classic() +
  labs(x = "Mean Functional Contribution",
       y = "Functional Effect in WT BG") +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggMarginal(p, type = "histogram", fill = "grey")

#### Again with 1:1 line

p = ggplot(d1_range_res, aes(x = pos_mean, y = wt_effect)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0) +
  geom_abline(intercept = log(1.5), lty = 2) +
  geom_abline(intercept = log(1/1.5), lty = 2) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  #ggtitle(paste0(round(cor * 100, 1), "% R-squared")) +
  theme_classic() +
  labs(x = "Mean Functional Contribution",
       y = "Functional Effect in WT BG") +
  xlim(-5,3) +
  ylim(-5,3) +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggMarginal(p, type = "histogram", fill = "#edae49", color = "#edae49", bins = 75)



dens_plot = d1_range_res %>%
  ggplot(aes(x = dev_mean)) +
  geom_density(fill = "grey", alpha = 0.5) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  labs(x = "log10(Deviation) of Mean Positional Effect from WT",
       y = "Density") +
  theme(text = element_text(size=22),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

dens_plot_d = ggplot_build(dens_plot)$data[[1]]

dens_plot + geom_area(data = subset(dens_plot_d, x>log10(1.5)), aes(x=x,y=y), fill = "#76AA95")

sum(d1_range_res$dev_median > log10(1.5)) / length(d1_range_res$dev_median)
sum(d1_range_res$dev_mean > log10(1.5)) / length(d1_range_res$dev_mean)

## Pie chart for how many are sign

prop = c(length(sign_unique_ids) / length(unique(d1$unique_id)), (1 - length(sign_unique_ids) / length(unique(d1$unique_id))))
labs = c("Crosses Sign", "Doesn't Cross Sign")

pie(c(prop[1], prop[2]), labels = c(prop[1], prop[2]), col = c("#ffe4a6", "#b8eaf5"), family = "sans")

legend("topright", labs, cex = 1.2, fill = c("#ffe4a6", "#b8eaf5"))


############



## Reads all csvs in collated folder and combines them

idio_d = do.call("rbind", lapply(list.files(pattern="idio_df", recursive = T), read_csv))

idio_d %>%
  ggplot(aes(x = mutations, y = idiosync)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  stat_summary(fun=mean, colour="black", geom="crossbar", width=0.2)

## Histogram of idiosynracies

idio_d %>%
  filter(mutations == "Order 1") %>%
  ggplot(aes(x = idiosync)) +
  geom_histogram(binwidth = 0.1, fill = "#edae49") +
  geom_vline(xintercept = log10(1.5), lty = 2, lwd = 1) +
  labs(x = "Idiosyncrasy Metric",
       y = "Genotype Count") +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

idio_first = idio_d %>% 
  filter(mutations == "Order 1") %>%
  pull(idiosync)

sum(idio_first > log10(1.5)) / length(idio_first)
paste(sum(idio_first > log10(1.5)), "/", length(idio_first))

sum(idio_first > log10(2)) / length(idio_first)
sum(idio_first > log10(5)) / length(idio_first)
sum(idio_first > log10(10)) / length(idio_first)

### Faceted idio_d for Order 2, 3, and 4

idio_d %>%
  filter(mutations == "Order 2" | mutations == "Order 3" | mutations == "Order 4") %>%
  ggplot(aes(x = idiosync)) +
  geom_histogram(binwidth = 0.1, fill = "#edae49") +
  geom_vline(xintercept = log10(1.5), lty = 2, lwd = 1) +
  labs(x = "Idiosyncrasy Metric",
       y = "Genotype Count") +
  facet_wrap(~ mutations) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

# What percentage of positions remain poorly characterized at this order

idio_d %>%
  mutate(idiosync = idiosync > log10(1.5)) %>%
  group_by(mutations) %>%
  summarise(idiosync_sig = sum(idiosync),
            idiosync_total = length(idiosync),
            idiosync = sum(idiosync)/length(idiosync) * 100
            )

idio_d %>%
  group_by(mutations) %>%
  summarise(range = max(idiosync) - min(idiosync))


#####
### Looking at MAE of linear model vs feed forward model ###
#####

mae_dif = aic_df %>%
  filter(Step == 1) %>%
  ggplot(aes(x = MAE, y = `Linear MAE`)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  #ggtitle(paste0(round(cor * 100, 1), "% R-squared")) +
  theme_classic() +
  labs(x = "Feed Forward Mean Absolute Error",
       y = "Linear Mean Absolute Error") +
  xlim(0, 4) +
  ylim(0, 4) +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggMarginal(mae_dif, type = "histogram", fill = "#edae49", color = "#edae49", bins = 40)

mae_dif = aic_df %>%
  ggplot(aes(x = MAE, y = `Linear MAE`)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = log10(1.5), lty = 2) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  #ggtitle(paste0(round(cor * 100, 1), "% R-squared")) +
  theme_classic() +
  labs(x = "Feed Forward Mean Absolute Error",
       y = "Linear Mean Absolute Error") +
  xlim(0, 4) +
  ylim(0, 4) +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggMarginal(mae_dif, type = "histogram", fill = "#edae49", color = "#edae49", bins = 40)

##

aic_df %>%
  group_by(Step) %>%
  summarize(mean_mae = mean(MAE),
            sd_mae = sd(MAE),
            mean_mae_lin = mean(`Linear MAE`),
            sd_mae_lin = sd(`Linear MAE`)) %>%
  pivot_longer(cols = c(mean_mae, mean_mae_lin)) %>%
  ggplot(aes(x = Step, y = value)) +
  geom_bar(aes(fill = name), position = "dodge", stat = 'identity') +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  scale_fill_manual(values = c(rgb(160/255, 160/255, 160/255), "#edae49")) +
  theme_classic()


aic_df %>%
  filter(Step <= 4) %>%
  group_by(Step) %>%
  summarize(mean_mae = mean(MAE),
            mean_mae_lin = mean(`Linear MAE`)) %>%
  pivot_longer(cols = c(mean_mae, mean_mae_lin)) %>%
  ggplot(aes(x = Step, y = value)) +
  geom_bar(aes(fill = name), position = "dodge", stat = 'identity') +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  scale_fill_manual(values = c(rgb(160/255, 160/255, 160/255), "#edae49")) +
  theme_classic()

aic_df %>%
  filter(Step <= 4) %>%
  group_by(Step) %>%
  summarize(mean_last_mae = mean(`Last MAE`),
            mean_lin_last_mae = mean(`Linear Last MAE`)) %>%
  pivot_longer(cols = c(mean_last_mae, mean_lin_last_mae)) %>%
  mutate() %>%
  ggplot(aes(x = Step, y = value)) +
  geom_line(aes(color = name), lwd = 1) +
  geom_point(aes(color = name), size = 5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(1.5), lty = 2) +
  scale_color_manual(values = c(rgb(160/255, 160/255, 160/255), "#edae49")) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")


### Idio by order

idio_d %>%
  ggplot(aes(x = idiosync)) +
  geom_histogram(binwidth = 0.1, fill = "#edae49") +
  geom_vline(xintercept = log10(1.5), lty = 2, lwd = 1) +
  labs(x = "Idiosyncrasy Metric",
       y = "Genotype Count") +
  facet_wrap(~ mutations) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))




########################
### Higher order DFs ###
########################

higher_order_list = list.files(pattern="higher_box", recursive = T)

higher_df = data.frame()
for(higher_file in 1:length(higher_order_list)){
  Condition = str_split(str_split(higher_order_list[higher_file], "_")[[1]][3], "/")[[1]][1]
  Measurement = str_split(higher_order_list[higher_file], "_")[[1]][2]
  Enzyme = str_split(higher_order_list[higher_file], "_")[[1]][1]
  appending_file = read_csv(higher_order_list[higher_file])
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


#### Higher order sign assignment

threshold = 1.5

types = c()
for(each in 1:dim(higher_all_sign_check)[1]){
  if(higher_all_sign_check$min_effect[each] < log10(1/threshold)){
    # Negative branch
    if(higher_all_sign_check$max_effect[each] < log10(1/threshold)) {
      # Full negative
      types = c(types, "Negative")
    } else if(higher_all_sign_check$max_effect[each] <= log10(threshold)) {
      # Negative Neutral
      types = c(types, "Neutral Negative")
    } else {
      # Negative Positive
      types = c(types, "Positive Negative")
    }
  } else {
    # Neutral or Positive
    if(higher_all_sign_check$min_effect[each] > log10(threshold)) {
      # Full positive
      types = c(types, "Positive")
    } else if(higher_all_sign_check$max_effect[each] > log10(threshold)) {
      # Neutral Positive
      types = c(types, "Neutral Positive")
    } else {
      # Neutral
      types = c(types, "Neutral")
    }
  }
}

higher_all_sign_check$type = types

higher_all_sign_check %>%
  ggplot(aes(x = unique_id)) +
  geom_errorbar(aes(ymin = min_effect, ymax = max_effect)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = log10(threshold), lty = 2) +
  geom_hline(yintercept = log10(1/threshold), lty = 2) +
  facet_grid(~ type, scale = "free_x") +
  theme_classic() +
  theme(text = element_text(size=18), axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

higher_all_sign_check$type = factor(higher_all_sign_check$type)

order_checker = higher_all_sign_check %>% 
  filter(mutations == "Order 4") %>%
  pull(type) 

summary(order_checker) 
sum(summary(order_checker)) 
summary(order_checker) / length(order_checker) * 100


###

# get ranges for each enzyme

## WHICH ORDER TO TEST?

test_order = "Order 4"

higher_range_res = higher_df %>%
  mutate(mutations = factor(paste("Order", mutations, sep = " "))) %>%
  filter(mutations == test_order) %>%
  mutate(Enzyme = factor(Enzyme)) %>%
  group_by(unique_id, Enzyme) %>%
  summarize(range_effects = abs(max(avg)) + abs(min(avg)), 
            wt_effect = avg[1],
            dev_median = abs(avg[1] - median(avg)),
            dev_mean = abs(avg[1] - mean(avg)),
            pos_mean = mean(avg)) 

# get wt residuals for each enzyme

higher_range_res$wt_res = (higher_df %>% 
                             mutate(mutations = factor(paste("Order", mutations, sep = " "))) %>%
                             filter(mutations == test_order) %>%
                             group_by(unique_id, Enzyme, mutations) %>% 
                             mutate(avg = avg + (-1*avg[1])) %>% # normalize to WT
                             summarize(wt_res = sum(abs(avg))) %>% pull(wt_res))

# Plot correlation between wt effect and mean per position, and add histogram

cor = summary(lm(wt_effect ~ pos_mean, higher_range_res))$r.squared

cor

p = ggplot(higher_range_res, aes(x = pos_mean, y = wt_effect)) +
  geom_point(size = 3) +
  geom_abline(intercept = 0) +
  geom_abline(intercept = log10(1.5), lty = 2) +
  geom_abline(intercept = log10(1/1.5), lty = 2) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  #ggtitle(paste0(round(cor * 100, 1), "% R-squared")) +
  theme_classic() +
  labs(x = "Mean Functional Contribution",
       y = "Functional Effect in WT BG") +
  xlim(-5,5) +
  ylim(-5,5) +
  theme(text = element_text(size=18), axis.text = element_text(size = 16),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))

ggMarginal(p, type = "histogram", fill = "#edae49", color = "#edae49", bins = 75)


sum(abs(higher_range_res$wt_effect - higher_range_res$pos_mean) > log10(1.5)) / length(abs(higher_range_res$wt_effect - higher_range_res$pos_mean) > log10(1.5))


#######
## Do higher order individual points deviate from mean
#######

higher_df_mean = higher_df %>%
  filter(mutations < 5) %>%
  group_by(unique_id) %>%
  summarise(avg = mean(avg))

devs = c()
cur_mutations = c()
for(i in 1:dim(higher_df_mean)[1]) {
  curr_id = higher_df_mean[i,1] %>% pull(unique_id)
  devs = c(devs, abs( (higher_df %>% filter(unique_id == curr_id) %>% pull(avg)) - (higher_df_mean[i,2] %>% pull(avg)) ))
  cur_mutations = c(cur_mutations, (higher_df %>% filter(unique_id == curr_id) %>% pull(mutations)))
}

devs_df = data.frame(devs = devs, muts = cur_mutations)

## WT heterogenity for higher_orders

devs_df %>%
  filter(muts > 1) %>%
  ggplot(aes(x = devs)) +
  geom_histogram(binwidth=0.1, fill = "#edae49") +
  geom_vline(xintercept = log10(1.5), lty = 2, lwd = 1) +
  facet_wrap(~ muts) +
  theme_classic() +
  theme(text = element_text(size=18), axis.text = element_text(size = 16, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())


devs_df %>%
  mutate(muts = factor(muts)) %>%
  group_by(muts) %>%
  summarise(percent = sum(devs > log10(1.5)) / length(devs),
            outlier = sum(devs > log10(1.5)),
            total = length(devs))
