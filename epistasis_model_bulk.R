### Library ###

library(tidyverse)
library(ggrepel) ## For repelling labels on ggplot
library(igraph) ## Network library
library(gtools) ## Permutation library
library(e1071) ## Hamming distance lib
library(svglite) ## SVG file export
library(ggpubr) ## Publication ready ggplot plotting
library(plotly) ## Interactive plotting
library(htmlwidgets)

### Update ###
# If updater TRUE, it will run normally
# If updater FALSE, it will only find files for which folders don't exist

updater = T

### Inputs ###

## legacy positive = '#ff4e45'
## legacy negative = '#45d7ff'

setwd("/Users/karolbuda/OneDrive - UBC/PhD/Epistasis_Lit/__VersionCtrl__/epistasis_analysis_v2")
inputs = list.files("./Input")

# Constants

log_base = 10
transition_threshold = log(1/1.5, log_base) #fold
traj_color = "grey"
color_high = '#ff4e45'
color_low = '#45d7ff'
par(family = "sans")

### Functions ###

data_loading = function() {
  
  d = read.csv(paste("./Input/", file, sep = ""))
  data = d[5:dim(d)[1],1:2]
  
  codes = data.frame()
  simplex_chart = data.frame()
  
  for(j in 1:dim(data)[1]){
    code = c()
    chr = character()
    for(i in 1:length(strsplit(d[1,1], "")[[1]])){
      if(strsplit(data[j,1], "")[[1]][i] == strsplit(d[1,1], "")[[1]][i]){
        code = c(code, -1)
        chr = paste0(chr, strsplit(data[j,1], "")[[1]][i])
      } else {
        code = c(code, 1)
        chr = paste0(chr, strsplit(data[j,1], "")[[1]][i])
      }
    }
    codes = rbind(codes, code)
    simplex_chart = rbind(simplex_chart, c(as.numeric(code), chr))
  }
  
  simplex_chart = unique(simplex_chart) ## Filter out repeats
  
  codes = cbind(codes, data[,2])
  colnames(codes) = c(paste("p", d[3,][!is.na(d[3,])], sep = ""), "effect")
  colnames(simplex_chart) = c(paste("p", d[3,][!is.na(d[3,])], sep = ""), "genotype")
  
  test_codes = gtools::permutations(2, (dim(codes)[2] - 1), c(-1, 1), repeats = TRUE)
  
  for(j in 1:dim(test_codes)[1]) {
    passed = F
    for(i in 1:dim(codes)[1]){
      if(all(test_codes[j,] == codes[-dim(codes)[2]][i,])) {
        passed = T
      }
    }
    if(!passed){
      print(paste0("Data for ", paste(test_codes[j,], collapse = " "), " is missing"))
    }
  }
  codes$effect = as.numeric(codes$effect)
  codes <<- codes
  simplex_chart <<- simplex_chart
}

epistasis_analysis = function() {
  
  order = data.frame()
  
  vars = paste(colnames(codes)[1:length(colnames(codes))-1], collapse="+")
  
  ## No interaction terms model
  current_model = lm(paste("effect ~ (",vars ,")", sep=""), codes)
  order = rbind(order, c("order1", summary(current_model)$r.squared, 
                         summary(current_model)$adj.r.squared, 
                         summary(current_model)$r.squared))
  
  rolling_mae = mean(abs(predict(current_model, codes[-dim(codes)[2]]) - codes$effect))
  
  rolling_last_mae = mean(abs(predict(current_model, codes[which(apply(codes[-dim(codes)[2]] == "1", 1, function(x) all(x))),][-dim(codes)[2]]) - 
                                codes[which(apply(codes[-dim(codes)[2]] == "1", 1, function(x) all(x))),]$effect))
  
  for(k in 2:(dim(codes)[2]-1)){
    ## Needed to paste with **k because lm doesn't like the power term introduced otherwise
    next_model = lm(paste("effect ~ (",vars ,")**", k, sep=""), codes)

    if(!is.na(anova(current_model, next_model)[6][2,])) {
      print(paste0("Next model p-value: ", anova(current_model, next_model)[6][2,]))
      current_model = next_model
      # Add order into order list
      order = rbind(order, c(paste("order", k, sep=""), summary(current_model)$r.squared, 
                             summary(current_model)$adj.r.squared, 
                             summary(current_model)$r.squared - as.numeric(order[k-1,2])))
      
      # Add MAE to the list
      rolling_mae = c(rolling_mae, mean(abs(predict(current_model, codes[-dim(codes)[2]]) - codes$effect)))
      
      rolling_last_mae = c(rolling_last_mae, mean(abs(predict(current_model, codes[which(apply(codes[-dim(codes)[2]] == "1", 1, function(x) all(x))),][-dim(codes)[2]]) - 
                                                        codes[which(apply(codes[-dim(codes)[2]] == "1", 1, function(x) all(x))),]$effect)))
      
    } else{
      print(paste0("Highest Model Order: ", k-1))
      print(paste0("Next model p-value: ", anova(current_model, next_model)[6][2,]))
      break
    }
    print(paste0("Setting next model as current model...: "))
    print(paste0("Analyzing Next Model Order: ", k))
  }
  
  ## Short term MAE global export
  
  rolling_mae <<- rolling_mae
  rolling_last_mae <<- rolling_last_mae
  
  ## Prepare output
  
  pos_out = as.data.frame(summary(current_model)$coef[,1]*2)
  pos_out[1,1] = summary(current_model)$coef[,1][1] ## Don't multiply intercept by 2
  
  pos_out$names = rownames(pos_out)
  rownames(pos_out) = c()
  colnames(pos_out) = c("effect", "indices")
  
  colnames(order) = c("model_order", "R2", "(Adjusted)", "Delta R2")
  
  ## Convert : and p characters
  pos_out$indices = gsub(':', '|', pos_out$indices)
  pos_out$indices = gsub('p', '', pos_out$indices)
  pos_out$indices[1] = "INTERCEPT"
  
  ## Genotype predicted out
  
  vars_sep = unlist(strsplit(vars, "+", fixed=TRUE))
  simplex_chart = cbind(apply(simplex_chart[-dim(simplex_chart)[2]], 2, as.numeric), simplex_chart[dim(simplex_chart)[2]])
  
  preds = c()
  for(i in 1:dim(simplex_chart)[1]) {
    preds = c(preds, predict(current_model, simplex_chart[-dim(simplex_chart)[2]][i,]))
  }
  
  pred_df = cbind(simplex_chart[dim(simplex_chart)[2]], preds)
  colnames(pred_df) = c("genotype", "predicted effect")
  
  ## Truncated First order model predicted out vs real data
  
  simplex_chart = simplex_chart %>%
    arrange_at(names(simplex_chart)[-dim(simplex_chart)[2]])
  
  ## Writing csvs
  
  write.csv(pos_out[,c(2,1)],"pos_out.csv", row.names = FALSE)
  write.csv(order,"model_order.csv", row.names = FALSE)
  write.csv(pred_df, "gen_out.csv", row.names = FALSE)
  #write.csv(pred_compare_df, "pred_out.csv", row.names = FALSE)
  
  ## Provide feedback in log.txt form if there are singularities because of missing data AND the p-value of next model if exists
  
  if(length(which(is.na(current_model$coef))) > 0) {
    write(c("-------------", paste("Variable", names(current_model$coef[which(is.na(current_model$coef))]), "= NA")), "log.txt")
  }
  if(!is.na(anova(current_model, next_model)[6][2,])) {
    write(c("-------------", paste0("Highest Model Order: ", k-1), paste0("Next model p-value: ", anova(current_model, next_model)[6][2,])), "log.txt", append = T)
  }
  
  ## Plotting
  
  pos_out$order = c(NA, paste0("Order ", str_count(pos_out$indices[-1], "[|]") + 1, " (", formatC(as.numeric(order[str_count(pos_out$indices[-1], "[|]") + 1, 2]), digits = 3), ")"))
  
  pos_out$indices = factor(pos_out$indices, levels = pos_out$indices)
  
  stat_plot = ggplot(pos_out[-1,], aes(x = indices, y = effect)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat="identity", fill = traj_color, color = "black") +
    #geom_text(aes(label=formatC(effect, digits=2), y = effect + 0.1*sign(effect)*abs(max(effect))), size = 3) +
    facet_wrap(~ order, scales = "free_x") + 
    labs(x = "Mutation(s)", y = "Fold Effect on Activity") +
    theme_classic() +
    theme(text = element_text(size=22), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  
  ggsave("epistatic_effect_no_lables.svg", plot = stat_plot, width = 12, height = 12)
  
  ## Combining points into averages 
  
  newer_codes = codes %>%
    group_by_at(names(codes)[-length(names(codes))]) %>%
    summarise(avg = log(signif(mean(log_base^effect), digits = 3), log_base), 
              cv = (sd(log_base^effect) / mean(log_base^effect))^2)
  
  ## This just stores newer_codes as a global variable
  
  newer_codes <<- newer_codes
  # Code above ensures that I don't take mean of the logs, first I transform them back with 10^x
  # I also use the signif function because too many digits makes log10(1) != 0 for some reason
  
  replace_occurence = function(x, replaced, value){
    x[x == replaced] = value
    return(x)
  }
  # Used to remove values in data set which don't have either a WT-like or mutant-like genotype... i.e. we can't show
  # x0001 if we don't have both 00001 and 10001 in our data, so we remove the incorrectly calculated x0001
  
  codes_effects = c()
  codes_identity = c()
  codes_position = c()
  
  subtract = function(x) {
    if(length(x) == 1) {
      return(NA)
    } else {
      return(x[2] - x[1])
    }
  } ## Very primitive function with little checks, however after taking the mean there should only be 1 or 2 values
  ## If there is 1 value we can't display it because there is no interaction to compare so we return an NA and remove it
  
  for(i in 1:(length(names(newer_codes)) - 1)) {
    
    current = newer_codes %>%
      group_by_at(names(newer_codes)[(1:(length(codes)-1))[(1:(length(codes)-1)) != i]]) %>%
      summarise(idk = subtract(avg)) %>%
      na.omit()
    
    if(dim(current)[1] > 0) {
      codes_effects = c(codes_effects, current %>% pull(idk))
      
      for(j in 1:(dim(current)[1])) {
        str_code = paste(as.character(replace_occurence(current[j, 1:(length(names(current)) - 1)], -1, 0)), collapse = "")
        codes_identity = c(codes_identity, paste(c(substr(str_code, 0, i-1), substr(str_code, i, (length(names(newer_codes))-2))), collapse="x"))
      }
      
      codes_position = c(codes_position, rep(names(newer_codes)[i], dim(current)[1]))
    }
    
  }
  
  new_codes = data.frame(positions = codes_position, effects = codes_effects, identity = codes_identity) 
  
  new_codes$mut = factor(str_count(new_codes$identity, "1") + 1)
  
  new_codes$positions = factor(new_codes$positions, levels = unique(new_codes$positions))
  
  new_codes %>%
    mutate(condition = paste(cond, str_remove(file, ".csv"), sep = "_")) %>%
    write.csv(file = paste(cond, str_remove(file, ".csv"), "box_plot.csv", sep = "_"), row.names = F)
  
  ## Genotype predicted out
  
  effects_vector = newer_codes %>% pull(avg)
  
  vars_sep = unlist(strsplit(vars, "+", fixed=TRUE))
  simplex_chart = cbind(apply(simplex_chart[-dim(simplex_chart)[2]], 2, as.numeric), simplex_chart[dim(simplex_chart)[2]])
  
  ## Arrange simplex chart the same way for all
  
  simplex_chart = simplex_chart %>%
    arrange_at(names(simplex_chart)[-dim(simplex_chart)[2]])
  
  preds = c()
  for(i in 1:dim(simplex_chart)[1]) {
    preds = c(preds, predict(current_model, simplex_chart[-dim(simplex_chart)[2]][i,]))
  }
  
  pred_df = cbind(simplex_chart, effects_vector, preds)
  
  ## Truncated First order model predicted out vs real data
  
  preds = c()
  for(i in 1:dim(simplex_chart)[1]) {
    preds = c(preds, predict(lm(paste("effect ~ (",vars ,")", sep=""), codes), simplex_chart[-dim(simplex_chart)[2]][i,]))
  }
  
  pred_df = cbind(pred_df, preds)
  
  ### WT vs average effect
  
  wt_effects = new_codes$effect[which(new_codes$identity %in% new_codes$identity[-grep("1", new_codes$identity)])]
  
  ## Adds NA to the places where variable is missing
  
  missing = which(!(vars_sep %in% as.character(new_codes[which(new_codes$identity %in% new_codes$identity[-grep("1", new_codes$identity)]),]$positions)))
  
  if(length(missing) > 0) {
    for(i in 1:length(missing)) {
      wt_effects = append(wt_effects, 0, after=(missing[i]-1))
    }
  }
  
  wt_vs_avg = data.frame(indices = pos_out[2:(1+length(wt_effects)),]$indices, average_effect = pos_out[2:(1+length(wt_effects)),]$effect, wt_effect = wt_effects)
  
  wt_vs_avg_plot = data.frame(index = rep(wt_vs_avg$indices, 2), effects = c(wt_vs_avg$average_effect, wt_vs_avg$wt_effect), 
                              type = c(rep("average", dim(wt_vs_avg)[1]), rep("wt", dim(wt_vs_avg)[1])))
  
  wt_codes = codes[-length(codes)]
  wt_codes[wt_codes == -1] = 0
  wt_codes = unique(wt_codes)
  
  wt_codes = wt_codes %>%
    arrange_at(names(wt_codes))
  
  wt_model = lm(paste("effect ~ (",vars ,")", sep=""), codes)
  
  wt_model$coefficients = c(0, wt_effects)
  names(wt_model$coefficients) = names(lm(paste("effect ~ (",vars ,")", sep=""), codes)$coefficients)
  
  preds = c()
  for(i in 1:dim(wt_codes)[1]) {
    preds = c(preds, predict(wt_model, wt_codes[i,]))
  }
  
  pred_df = cbind(pred_df, preds)
  
  colnames(pred_df) = c(colnames(simplex_chart), "observed effect", "epistatic effect", "average effect", "WT effect")
  
  if(length(missing) > 0) {
    for(i in 1:length(missing)) {
      pred_df$`WT effect`[which(pred_df[,missing[i]] == 1)] = NA
    }
  }
  
  ## Add in cvs to pred_df
  
  pred_df$cv = newer_codes$cv
  
  mutations = c()
  for(i in 1:dim(pred_df)[1]){
    mutations = c(mutations, length(which(pred_df[i,1:dim(codes)[2]-1] == 1)))
  }
  
  pred_df$mutations = as.factor(mutations)
  
  pred_df = pred_df %>% add_column(numbers = pred_df %>% 
                                     select(unlist(strsplit(vars, "+", fixed=TRUE))) %>% 
                                     unite("numbers", 1:dim(.)[2], sep = "") %>% 
                                     mutate(numbers = str_replace_all(numbers, "-1", "0")) %>% 
                                     pull(numbers), .before = "genotype")
  
  rss_pred_df = pred_df %>%
    mutate(mutations = as.numeric(mutations) - 1)
  
  ## Mean Absolute Error (MAE) calculation
  
  ff_prediction_df = data.frame(mae = c(), mae_last = c(), step = c())
  
  last_one = rss_pred_df[which(apply(rss_pred_df[,1:dim(codes[-dim(codes)[2]])[2]] == "1", 1, function(x) all(x))),]
  
  mean_resi = mean(as.numeric(unlist(abs(rss_pred_df$`observed` - rss_pred_df$`WT effect`))))
  mean_last_resi = mean(as.numeric(unlist(abs(last_one$`observed` - last_one$`WT effect`))))
  
  ff_prediction_df = rbind(ff_prediction_df, c(mean_resi, mean_last_resi, 1))
  
  ##################################################
  ## Quick ratio calculation for pymol csv output ##
  ###################################################
  
  pred_df <<- pred_df
  new_codes <<- new_codes
  wt_effects <<- wt_effects
  ff_prediction_df <<- ff_prediction_df
}

blanket_analysis = function() {
  all_stats = new_codes %>%
    group_by(positions) %>%
    summarise(mean = mean(effects),
              median = median(effects),
              top = max(effects), 
              bottom = min(effects),
              range = max(effects) - min(effects),
              skew = (3*(mean(effects) - median(effects)))/sd(effects),
              pos_pct = sum(effects > 0) / length(effects),
              neg_pct = sum(effects < 0) / length(effects))
  
  wt_res = c()
  
  for(i in 1:length(unique(new_codes$positions))) {
    wt_res = c(wt_res, sum(((new_codes %>%
                               filter(positions == unique(new_codes$positions)[i]) %>%
                               filter(mut == 1) %>%
                               pull(effects)) - (new_codes %>% filter(positions == unique(new_codes$positions)[i]) %>% pull(effects)))^2))
  }
  
  all_stats = bind_cols(all_stats, wt_res = wt_res)
  
  ## Pymol ratio export
  
  pymol_csv = pred_df[1:(length(simplex_chart)-1)]
  
  ind = c()
  for(i in 1:dim(pymol_csv)[1]) {
    if(length(colnames(pymol_csv[i,][which(pymol_csv[i, ] == 1)])) > 0) {
      ind = c(ind, paste(colnames(pymol_csv[i,][which(pymol_csv[i, ] == 1)]), collapse = "|"))
    } else {
      ind = c("INTERCEPT")
    }
  }
  
  ind = gsub('p', '', ind)
  
  ### CURRENTLY NOT COMPATIBLE WITH PIE CHART
  
  pymol_csv = data.frame(indices = ind, 
                         magnitude = log_base^(abs(pred_df$`observed effect` - pred_df$`WT effect`)),
                         sign = pred_df$`observed effect`*pred_df$`WT effect`,
                         negative = pred_df$`observed effect` < pred_df$`WT effect`) %>%
    mutate(effect = log(magnitude, log_base)*c(1,-1)[as.numeric(negative) + 1],
           magnitude = magnitude > 1.5,
           sign = as.numeric(sign < 0)) %>%
    filter(magnitude == T | sign == T) %>%
    select(-c(magnitude, negative))
  
  write.csv(pymol_csv, "pymol_csv.csv", row.names = F)
  
  ## WT ratio plots
  
  test = c()
  
  for(i in 1:length(unique(new_codes$positions))) {
    if(!dim(new_codes %>%
            filter(positions == unique(new_codes$positions)[i]) %>%
            filter(mut == 1))[1] > 0) {
      test = c(test, rep(NA, dim(new_codes %>% filter(positions == unique(new_codes$positions)[i]))[1]))
    } else {
      test = c(test, ((new_codes %>% filter(positions == unique(new_codes$positions)[i]) %>% pull(effects)) - new_codes %>%
                        filter(positions == unique(new_codes$positions)[i]) %>%
                        filter(mut == 1) %>%
                        pull(effects)))
    }
  }
  
  test <<- test
}

network_analysis = function(eff = "observed") {
  pred_df = pred_df[order(pred_df$mutations),]
  
  mut = length(strsplit(pred_df$genotype[1], "")[[1]]) ## Get length of WT as indicator of mutations
  
  ids = pred_df[,1:mut]
  ids[ids == -1] = 0
  ids = apply(ids, 1, function(x) paste(as.character(x), collapse = ""))
  
  if(eff == "WT") {
    codes = data.frame(id = ids, effect = pred_df$`WT effect`)
  } else if(eff == "average") {
    codes = data.frame(id = ids, effect = pred_df$`average effect`)
  } else if(eff == "epistatic") {
    codes = data.frame(id = ids, effect = pred_df$`epistatic effect`)
  } else if(eff == "observed") {
    codes = data.frame(id = ids, effect = pred_df$`observed effect`)
  } else if(eff == "ratio") {
    codes = data.frame(id = ids, effect = log_base^(pred_df$`observed effect`)/log_base^(pred_df$`WT effect`))
  }
  
  muts = c()
  for(i in 1:length(codes$id)) {
    muts = c(muts, sum(as.numeric(str_split(codes$id[i], "")[[1]])))
  }
  
  codes = codes[order(muts),] # order by mutation level
  codes$muts = muts[order(muts)]
  
  #Adds in the cvs to be there
  
  codes$cv = pred_df$cv
  
  ##
  
  codes_perm = data.frame()
  
  for(i in 1:length(codes$id)) {
    codes_perm = rbind(codes_perm, as.numeric(str_split(codes$id[i], "")[[1]]))
  }
  
  codes_perm = as.matrix(codes_perm)
  
  from = c()
  to = c()
  
  for(i in 1:dim(codes_perm)[1]) {
    for(j in 1:dim(codes_perm)[1]) {
      if(hamming.distance(codes_perm[i, ], codes_perm[j, ]) == 1 & sum(codes_perm[i, ]) < sum(codes_perm[j, ])) {
        ## Hamming distance ensures only 1 change and sum ensures no back tracking
        from = c(from, codes$id[i])
        to = c(to, codes$id[j])
      }
    }
  }
  
  links = data.frame(from = from, to = to)
  
  ## change effects with data.frame that contains IDs and effect
  
  pascalTriangle <- function(h) {
    lapply(0:h, function(i) choose(i, 0:i))
  }
  
  # Define max vector
  
  weights = c()
  magnitude = c()
  
  pascal_terms = tail(pascalTriangle(mut),1)[[1]][-length(tail(pascalTriangle(mut),1)[[1]])]
  
  ## Function which decides whether link is possible or not based on transition threshold
  
  for(i in 1:dim(links)[1]) {
    
    dif = codes$effect[which(codes$id %in% links[i,2])] - codes$effect[which(codes$id %in% links[i,1])]
    magnitude = c(magnitude, dif)
    
    if(is.na(dif)) {
      weights = c(weights, 1)
    } else if(dif >= transition_threshold & eff != "ratio") {
      weights = c(weights, 2)
    } else {
      weights = c(weights, 1)
    }
  }
  
  links = cbind(links, weights, magnitude)
  
  ## Currently spacing is based on pascal terms instead of existing terms... may need to change that for missing data
  
  ## x coordinates of nodes
  
  test = c(0)
  
  for(i in unique(codes$muts)[-1]) {
    test = c(test, rep(i, length(codes$muts[codes$muts == i])))
  }
  
  
  ## y coordinates of nodes
  
  pasc_terms = c()
  for(i in unique(codes$muts)) {
    pasc_terms = c(pasc_terms, length(codes$muts[codes$muts == i]))
  }
  
  pasc_terms = pasc_terms[-1]
  
  test_2 = c(max(pasc_terms)/2)
  
  for(x in 1:length(pasc_terms)) {
    test_2 = c(test_2, max(pasc_terms)-(max(pasc_terms)/(pasc_terms[x]+1))*(1:(pasc_terms[x])))
  }
  
  l <- matrix(c(test, test_2), nrow = (tail(cumsum(pasc_terms), 1) + 1), ncol = 2)
  
  net_test = graph.data.frame(links, codes, directed = T)
  
  ## Calculate shortest paths
  
  if(eff == "observed") {
    likely_path = shortest_paths(net_test, from = codes$id[1], to = codes$id[which.max(codes$effect)], 
                                 weights = 1/(E(net_test)$magnitude + max(c(abs(min(E(net_test)$magnitude)), abs(max(E(net_test)$magnitude)))) + 1))
    
    derived_likely_path = shortest_paths(net_test, from = codes$id[1], to = codes$id[length(codes$id)], 
                                         weights = 1/(E(net_test)$magnitude + max(c(abs(min(E(net_test)$magnitude)), abs(max(E(net_test)$magnitude)))) + 1))
    
    if(codes$id[which.max(codes$effect)] != codes$id[1]) {
      accessible_path = c(codes$id[1])
      current_code = codes$id[1]
      previous_code = codes$id[1]
      
      while(current_code != codes$id[which.max(codes$effect)]) {
        if(any(links$magnitude[which(links$from == previous_code)]) > 0) {
          current_code = links$to[which(links$from == previous_code)[which.max(links$magnitude[which(links$from == previous_code)])]]
          
          accessible_path = c(accessible_path, current_code)
          
          previous_code = current_code
        } else {
          break
        }
      }
      
      accessible_path <<- accessible_path
    } 
    
    likely_path <<- likely_path$vpath[[1]]$name
    derived_likely_path <<- derived_likely_path$vpath[[1]]$name
    
  }
  
  
  if(eff != "ratio") {
    
    net_test_pos = graph.data.frame(links[links$weights > 1,], codes, directed = T)
    
    paths_pos = all_simple_paths(net_test_pos, from = codes$id[1], to = codes$id[which.max(codes$effect)])
    
    # Vector in order of edges that don't lead to dead ends
    
    if(length(paths_pos) > 0) {
      pos = rep(1, dim(links)[1])
      for(i in 1:length(paths_pos)) {
        paths_pos_cur = as_ids(paths_pos[[i]])
        for(j in 1:dim(links)[1]){
          for(i in 1:(length(paths_pos_cur) - 1)) {
            if(links[j,1] == paths_pos_cur[i] & links[j,2] == paths_pos_cur[i+1]) {
              pos = replace(pos, j, 2)
            }
          }
        }
      }
    }
    
  }
  ## Creates gradient of colors
  
  rbPal = colorRampPalette(c(color_low,'white', color_high))
  
  ## This sets 0 in the center with max and min being equal to pos and neg of the largest value
  idk = max(c(abs(max(codes$effect, na.rm = T)), abs(min(codes$effect, na.rm = T))))
  col = rbPal(length(seq(-idk, idk, idk/50)))[as.numeric(cut(seq(-idk, idk, idk/50), breaks = length(seq(-idk, idk, idk/50))))]
  
  cols = c()
  #val = c()
  for(i in 1:length(codes$effect)) {
    if(!is.na(codes$effect[i])){
      cols = c(cols, col[which.min(abs(seq(-idk, idk, idk/50) - codes$effect[i]))])
      #val = c(val, codes$effect[which.min(abs(seq(-idk, idk, idk/50) - codes$effect[i]))])
    } else {
      cols = c(cols, 'black')
    }
  }
  
  codes$colors = cols
  
  ## Colors for ratio plot
  
  if(eff == "ratio") {
    
    ## Creates feed forward ratio plot
    
    ratio_codes = codes
    
    for(i in 1:dim(filter(ratio_codes, muts > 2))[1]) {
      
      node = filter(ratio_codes, muts > 2)$id[i]
      curr_mut = filter(ratio_codes, muts > 2)$muts[i]
      
      prev = c()
      nodes = c()
      new_prev = c()
      for(j in 1:(curr_mut-2)) {
        if(j > 1) {
          for(k in 1:length(nodes)) {
            prev = c(prev, filter(links, to == nodes[k]) %>% pull(from))
            new_prev = c(new_prev, filter(links, to == nodes[k]) %>% pull(from))
          }
          nodes = new_prev
          new_prev = c()
        } else {
          prev = c(prev, filter(links, to == node) %>% pull(from))
          nodes = prev
        }
      }
      
      curr_prod = prod(ratio_codes[which(ratio_codes$id %in% prev), ]$effect)
      curr_error = sum(ratio_codes[which(ratio_codes$id %in% prev), ]$cv)
      
      ratio_codes[ratio_codes$id == node,]$effect = ratio_codes[ratio_codes$id == node,]$effect / curr_prod
      ratio_codes[ratio_codes$id == node,]$cv = ratio_codes[ratio_codes$id == node,]$cv + curr_error
    }
    
    # Adds idiosinracy value
    
    ratio_codes <<- ratio_codes
    
    ratio_effects = codes$effect
    
    no_imp = which(is.na(pred_df$`WT effect`))[which(which(is.na(pred_df$`WT effect`)) %in% which(is.na(pred_df$`observed effect`)))]
    
    ratio_effects[no_imp] = 1e-4 # Place holder value to definitely make these negative and positive in error checking
    ratio_effects[is.na(ratio_effects)] = 1e4
    
    cols = c()
    
    for(i in 1:length(ratio_effects)) {
      if(ratio_effects[i] > 1) {
        cols = c(cols, color_high)
      } else if(ratio_effects[i] < 1) {
        cols = c(cols, color_low)
      } else {
        cols = c(cols, 'white')
      }
    }
    
    ratio_codes_effects = ratio_codes$effect
    
    ratio_cols = c()
    
    for(i in 1:length(ratio_codes_effects)) {
      if(is.na(ratio_codes_effects[i])) {
        ratio_cols = c(ratio_cols, 'black')
      } else if(ratio_codes_effects[i] > 1) {
        ratio_cols = c(ratio_cols, color_high)
      } else if(ratio_codes_effects[i] < 1) {
        ratio_cols = c(ratio_cols, color_low)
      } else {
        ratio_cols = c(ratio_cols, 'white')
      }
    }
    
    ratio_codes$colors = ratio_cols  
  }
  
  V(net_test)$frame.color <- "black"
  if(eff != "ratio") {
    V(net_test)$color <- codes$colors
  } else {
    V(net_test)$color <- cols
  }
  V(net_test)$size <- 15
  E(net_test)$lty <- c(2, 1)[links$weights]
  E(net_test)$width <- (E(net_test)$weights / 2)
  E(net_test)$arrow.mode <- 0
  
  E(net_test)$color <- "#C1C1C1"
  
  svglite(file = paste("traj_", eff, ".svg", sep=""), width = 12, height = 12)
  
  if(length(names(simplex_chart)) - 1 > 5) {
    if(eff == "ratio") {
      
      ratio_effects = as.character(signif(codes$effect, 2))
      
      ratio_effects[no_imp] = '-'
      ratio_effects[is.na(ratio_effects)] = '+'
      
      plot(net_test, layout = l, vertex.label = c(ratio_effects, codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(0.7, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])),
           vertex.size = 7, family = "sans")
      
      dev.off()
      
      svglite(file = paste("traj_", eff, "_feedforward.svg", sep=""), width = 12, height = 12)
      
      V(net_test)$color <- ratio_codes$colors
      
      plot(net_test, layout = l, vertex.label = c(signif(ratio_codes_effects, 2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(0.7, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])),
           vertex.size = 7, family = "sans")
      
      dev.off()
      
    } else {
      plot(net_test, layout = l, vertex.label = c(signif(log_base^codes$effect,2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(0.7, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])),
           vertex.size = 7, family = "sans")
      
      dev.off()
      
      tiff(file = paste("traj_", eff, ".tiff", sep=""), units = "in", width = 7, height = 7, res = 300)
      
      plot(net_test, layout = l, vertex.label = c(signif(log_base^codes$effect,2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(0.7, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])),
           vertex.size = 7, family = "sans", vertex.label.cex=0.5)
      
      dev.off()
      
    }
  } else {
    if(eff == "ratio") {
      
      ratio_effects = as.character(signif(codes$effect, 2))
      
      ratio_effects[no_imp] = '-'
      ratio_effects[is.na(ratio_effects)] = '+'
      
      plot(net_test, layout = l, vertex.label = c(ratio_effects, codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(1.3, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])), family = "sans")
      
      dev.off()
      
      svglite(file = paste("traj_", eff, "_feedforward.svg", sep=""), width = 12, height = 12)
      
      V(net_test)$color <- ratio_codes$colors
      
      plot(net_test, layout = l, vertex.label = c(signif(ratio_codes_effects, 2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(1.3, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])), family = "sans")
      
      dev.off()
      
    } else {
      plot(net_test, layout = l, vertex.label = c(signif(log_base^codes$effect,2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(1.3, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])), family = "sans")
      
      dev.off()
      
      tiff(file = paste("traj_", eff, ".tiff", sep=""), units = "in", width = 7, height = 7, res = 300)
      
      plot(net_test, layout = l, vertex.label = c(signif(log_base^codes$effect,2), codes$id), 
           vertex.label.dist=c(rep(0, dim(codes)[1]), rep(1.3, dim(codes)[1])), 
           vertex.label.degree=pi/2, vertex.label.color="black", 
           vertex.label.font=c(rep(2, dim(codes)[1]), rep(1,dim(codes)[1])), family = "sans",
           vertex.label.cex=0.5)
      
      dev.off()
      
    }
  }
  
  
  if(eff == "observed") {
    net_codes <<- codes
    obs_links <<- links
  } else if(eff == "ratio") {
    ratio_export = codes
    ratio_export <<- ratio_export
  }
  
  #legend_image = as.raster(rbPal(100), ncol=1)
  #plot(c(0,2),c(0,1),type = 'n', axes = F, xlab = '', ylab = '')
  #text(x=1.3, y = seq(0,1,l=5), labels = seq(0,1,l=5))
  #rasterImage(legend_image, 0, 0, 1,1)
  net_test <<- net_test
  net_test_pos <<- net_test_pos
  links <<- links
  
}

feed_forward = function() {
  
  #####################################
  ## FEED FORWARD COEFFICIENTS PLOT ###
  #####################################
  
  ratio_codes_ff = ratio_codes
  ratio_codes_ff[ratio_codes_ff$muts == 1,]$effect = log_base^(pred_df %>% filter(mutations == 1) %>% pull(`observed effect`))
  
  temp_plot_2 = ratio_codes_ff %>%
    filter(muts > 0) %>%
    mutate(muts = paste("Step ", muts, sep = ""), 
           effect = log(effect, log_base), 
           id = factor(id, rev(unique(id)))) %>%
    ggplot(aes(x = id, y = effect)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat = 'identity', color = "black") +
    facet_wrap(~ muts, scale = "free_x") +
    xlab("Genotype") +
    ylab("log10(Epistatic Coefficient)") +
    theme_classic() +
    theme(text = element_text(size=22), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  
  ggsave("feed_forward_coefficients.svg", plot = temp_plot_2, width = 12, height = 12)
  
  ##
  
  write.csv(ratio_codes, "feedforward_values.csv", row.names = F)
  write.csv(ratio_export, "ratio_export.csv", row.names = F)
  
  new_codes %>%
    select(-c(effects)) %>%
    bind_cols(effects = test, enzyme = enzyme_name, type = measure_type, cond = cond) %>%
    write_csv(., paste(c(str_remove(file, ".csv"), "_wt_rel.csv"), collapse = ""), col_names = F)
  
  new_codes %>%
    select(-c(effects)) %>%
    bind_cols(effects = new_codes$effects, enzyme = enzyme_name, type = measure_type, cond = cond) %>%
    write_csv(., paste(c(str_remove(file, ".csv"), "_2d_box.csv"), collapse = ""), col_names = F)
  
  # Truncated model calculation
  
  vars = paste(colnames(codes)[1:length(colnames(codes))-1], collapse="+")
  
  ff_codes = newer_codes %>% select(-c(avg, cv))
  ff_codes[ff_codes == -1] = 0
  ff_codes$effect = newer_codes$avg
  
  # I can just construct a linear model with coefficients
  
  coefs = c()
  
  for(k in 2:(dim(ff_codes)[2] - 1)) {
    
    all = c(0, wt_effects)
    
    trunc_m = lm(paste("effect ~ (",vars ,")**", k, sep=""), ff_codes)
    coefs = c(coefs, ratio_codes %>% filter(muts == k) %>% map_df(rev) %>% mutate(effect = log(effect, log_base)) %>% pull(effect))
    all = c(all, coefs)
    
    trunc_m$coefficients = all
    names(trunc_m$coefficients) = names(lm(paste("effect ~ (",vars ,")**", k, sep=""), ff_codes)$coefficients)
    
    preds = c()
    for(i in 1:dim(ff_codes)[1]) {
      preds = c(preds, predict(trunc_m, ff_codes[i,]))
    }
    
    pred_df = cbind(pred_df, preds)
    if(k == (dim(ff_codes)[2] - 1)) {
      names(pred_df)[dim(pred_df)[2]] = "final step"
    } else {
      names(pred_df)[dim(pred_df)[2]] = paste(k, "step", sep = " ")
    }
    
    ## Get residuals from 1:1 line vs null slope line just to the points I'm not fitting directly:
    # Maybe I should do average predictability instead?
    
    rss_pred_df = pred_df %>%
      mutate(mutations = as.numeric(mutations) - 1)
    
    last_one = rss_pred_df[which(apply(rss_pred_df[,1:dim(codes[-dim(codes)[2]])[2]] == "1", 1, function(x) all(x))),]
    
    mean_resi = mean(pull(abs(rss_pred_df$`observed` - rss_pred_df[dim(pred_df)[2]])))
    mean_last_resi = mean(pull(abs(last_one$`observed` - last_one[dim(pred_df)[2]])))
    sd_resi = sd(pull(abs(rss_pred_df$`observed` - rss_pred_df[dim(pred_df)[2]])))
    
    ff_prediction_df = rbind(ff_prediction_df, c(mean_resi, mean_last_resi, k))
    
  }
  
  pred_df_ff <<- pred_df
  ff_prediction_df <<- ff_prediction_df
  
  write_csv(pred_df, paste("pred_df_", cond, "_", str_replace(file, ".csv", ""), ".csv", sep = ""))
  
}

higher_order_box = function() {
  
  
  subtract = function(x) {
    if(length(x) == 1) {
      return(NA)
    } else {
      return(x[2] - x[1])
    }
  } 
  
  
  `%notin%` <- Negate(`%in%`) # creates a quick "not in" piper
  
  position_names = names(codes)[1:(length(names(codes))-1)]
  final_work = tibble()
  int_work = tibble()
  next_work = tibble()
  
  ###############
  # Iteration 1 #
  ###############
  
  # Some for loop to change combos to be 1 shorted, and "this" to be the combinations
  
  places = dim(codes)[2] - 1
  
  for(counter in 1:places) {
    combos = 1:(places - (counter - 1))
    if(counter == 1) {
      
      for(combo in 1:length(combos)) {
        curr_work = newer_codes[unique(c(which(apply(newer_codes[, combos[combo]] == -1, 1, all)), which(apply(newer_codes[, combos[combo]] == 1, 1, all)))),]
        curr_work = curr_work %>% 
          group_by_at(which(1:places %notin% combos[combo])) %>%
          # Remove chosen positions to allow for summarizing with the subtract function 
          select(-c(combos[combo])) %>%
          summarise(avg = subtract(avg))
        # Artifically add in columns with Xs to account for where the positions were removed
        curr_work = add_column(curr_work, place_holder = rep("x", dim(curr_work)[1]), .before = combos[combo], .name_repair = "minimal")
        names(curr_work)[which(names(curr_work) == "place_holder")] = combo
        
        curr_work = curr_work %>%
          unite("genotype", 1:places, sep = "") %>%
          mutate(genotype = str_replace_all(genotype, "-1", "0"))
        
        int_work = bind_rows(int_work, curr_work)
        
      }
      
    } else {
      this = gtools::combinations(places, counter - 1, 1:places, repeats = F)
      
      # This is 2nd order
      
      for(i in 1:dim(this)[1]) {
        for(j in 1:dim(this)[2]) {
          if(j == 1) {
            dummy_work = int_work[which(str_sub(int_work$genotype, this[i, j], this[i, j]) == "x"),]
          } else {
            dummy_work = dummy_work[which(str_sub(dummy_work$genotype, this[i, j], this[i, j]) == "x"),]
          }
        }
        dummy_work = dummy_work %>% separate(genotype, letters[1:(places+1)], "") %>% select(-c(1))
        colnames(dummy_work) = c(position_names, "avg")
        
        for(del in 1:length(this[i,])) {
          that = this[i,][order(this[i,], decreasing = T)]
          dummy_work = dummy_work %>% select(-c(that[del]))
        }
        
        dummy_work[cbind(dummy_work[1:length(dummy_work)-1] == "0", avg = FALSE)] <- "-1"
        dummy_work = dummy_work %>% summarise(across(1:(places - counter + 2), as.numeric))
        
        for(combo in 1:length(combos)) {
          curr_work = dummy_work[unique(c(which(apply(dummy_work[, combos[combo]] == -1, 1, all)), which(apply(dummy_work[, combos[combo]] == 1, 1, all)))),]
          curr_work = curr_work %>% 
            #ungroup()
            # Arrange by chosen positionsWW
            #arrange(combos[i,][1], combos[i,][2]) %>%
            # Group by the rest of position so we can summarize them 
            group_by_at(which(1:(places - counter + 1) %notin% combos[combo])) %>%
            # Remove chosen positions to allow for summarizing with the subtract function 
            select(-c(combos[combo])) %>%
            summarise(avg = subtract(avg))
          # Artifically add in columns with Xs to account for where the positions were removed
          
          # Adds the x back the combo you just took out
          curr_work = add_column(curr_work, past_var = rep("x", dim(curr_work)[1]), .before = combos[combo], .name_repair = "minimal")
          
          for(j in 1:length(this[i,])) {
            curr_work = add_column(curr_work, place_holder = rep("x", dim(curr_work)[1]), .before = this[i,j], .name_repair = "minimal")
            names(curr_work)[which(names(curr_work) == "place_holder")] = j
          }
          
          curr_work = curr_work %>%
            unite("genotype", 1:places, sep = "") %>%
            mutate(genotype = str_replace_all(genotype, "-1", "0"))
          
          next_work = bind_rows(next_work, curr_work)
          
        }
      }
      
      next_work = next_work %>% distinct(genotype, .keep_all = T)
      next_work$avg = next_work$avg ## I used to divide by 2 to make it compatible with linear model but I removed that
      int_work = next_work
      
      next_work = tibble()  
      
    }
    
    final_work = bind_rows(final_work, int_work)
    
  }
  
  
  #### Determine positional distribution and plot
  
  final_work
  
  split_genos = str_split(final_work$genotype, "")
  
  format_pos = c()
  muts = c()
  for(each in 1:length(split_genos)) {
    muts = c(muts, length(which(split_genos[[each]] == "x")))
    format_pos = c(format_pos, paste(position_names[which(split_genos[[each]] == "x")], collapse = "|"))
  }
  
  final_work$pos = format_pos
  final_work$mutations = muts
  
  first_order = final_work %>%
    filter(mutations == 1) %>%
    mutate(pos = factor(pos, levels = unique(pos)),
           mutations = paste("Order", mutations)) %>%
    ggplot(aes(x = pos, y = avg, label = genotype)) +
    geom_boxplot(coef = 6) +
    geom_jitter(width = 0.2) +
    geom_hline(yintercept = 0, lwd = 0.5) +
    geom_hline(yintercept = log10(1.5), lwd = 0.5, lty = 2) +
    geom_hline(yintercept = log10(1/1.5), lwd = 0.5, lty = 2) +
    theme_classic() +
    labs(x = "Position",
         y = expression("Functional Contribution ("*Delta*italic("F")*")")) +
    stat_summary(fun=mean, colour="darkred", geom="crossbar", width=0.2) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(axis.line = element_line(size = 0.2, color = "black"), 
          axis.ticks = element_line(size = 0.2, color = "black"), 
          text = element_text(size = 9), 
          axis.text = element_text(size = 8, color = "black"))
  
  ggsave("first_order_box_plot.svg", plot = first_order, width = 180, height = 247/1.5, units = "mm")
  
  higher_order = final_work %>%
    filter(mutations > 1) %>%
    mutate(pos = factor(pos, levels = unique(pos)),
           mutations = paste("Order", mutations)) %>%
    ggplot(aes(x = pos, y = avg, label = genotype)) +
    geom_boxplot(coef = 6) +
    geom_jitter(width = 0.2, size = 1) +
    geom_hline(yintercept = 0, lwd = 0.5) +
    geom_hline(yintercept = log10(1.5), lwd = 0.5, lty = 2) +
    geom_hline(yintercept = log10(1/1.5), lwd = 0.5, lty = 2) +
    facet_wrap(~ mutations, scales = "free_x") +
    theme_classic() +
    labs(x = "Combination",
         y = expression("Epistasis ("*epsilon*")")) +
    stat_summary(fun=mean, colour="darkred", geom="crossbar", width = 0.2) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(axis.line = element_line(size = 0.2, color = "black"), 
          axis.ticks = element_line(size = 0.2, color = "black"), 
          text = element_text(size = 10), 
          axis.text = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 5))
  
  ggsave("higher_order_box_plot.svg", plot = higher_order, width = 180, height = 247/1.5, units = "mm")
  
  #p = ggplotly(higher_order, tooltip = c("label", "avg"))
  #saveWidget(p, "higher_order_box_widget.html")
  
  ## Idiosyncrasy index
  
  idio_index = final_work %>%
    mutate(pos = factor(pos, levels = unique(pos)),
           mutations = paste("Order", mutations)) %>%
    group_by(pos, mutations) %>%
    summarise(idiosync = 2*sd(avg)) %>%
    na.exclude()
  
  write_csv(idio_index, paste("idio_df_", cond, "_", str_replace(file, ".csv", ""), ".csv", sep = ""))
  
  final_work <<- final_work
  idio_index <<- idio_index
  
}

## Analysis

# Make Empty Output directory if it doesn't already exist

if(!dir.exists("Output")) {
  dir.create("Output")
}

for(cur_input in 1:length(inputs)) {
  # Reset back to home directory
  setwd("/Users/karolbuda/OneDrive - UBC/PhD/Epistasis_Lit/__VersionCtrl__/epistasis_analysis_v2")
  file = inputs[cur_input]
  
  if(dir.exists(paste(c("./Output/", str_remove(file, ".csv")), collapse = "")) & !updater) {
    next
  }
  
  enzyme_name = str_split(str_remove(file, ".csv"), "_")[[1]][1]
  measure_type = str_split(str_remove(file, ".csv"), "_")[[1]][2]
  cond = str_split(str_remove(file, ".csv"), "_")[[1]][3]
  
  outputs = data_loading()
  if(!dir.exists(paste(c("./Output/", str_remove(file, ".csv")), collapse = ""))) {
    dir.create(paste(c("./Output/", str_remove(file, ".csv")), collapse = ""))
  }
  setwd(paste(c("./Output/", str_remove(file, ".csv")), collapse = ""))
  epistasis_analysis()
  blanket_analysis()
  network_analysis()
  network_analysis("WT")
  network_analysis("ratio")
  feed_forward()
  higher_order_box()
  
  net_codes$enz = rep(enzyme_name, dim(net_codes)[1])
  write.csv(net_codes, "observed_values.csv", row.names = F)
  
  ### Diminishing returns pattern output
  
  transition_df = tibble()
  starting_points = unique(obs_links$from)
  
  for(i in 1:length(starting_points)) {
    jumps = obs_links %>% 
      filter(from == starting_points[i]) %>%
      pull(magnitude)
    
    transition_df = rbind(transition_df, tibble(mean = mean(jumps), max = max(jumps)))
  }
  
  transition_df$obs = pred_df[match(starting_points, pred_df$numbers),]$`observed effect`
  
  transition_df$geno = starting_points
  
  transition_df = transition_df %>%
    mutate(enzyme_name = enzyme_name,
           measure_type = measure_type,
           cond = cond)
  
  write_csv(transition_df, paste("transition_df_", cond, "_", str_replace(file, ".csv", ""), ".csv", sep = ""))
  
  ### Network output
  
  if(exists("accessible_path")) {
    
    if(all(derived_likely_path %in% accessible_path)) {
      
    } else {
      write("Accessible path is NOT to the last variant", file = "last_variant.txt")
    }
    
    epi_in_traj = final_work[!str_detect(final_work$genotype, "1"), ] %>%
      mutate(likely = ifelse(genotype %in% str_replace_all(accessible_path, "1", "x"), TRUE, FALSE),
             enzyme_name = enzyme_name,
             measure_type = measure_type,
             cond = cond)
    
    write_csv(epi_in_traj, paste("traj_epi_", cond, "_", str_replace(file, ".csv", ""), ".csv", sep = ""))
    
  }
  
  ## Feed-Forward Output 
  
  ff_prediction_df$`Enzyme` = enzyme_name
  ff_prediction_df$`Linear MAE` = rolling_mae[1:(dim(ff_prediction_df)[1])]
  ff_prediction_df$`Linear Last MAE` = rolling_last_mae[1:(dim(ff_prediction_df)[1])]
  colnames(ff_prediction_df) = c("MAE", "Last MAE", "Step", "Enzyme", "Linear MAE", "Linear Last MAE")
  ff_prediction_df = ff_prediction_df %>% na.omit() 
  
  write_csv(ff_prediction_df, paste("aic_df_", cond, "_", str_replace(file, ".csv", ""), ".csv", sep = ""))
  write_csv(final_work, paste("higher_box_df_", cond, "_", str_replace(file, ".csv", ""), ".csv", sep = ""))
}

##### FOR DEBUGGING COMMENT NEXT BIT OUT