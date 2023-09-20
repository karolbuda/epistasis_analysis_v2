# R 4.1.2

library(tidyverse)
## ✔ dplyr     1.1.3     ✔ readr     2.1.4
## ✔ forcats   1.0.0     ✔ stringr   1.5.0
## ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
## ✔ purrr     1.0.2
library(writexl) # 1.4.2

###

### Supplementary 1 is just a xlsx export of fit_land

fit_land = do.call("rbind", lapply(paste0("Non-linear Transformation/trans_Output/", list.files(path = "Non-linear Transformation/trans_Output", pattern="observed_values.csv", recursive = T)), read_csv, skip = 1, show_col_types=F, col_names = c("id", "effect", "muts", "junk1", "junk2", "enz")))

## Make fit_land have IDs by finding where muts are 0, i.e. start of new dataset, and multiplying
## Those values by the folder names to generate ID vector

fit_land_spots = c(which(fit_land$muts == 0), (dim(fit_land)[1] + 1)) # Add end of dataset
fit_land_counts = c()

for(i in 1:(length(fit_land_spots) - 1)) {
  
  if(i == 1) {
    fit_land_counts = c(fit_land_counts, fit_land_spots[i + 1] - 1)
  } else {
    fit_land_counts = c(fit_land_counts, fit_land_spots[i + 1] - fit_land_spots[i])
  }
  
}

## This assumes folders are read in the order that they are listed

fit_land_ids = rep(str_remove_all(list.files(path = "Non-linear Transformation/trans_Output", pattern="observed_values.csv", recursive = T), "/observed_values.csv"), fit_land_counts)

fit_land$unique_id = fit_land_ids

fit_land = fit_land %>% 
  filter(muts != 0) %>% 
  dplyr::select(-c("junk1", "junk2"))

names(fit_land) = c("Genotype", "F", "Mutations", "Enzyme", "Unique ID")

writexl::write_xlsx(fit_land, "Supplementary File 1.xlsx")

### Supplementary 2 is just a xlsx export of d1

d = do.call("rbind", lapply(paste0("Non-linear Transformation/trans_Output/", list.files(path = "Non-linear Transformation/trans_Output", pattern="*2d_box.csv", recursive = T)), read_csv, show_col_types=F, col_names = c("positions", "identity", "mut", "effects", "enzyme", "type", "cond")))

d1 = d %>% filter_all(any_vars(!is.na(.))) # removes rows with all NA before making unique ID

d1 = d1 %>%
  unite("unique_id", c(positions, enzyme, type, cond), remove = F)

d1 = d1 %>%
  unite("partial_id", c(enzyme, type, cond), remove = F) %>%
  dplyr::select(-c("type", "cond"))

names(d1) = c("Unique ID", "Position", "Genotype", "Mutations", "SME", "Partial ID", "Enzyme")

writexl::write_xlsx(d1, "Supplementary File 2.xlsx")

### Supplementary 3 is just a xlsx export of higher_df with SMEs removed

higher_order_list = list.files(path = "Non-linear Transformation/trans_Output", pattern="higher_box", recursive = T)

higher_df = data.frame()
for(higher_file in 1:length(higher_order_list)){
  Condition = str_split(str_split(higher_order_list[higher_file], "_")[[1]][3], "/")[[1]][1]
  Measurement = str_split(higher_order_list[higher_file], "_")[[1]][2]
  Enzyme = str_split(higher_order_list[higher_file], "_")[[1]][1]
  appending_file = read_csv(paste0("Non-linear Transformation/trans_Output/", higher_order_list[higher_file]), show_col_types=F)
  appending_file$Condition = Condition
  appending_file$Measurement = Measurement
  appending_file$Enzyme = Enzyme
  higher_df = bind_rows(higher_df, appending_file)
}


higher_df = higher_df %>%
  filter(mutations > 1) %>%
  unite("unique_id", c(pos, Condition, Measurement, Enzyme), remove = F) %>%
  unite("partial_id", c(Enzyme, Measurement, Condition), remove = F) %>%
  dplyr::select(-c(Measurement, Condition))

names(higher_df) = c("Genotype", "Epistasis", "Unique ID", "Combination", "Partial ID", "Enzyme")

writexl::write_xlsx(higher_df, "Supplementary File 3.xlsx")
