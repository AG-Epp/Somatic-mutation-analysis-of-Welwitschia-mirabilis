library(dplyr)
library(ggplot2)

### Distance of variant sites to each-other ###


# import bed files containing variant positions
R010_WW01_var_sites <- read.table("Welwitschia output files/240508.R010.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all_Chr.allsites.no_repeats.WW01_variant_sites_passed_filters.bed") 
R010_WW04_var_sites <- read.table("Welwitschia output files/240508.R010.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all_Chr.allsites.no_repeats.WW04_variant_sites_passed_filters.bed") 
R010_ME01_var_sites <- read.table("Welwitschia output files/240508.R010.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all_Chr.allsites.no_repeats.ME01_variant_sites_passed_filters.bed") 

R011_WW01_var_sites <- read.table("Welwitschia output files/240508.R011.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all_Chr.allsites.WW01_variant_sites_passed_filters.bed")
R011_WW04_var_sites <- read.table("Welwitschia output files/240508.R011.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all_Chr.allsites.WW04_variant_sites_passed_filters.bed")
R011_ME01_var_sites <- read.table("Welwitschia output files/240508.R011.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all_Chr.allsites.ME01_variant_sites_passed_filters.bed")

# calculate distance to previous variant
R010_WW01_var_sites$Distance <- with(R010_WW01_var_sites, ave(R010_WW01_var_sites$V2, R010_WW01_var_sites$V1, FUN = function(x) c(NA, diff(x))))
R010_WW04_var_sites$Distance <- with(R010_WW04_var_sites, ave(R010_WW04_var_sites$V2, R010_WW04_var_sites$V1, FUN = function(x) c(NA, diff(x))))
R010_ME01_var_sites$Distance <- with(R010_ME01_var_sites, ave(R010_ME01_var_sites$V2, R010_ME01_var_sites$V1, FUN = function(x) c(NA, diff(x))))
R011_WW01_var_sites$Distance <- with(R011_WW01_var_sites, ave(R011_WW01_var_sites$V2, R011_WW01_var_sites$V1, FUN = function(x) c(NA, diff(x))))
R011_WW04_var_sites$Distance <- with(R011_WW04_var_sites, ave(R011_WW04_var_sites$V2, R011_WW04_var_sites$V1, FUN = function(x) c(NA, diff(x))))
R011_ME01_var_sites$Distance <- with(R011_ME01_var_sites, ave(R011_ME01_var_sites$V2, R011_ME01_var_sites$V1, FUN = function(x) c(NA, diff(x))))

# calculate distance to next variant
R010_WW01_var_sites$Distance_to_next <- with(R010_WW01_var_sites, ave(R010_WW01_var_sites$V2, R010_WW01_var_sites$V1, FUN = function(x) c(diff(x), NA)))
R010_WW04_var_sites$Distance_to_next <- with(R010_WW04_var_sites, ave(R010_WW04_var_sites$V2, R010_WW04_var_sites$V1, FUN = function(x) c(diff(x), NA)))
R010_ME01_var_sites$Distance_to_next <- with(R010_ME01_var_sites, ave(R010_ME01_var_sites$V2, R010_ME01_var_sites$V1, FUN = function(x) c(diff(x), NA)))
R011_WW01_var_sites$Distance_to_next <- with(R011_WW01_var_sites, ave(R011_WW01_var_sites$V2, R011_WW01_var_sites$V1, FUN = function(x) c(diff(x), NA)))
R011_WW04_var_sites$Distance_to_next <- with(R011_WW04_var_sites, ave(R011_WW04_var_sites$V2, R011_WW04_var_sites$V1, FUN = function(x) c(diff(x), NA)))
R011_ME01_var_sites$Distance_to_next <- with(R011_ME01_var_sites, ave(R011_ME01_var_sites$V2, R011_ME01_var_sites$V1, FUN = function(x) c(diff(x), NA)))

# add RunID
R010_WW01_var_sites$RunID <- "R010"
R010_WW04_var_sites$RunID <- "R010"
R010_ME01_var_sites$RunID <- "R010"
R011_WW01_var_sites$RunID <- "R011"
R011_WW04_var_sites$RunID <- "R011"
R011_ME01_var_sites$RunID <- "R011"

# add dataset
R010_WW01_var_sites$dataset <- "no_repeats"
R010_WW04_var_sites$dataset <- "no_repeats"
R010_ME01_var_sites$dataset <- "no_repeats"
R011_WW01_var_sites$dataset <- "all_sites"
R011_WW04_var_sites$dataset <- "all_sites"
R011_ME01_var_sites$dataset <- "all_sites"

# add individual
R010_WW01_var_sites$Indiv <- "WW01"
R010_WW04_var_sites$Indiv <- "WW04"
R010_ME01_var_sites$Indiv <- "ME01"
R011_WW01_var_sites$Indiv <- "WW01"
R011_WW04_var_sites$Indiv <- "WW04"
R011_ME01_var_sites$Indiv <- "ME01"


# Combine all six dataframes into one
combined_df <- rbind(R010_WW01_var_sites, R010_WW04_var_sites, R010_ME01_var_sites,
                         R011_WW01_var_sites, R011_WW04_var_sites, R011_ME01_var_sites)

# Calculate breaks based on the combined dataframe
breaks_2 <- c(0, 10, 150, 1000, Inf)  # Adjust as needed
labels_2 <- c("0-10", "11-150", "> 1000")

# Calculate distance distribution and percentage for each dataframe
combined_df_breaks_2 <- combined_df %>%
  mutate(Distance_to_prev_Bin = cut(Distance, breaks = breaks_2)) %>%
  group_by(Indiv, dataset, Distance_to_prev_Bin) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

combined_df_breaks_2$Indiv <- factor(combined_df_breaks_2$Indiv, levels =c("WW01", "WW04", "ME01"))
combined_df_breaks_2$dataset <- factor(combined_df_breaks_2$dataset, levels=c("no_repeats", "all_sites"))


######### Distance to previous and next variant combined ################

# Calculate the bins for both Distance to previous and Distance to next
combined_df_bothways <- combined_df %>%
  mutate(Distance_to_prev_Bin = cut(Distance, breaks = breaks_2),
         Distance_to_next_Bin = cut(Distance_to_next, breaks = breaks_2))

# Combine bins based on the condition: if at least one of the distances is in a lower bin
combined_df_conditional <- combined_df_bothways %>%
  rowwise() %>%
  mutate(Combined_Bin = case_when(
    Distance_to_prev_Bin == "(0,10]" | Distance_to_next_Bin == "(0,10]" ~ "(0,10]",
    Distance_to_prev_Bin == "(10,150]" | Distance_to_next_Bin == "(10,150]" ~ "(10,150]",
    Distance_to_prev_Bin == "(150,1000]" | Distance_to_next_Bin == "(150,1000]" ~ "(150,1000]",
    TRUE ~ "(1000,Inf]"
  )) %>%
  ungroup()

# Calculate total count for each Indiv and dataset combination
total_counts <- combined_df_conditional %>%
  group_by(Indiv, dataset) %>%
  summarise(Total_Count = n(), .groups = 'drop')


# Group by Indiv, dataset, and the combined bin, then summarize
combined_df_summary <- combined_df_conditional %>%
  group_by(Indiv, dataset, Combined_Bin) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  left_join(total_counts, by = c("Indiv", "dataset")) %>%
  mutate(Percentage = (Count / Total_Count) * 100) %>%
  select(-Total_Count)


combined_df_summary$Indiv <- factor(combined_df_summary$Indiv, levels =c("WW01", "WW04", "ME01"))
combined_df_summary$dataset <- factor(combined_df_summary$dataset, levels=c("no_repeats", "all_sites"))

# Plot the distance distribution
ggplot(combined_df_summary %>% filter(!is.na(Combined_Bin)), 
       aes(x = as.factor(Combined_Bin), y = Percentage, fill = Combined_Bin)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = labels_2)+
  ylim(0, 90) + # Set y-axis limits
  geom_text(aes(label = paste(round(Percentage, 2), "%")), 
            position = position_stack(vjust = 1.0),  # Center vertically
            size = 3.5, 
            color = "black",
            vjust = -0.5) +  # Adjust vertical alignment (subtract a small value)
  labs(x = "Distance [bp] to closest variant", y = "% of variants") +
  theme_minimal() +
  facet_grid(Indiv ~ dataset, scales = "free" )+
  scale_fill_brewer(palette = "Set2")


####### log breaks ###########

library(dplyr)
library(ggplot2)

# Define logarithmic breaks and corresponding labels
log_breaks <- c(0, 10, 100, 1000, 10000, 100000, Inf)
bin_labels <- paste0("(", log_breaks[-length(log_breaks)], ",", log_breaks[-1], "]")

# Calculate the bins for both Distance to previous and Distance to next
combined_df_bothways_log <- combined_df %>%
  mutate(Distance_to_prev_Bin = cut(Distance, breaks = log_breaks, labels = bin_labels, include.lowest = TRUE),
         Distance_to_next_Bin = cut(Distance_to_next, breaks = log_breaks, labels = bin_labels, include.lowest = TRUE))

# Combine bins based on the condition: if at least one of the distances is in a lower bin
combined_df_conditional_log <- combined_df_bothways_log %>%
  rowwise() %>%
  mutate(Combined_Bin = case_when(
    Distance_to_prev_Bin == "(0,10]" | Distance_to_next_Bin == "(0,10]" ~ "(0,10]",
    Distance_to_prev_Bin == "(10,100]" | Distance_to_next_Bin == "(10,100]" ~ "(10,100]",
    Distance_to_prev_Bin == "(100,1000]" | Distance_to_next_Bin == "(100,1000]" ~ "(100,1000]",
    Distance_to_prev_Bin == "(1000,10000]" | Distance_to_next_Bin == "(1000,10000]" ~ "(1000,10000]",
    Distance_to_prev_Bin == "(10000,100000]" | Distance_to_next_Bin == "(10000,100000]" ~ "(10000,100000]",
    TRUE ~ "(100000,Inf]"
  )) %>%
  ungroup()

# Calculate total count for each Indiv and dataset combination
total_counts_log <- combined_df_conditional_log %>%
  group_by(Indiv, dataset) %>%
  summarise(Total_Count = n(), .groups = 'drop')

# Group by Indiv, dataset, and the combined bin, then summarize
combined_df_summary_log <- combined_df_conditional_log %>%
  group_by(Indiv, dataset, Combined_Bin) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  left_join(total_counts_log, by = c("Indiv", "dataset")) %>%
  mutate(Percentage = (Count / Total_Count) * 100) %>%
  select(-Total_Count)

# View the resulting dataframe
print(combined_df_summary_log)

combined_df_summary_log$Indiv <- factor(combined_df_summary_log$Indiv, levels =c("WW01", "WW04", "ME01"))
combined_df_summary_log$dataset <- factor(combined_df_summary_log$dataset, levels=c("no_repeats", "all_sites"))

# Plot the distance distribution
ggplot(combined_df_summary_log, 
       aes(x = Combined_Bin, y = Percentage, fill = Combined_Bin)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")), 
            position = position_stack(vjust = 1.0),  # Center vertically
            size = 3.5, 
            color = "black",
            vjust = -0.5) +  # Adjust vertical alignment (subtract a small value)
  labs(x = "Distance [bp] to closest variant", y = "% of variants", title = "Distance Distribution") +
  theme_minimal() +
  facet_grid(Indiv ~ dataset, scales = "free_y") +
  scale_y_continuous(limits = c(0, 90))  # Adjust y-axis limits if necessary

