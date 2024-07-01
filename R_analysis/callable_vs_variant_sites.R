
# load necessary packages
library(ggplot2)
library(dplyr)
library(viridis)


# Visualization of callable vs variant sites

# SUMMARIES OF VARIANT COUNTS (I:E: FOR ENTIRE INDIVIDUALS)
#summary tables were made using excel sheets lol why did I do that
R004_summary <- read.table("Variant rates summary files/R004-summary.txt", header = TRUE)
R005_summary <- read.table("Variant rates summary files/R005-summary.txt", header = TRUE)
R007_summary <- read.table("Variant rates summary files/R007-summary.txt", header = TRUE)
R008_summary <- read.table("Variant rates summary files/R008-summary.txt", header = TRUE)
R010_summary <- read.table("Variant rates summary files/R010-summary.txt", header = TRUE)
R011_summary <- read.table("Variant rates summary files/R011-summary.txt", header = TRUE)

# assign RunIDs, filtering criteria and dataset info to be able to plot

R004_summary$RunID <- "R004" 
R005_summary$RunID <- "R005"
R007_summary$RunID <- "R007"
R008_summary$RunID <- "R008"
R010_summary$RunID <- "R010"
R011_summary$RunID <- "R011"


R004_summary$minGQ <- 90
R005_summary$minGQ <- 90
R007_summary$minGQ <- 50
R008_summary$minGQ <- 50

R004_summary$minRGQ <- 90
R005_summary$minRGQ <- 90
R007_summary$minRGQ <- 50
R008_summary$minRGQ <- 50

R004_summary$dataset <- "no_repeats"
R005_summary$dataset <- "all_sites"
R007_summary$dataset <- "no_repeats"
R008_summary$dataset <- "all_sites"
R010_summary$dataset <- "no_repeats"
R011_summary$dataset <- "all_sites"



# combines all summary files (runs R004, R005, R007 and R008)
all_summaries <- rbind(R004_summary, R005_summary,R007_summary, R008_summary)

# puts individuals in order from smallest (youngest) to largest (oldest)
all_summaries$INDIV <- factor(all_summaries$INDIV, levels = c("WW01", "WW04", "ME01"))

# assign levels to datasets for plotting in the desired order
all_summaries$dataset <- factor(all_summaries$dataset, levels = c("no_repeats", "all_sites"))


# combines all summary files (runs R010, R011)

all_summaries_2 <- rbind(R010_summary, R011_summary)

# puts individuals in order from smallest (youngest) to largest (oldest)
all_summaries_2$INDIV <- factor(all_summaries_2$INDIV, levels = c("WW01", "WW04", "ME01"))

# assign levels to datasets for plotting in the desired order
all_summaries_2$dataset <- factor(all_summaries_2$dataset, levels = c("no_repeats", "all_sites"))

all_summaries_2$minGQ <- 90
all_summaries_2$minRGQ <- 40


### plots###

## since R010 and R011 need to be scaled differently from the other runs, plots were created separately and combined using Gimp v. 2.10.36

# plots using the combined dataset
ggplot(data = all_summaries, aes(x=INDIV)) +
  geom_bar(aes(y = callable_sites, fill ="Callable sites"), stat = "identity") +
  geom_bar(aes(y = var_sites *10000, fill ="Variant sites"), stat = "identity") +
  scale_fill_manual(values = c("Callable sites" = "#482878FF", "Variant sites" = "#B4DE2CFF")) +
  labs(x = "Individual", y = "Callable sites", fill = "Legend") +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~./10000, name = "Variant sites")) +   # Secondary y-axis scaled by dividing by 10000
  facet_wrap(dataset ~ minRGQ + minGQ, scales = "free_y", axes = "all", axis.labels = "all", labeller=label_both, nrow = 2)+
  theme(
    text = element_text(size = 12),  # Font size for main text
    axis.title = element_text(size = 14),  # Font size for axis titles
    axis.text = element_text(size = 12),  # Font size for axis labels
    axis.title.y = element_text(size = 14)  # Font size for y-axis title
  )



ggplot(data = all_summaries_2, aes(x=INDIV)) +
  geom_bar(aes(y = callable_sites, fill ="Callable sites"), stat = "identity") +
  geom_bar(aes(y = var_sites *100000, fill ="Variant sites"), stat = "identity") +
  scale_fill_manual(values = c("Callable sites" = "#482878FF", "Variant sites" = "#B4DE2CFF")) +
  labs(x = "Individual", y = "Callable sites", fill = "Legend") +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~./100000, name = "Variant sites")) +   # Secondary y-axis scaled by dividing by 100000
  facet_wrap(dataset ~ minRGQ + minGQ, scales = "free_y", axes = "all", axis.labels = "all", labeller=label_both, nrow = 2)+
  theme(
    text = element_text(size = 12),  # Font size for main text
    axis.title = element_text(size = 14),  # Font size for axis titles
    axis.text = element_text(size = 12),  # Font size for axis labels
    axis.title.y = element_text(size = 14)  # Font size for y-axis title
  )
