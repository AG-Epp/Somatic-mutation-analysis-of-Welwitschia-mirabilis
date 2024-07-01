
# load necessary packages
library(ggplot2)
library(dplyr)
library(viridis)


## CALCULATION AND VISUALIZATION OF SMR PER CHROMOSOME

# Import of summarized variant count files
#minGQ = min RGQ = 90
R004 <- read.table("Welwitschia output files/240423.R004.240423_variant_counting_v16_no-bed.all_Chr.allsites.no_repeats.variant_counts.tab")
R005 <- read.table("Welwitschia output files/240423.R005.240423_variant_counting_v16_no-bed.all_Chr.allsites.variant_counts.tab")

#minGQ = minRGQ = 50
R007 <- read.table("Welwitschia output files/240426.R007.240425_variant_counting_v17_yes_bed.all_Chr.allsites.no_repeats.variant_counts.tab")
R008 <- read.table("Welwitschia output files/240426.R008.240425_variant_counting_v17_yes_bed.all_Chr.allsites.variant_counts.tab")

#minRGQ = 40, min GQ = 90
R010 <- read.table("Welwitschia output files/240508.R010.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all-Chr.allsites.no_repeats.variant_counts.tab")
R011 <- read.table("Welwitschia output files/240508.R011.variant_counting_v18_yes_bed_minRGQ40_minGQ90.all-Chr.allsites.variant_counts.tab")

# add a sixth column that calculates the SMR
R004$V6 = R004$V4/R004$V3
R005$V6 = R005$V4/R005$V3
R007$V6 = R007$V4/R007$V3
R008$V6 = R008$V4/R008$V3
R010$V6 = R010$V4/R010$V3
R011$V6 = R011$V4/R011$V3

# rename all columns
names(R004) <- c("CHROM", "INDIV", "callable_sites", "var_sites", "hom_var_sites", "var_sites/callable_sites")
names(R005) <- c("CHROM", "INDIV", "callable_sites", "var_sites", "hom_var_sites", "var_sites/callable_sites")
names(R007) <- c("CHROM", "INDIV", "callable_sites", "var_sites", "hom_var_sites", "var_sites/callable_sites")
names(R008) <- c("CHROM", "INDIV", "callable_sites", "var_sites", "hom_var_sites", "var_sites/callable_sites")


# add run number and filtering parameters to the dataframes
R004$V7 = "R004"
R005$V7 = "R005"
R007$V7 = "R007"
R008$V7 = "R008"
R010$V7 = "R010"
R011$V7 = "R011"

R004$dataset <- "no_repeats"
R005$dataset <- "all_sites"
R007$dataset <- "no_repeats"
R008$dataset <- "all_sites"
R010$dataset <- "no_repeats"
R011$dataset <- "all_sites"


# concatenate to a single dataframe

all_variant_counts <- rbind(R004, R005, R007, R008)

names(all_variant_counts) <- c("CHROM", "INDIV", "callable_sites", "var_sites", "hom_var_sites", "var_sites/callable_sites", "RunID", "dataset" )


all_variant_counts$minRGQ <- sapply(all_variant_counts$RunID, function(run_id) {
  switch(run_id,
         "R004" = "90",
         "R005" = "90",
         "R007" = "50",
         "R008" = "50",
         "R010" = "40",
         "R011" = "40",
         "default_result")})
all_variant_counts$minGQ <- sapply(all_variant_counts$RunID, function(run_id_2) {
  switch(run_id_2,
         "R004" = "90",
         "R005" = "90",
         "R007" = "50",
         "R008" = "50",
         "R010" = "90",
         "R011" = "90",
         "default_result")})

all_variant_counts$INDIV <- factor(all_variant_counts$INDIV, levels = c("WW01", "WW04", "ME01"))
all_variant_counts$dataset <- factor(all_variant_counts$dataset, levels = c("no_repeats", "all_sites"))


all_variant_counts_2 <- rbind(R010, R011)
names(all_variant_counts_2) <- c("CHROM", "INDIV", "callable_sites", "var_sites", "hom_var_sites", "var_sites/callable_sites", "RunID", "dataset" )


all_variant_counts_2$minRGQ = "40"
all_variant_counts_2$minGQ = "90"

all_variant_counts_2$dataset <- factor(all_variant_counts_2$dataset, levels = c("no_repeats", "all_sites"))
all_variant_counts_2$INDIV <- factor(all_variant_counts_2$INDIV, levels = c("WW01", "WW04", "ME01"))


per_chr_SMR_summary <- rbind(all_variant_counts, all_variant_counts_2)

# plot of SMR per chromosome
ggplot(data = per_chr_SMR_summary, aes(x=CHROM, y=var_sites/callable_sites, color=INDIV)) +
  geom_point() +
  theme_minimal() +
  labs(y="Somatic mutation rate", x = "", color="Individual") +
  scale_fill_viridis_b(option = "rainbow") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
  facet_grid(minRGQ + minGQ ~ dataset, scales = "free_y", space = "fixed", labeller=label_both) +
  theme(
    text = element_text(size = 14),  # Font size for main text
    axis.title = element_text(size = 12),  # Font size for axis titles
    axis.text = element_text(size = 12),  # Font size for axis labels
    axis.title.y = element_text(size = 14)  # Font size for y-axis title
  )



