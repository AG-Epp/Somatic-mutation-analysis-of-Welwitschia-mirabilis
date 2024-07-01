library(ggplot2)
library(ggbreak)

################### ON THE SERVER - source files not available on GitHub ###################
# Visualization of GQ and RGQ
Chr21_VarToTab_GQ_240428 <- read.table("240428.GATK.VarToTab.all_samples.welwitschia.Chr21.genotypeGVCFs.allsites.GF_GQ.table", header = TRUE)

Chr21_VarToTab_RGQ_240428 <- read.table("240428.GATK.VarToTab.all_samples.welwitschia.Chr21.genotypeGVCFs.allsites.GF_RGQ.table", header = TRUE)

Chr21_VarToTab_DP_240428 <- read.table("240428.GATK.VarToTab.all_samples.welwitschia.Chr21.genotypeGVCFs.allsites.GF_DP.table", header = TRUE)

#NOTE: I am using sample ME01_1 for all analyses

# Create a histogram
#hist_obj <- hist(data)
#GQ Histogram: done and dusted
Chr21_GQ_hist <-hist(Chr21_VarToTab_GQ_240428$ME01_1.GQ, breaks = 100, na.rm = TRUE)

#RGQ histogram: needs adjustment if data distribution
Chr21_RGQ_hist <- hist(Chr21_VarToTab_RGQ_240428$ME01_1.RGQ, breaks = 100, na.rm = TRUE)


 range(Chr21_VarToTab_RGQ_240428$ME01_1.RGQ, na.rm = TRUE)
[1]    0 6694 
 
 
# Ensure that the data is numeric
Chr21_VarToTab_ME01_1.RGQ.num <- as.numeric(Chr21_VarToTab_RGQ_240428$ME01_1.RGQ)

# Cap values at 99
Chr21_VarToTab_ME01_1.RGQ.num[Chr21_VarToTab_ME01_1.RGQ.num > 99] <- 99


# Create histogram
Chr21_RGQ_hist <- hist(Chr21_VarToTab_ME01_1.RGQ.num, breaks = 100, na.rm = TRUE)

#DP Histogram

range(Chr21_VarToTab_DP_240428$ME01_1.DP, na.rm = TRUE)
[1]    0 8277

# Ensure that the data is numeric
Chr21_VarToTab_ME01_1.DP.num <- as.numeric(Chr21_VarToTab_DP_240428$ME01_1.DP)

# Cap values at 99
Chr21_VarToTab_ME01_1.DP.num[Chr21_VarToTab_ME01_1.DP.num > 100] <- 100
Chr21_DP_hist <- hist(Chr21_VarToTab_ME01_1.DP.num, breaks = 100, na.rm = TRUE)

# Save an object to a file
saveRDS(object, file = "my_data.rds")

# Restore the object
readRDS(file = "my_data.rds")

saveRDS(Chr21_DP_hist, file = "Chr21_ME01.1_DP_hist.rds")

saveRDS(Chr21_GQ_hist, file = "Chr21_ME01.1_DP_hist.rds")

saveRDS(Chr21_RGQ_hist, file = "Chr21_ME01.1_DP_hist.rds")

#############LOCAL -  files available on GitHub #############

Chr21_ME01.1_DP_hist <- readRDS("GQ-RGQ-DP-histograms/Chr21_ME01.1_DP_hist.rds")
Chr21_ME01.1_GQ_hist <- readRDS("GQ-RGQ-DP-histograms/Chr21_ME01.1_GQ_hist.rds")
Chr21_ME01.1_RGQ_hist <- readRDS("GQ-RGQ-DP-histograms/Chr21_ME01.1_RGQ_hist.rds")


##DP##
# Extract bin boundaries and counts
bin_boundaries_DP <- Chr21_ME01.1_DP_hist$breaks
bin_counts_DP <- Chr21_ME01.1_DP_hist$counts

# Calculate midpoints of bins
bin_midpoints_DP <- (bin_boundaries_DP[-1] + bin_boundaries_DP[-length(bin_boundaries_DP)]) / 2

# Create a data frame
Chr21_ME01.1_DP_hist_df <- data.frame(midpoint_DP = bin_midpoints_DP, count_DP = bin_counts_DP)

# Plot histogram using ggplot2
DP_plot <- ggplot(Chr21_ME01.1_DP_hist_df, aes(x = midpoint_DP, y = count_DP)) +
  geom_bar(stat = "identity") +
  labs(x = "DP", y = "no. of sites")+
  theme_minimal()+
  xlim(0,50)+
  scale_y_break(c(5e+06 , 3.5e+07))+
  scale_color_viridis(discrete=TRUE)

DP_plot

##GQ##
# Extract bin boundaries and counts
bin_boundaries_GQ <- Chr21_ME01.1_GQ_hist$breaks
bin_counts_GQ <- Chr21_ME01.1_GQ_hist$counts

log_bin_counts_GQ <- log(bin_counts_GQ)

# Calculate midpoints of bins
bin_midpoints_GQ <- (bin_boundaries_GQ[-1] + bin_boundaries_GQ[-length(bin_boundaries_GQ)]) / 2


# Create a data frame for GQ
Chr21_ME01.1_GQ_hist_df <- data.frame(midpoint_GQ = bin_midpoints_GQ, count_GQ = bin_counts_GQ)

#plot of GQ
GQ_plot <- ggplot(Chr21_ME01.1_GQ_hist_df, aes(x = midpoint_GQ, y = count_GQ)) +
  geom_bar(stat = "identity") +
  #labs(x = "GQ", y = "Count", title = "Genotype quality distribution, ME01_1, Chr21")+
  labs(x = "GQ", y = "no. of sites")+
  theme_minimal()+
  scale_color_viridis(discrete=TRUE)+
  scale_y_continuous(labels = scales::scientific_format())+
  scale_y_break(c(450000,1000000))


##RGQ
# Extract bin boundaries and counts
bin_boundaries_RGQ <- Chr21_ME01.1_RGQ_hist$breaks
bin_counts_RGQ <- Chr21_ME01.1_RGQ_hist$counts

log_bin_counts_RGQ <- log(bin_counts_RGQ)

# Calculate midpoints of bins
bin_midpoints_RGQ <- (bin_boundaries_RGQ[-1] + bin_boundaries_RGQ[-length(bin_boundaries_RGQ)]) / 2


# Create a data frame for RGQ
Chr21_ME01.1_RGQ_hist_df <- data.frame(midpoint_RGQ = bin_midpoints_RGQ, count_RGQ = bin_counts_RGQ)

#plot of RGQ
RGQ_plot <- ggplot(Chr21_ME01.1_RGQ_hist_df, aes(x = midpoint_RGQ, y = count_RGQ)) +
  geom_bar(stat = "identity") +
  #labs(x = "RGQ", y = "Count", title = "Reference genotype quality distribution, ME01_1, Chr21")+
  labs(x = "RGQ", y = "no. of sites")+
  theme_minimal()+
  scale_color_viridis(discrete=TRUE)+
  scale_y_continuous(labels = scales::scientific_format())+
  scale_y_break(c(7.5e+06, 3e+07))

(RGQ_plot / GQ_plot)


