R
# Load required libraries
library(dplyr)       # For data manipulation
library(ggplot2)     # For creating visualizations
library(ggnewscale)  # For adding multiple scales for color or fill in ggplot2

# Set working directory (Add the path here if necessary)

# Read baseline data
data <- read.table("RDABaseforplot.txt", header = TRUE, sep = "\t")  # Load the main RDA dataset
data$CHROM <- as.factor(data$CHR)  # Convert chromosome column to factor for categorical handling

# Read outliers LFMM and RDA data
outliers_lfmm <- read.table("BL_LMM_OUTLIERS_POS_CHROM_SPLIT.txt", header = TRUE, sep = "\t")  # LFMM outliers
str(outliers_lfmm)  # Check structure of LFMM outliers data

outliers_rda <- read.table("BL_RDA_OUTLIERS_POS_CHROM_SPLIT.txt", header = TRUE, sep = "\t")  # RDA outliers
str(outliers_rda)  # Check structure of RDA outliers data

# Convert CHROM columns to factors in both datasets
outliers_lfmm$CHROM <- as.factor(outliers_lfmm$CHROM)
outliers_rda$CHROM <- as.factor(outliers_rda$CHROM)

# Create a merged column to uniquely identify positions in all datasets
data$MergedColumn <- paste(data$CHROM, data$POS, sep = "_")
outliers_lfmm$MergedColumn <- paste(outliers_lfmm$CHROM, outliers_lfmm$POS, sep = "_")
outliers_rda$MergedColumn <- paste(outliers_rda$CHROM, outliers_rda$POS, sep = "_")

# Identify common and unique outliers
common_outliers <- inner_join(outliers_lfmm, outliers_rda, by = c("CHROM", "POS"))  # Find common outliers
str(common_outliers)  # Check structure of common outliers data

common_outliers$MergedColumn <- paste(common_outliers$CHROM, common_outliers$POS, sep = "_")  # Add merged column

# Identify unique outliers in each dataset
unique_lfmm <- outliers_lfmm$MergedColumn[!(outliers_lfmm$MergedColumn %in% outliers_rda$MergedColumn)]
str(unique_lfmm)  # Check structure of unique LFMM outliers
unique_rda <- outliers_rda$MergedColumn[!(outliers_rda$MergedColumn %in% outliers_lfmm$MergedColumn)]
str(unique_rda)  # Check structure of unique RDA outliers

# Assign colors to points based on outlier membership
data$color <- ifelse(data$MergedColumn %in% common_outliers$MergedColumn, "common",
                     ifelse(data$MergedColumn %in% unique_lfmm, "lfmm", "grey"))

# Check distribution of assigned color categories
print(table(data$color))  # Frequency of each color category
str(data)

# Merge grey points with CHROM to differentiate
factor_to_merge <- "grey"
data$color <- ifelse(data$color == factor_to_merge, paste(data$color, data$CHROM), data$color)

# Filter data to include only specific chromosomes
chromosome_pattern <- paste(c("OZ075377.1", "OZ075378.1", "OZ075379.1", "OZ075380.1", 
                              "OZ075381.1", "OZ075382.1", "OZ075383.1", "OZ075384.1", 
                              "OZ075385.1", "OZ075386.1", "OZ075387.1", "OZ075388.1", 
                              "OZ075389.1", "OZ075390.1", "OZ075391.1", "OZ075392.1", 
                              "OZ075393.1", "OZ075394.1"), collapse="|")
data <- subset(data, grepl(chromosome_pattern, CHROM))

# Recode chromosome names into numeric values
data <- data %>%
  mutate(CHROM = recode(CHROM,
                        "OZ075377.1"=1, "OZ075378.1"=2, "OZ075379.1"=3, "OZ075380.1"=4, 
                        "OZ075381.1"=5, "OZ075382.1"=6, "OZ075383.1"=7, "OZ075384.1"=8,
                        "OZ075385.1"=9, "OZ075386.1"=10, "OZ075387.1"=11, "OZ075388.1"=12, 
                        "OZ075389.1"=13, "OZ075390.1"=14, "OZ075391.1"=15, "OZ075392.1"=16, 
                        "OZ075393.1"=17, "OZ075394.1"=18))

# Read additional selection sweep data from selscan
sweepsppehh <- read.table("Xpehh_outliers.txt", header = TRUE, sep = "\t")
str(sweepsppehh)  # Check structure of sweep data

# Recode chromosomes in sweep data
sweepsppehh <- sweepsppehh %>%
  mutate(chrom = recode(chrom,
                        "OZ075377.1"=1, "OZ075378.1"=2, "OZ075379.1"=3, "OZ075380.1"=4,
                        "OZ075381.1"=5, "OZ075382.1"=6, "OZ075383.1"=7, "OZ075384.1"=8, 
                        "OZ075385.1"=9, "OZ075386.1"=10, "OZ075387.1"=11, "OZ075388.1"=12, 
                        "OZ075389.1"=13, "OZ075390.1"=14, "OZ075391.1"=15, "OZ075392.1"=16, 
                        "OZ075393.1"=17, "OZ075394.1"=18))

# Create a dataframe for midpoints of sweep windows
sweeps_dfXPEHH <- data.frame(
  x = (sweepsppehh$BP1 + sweepsppehh$BP2) / 2,  # Midpoint of sweep regions
  y = 0,  # Set y-value for visual separation
  CHROM = sweepsppehh$chrom  # Corresponding chromosome
)

# Visualize data with ggplot2
p <- ggplot(data, aes(x = POS, y = ABSRDA1, color = color)) +
  geom_segment(data = sweeps_dfXPEHH, aes(x = x, xend = x, y = -Inf, yend = Inf), 
               inherit.aes = FALSE, color = "blue", size = 0.5, alpha = 0.5) +  # Background blue lines for sweeps
  geom_point(size=0.5) +  # Main points
  scale_color_manual(values = c("lfmm" = "#fee090", "common" = "red", 
                                "grey OZ075377.1" = "#bababa", "grey OZ075378.1" = "#878787", 
                                # Repeat for other grey categories
                                "grey OZ075394.1" = "#878787")) +
  labs(x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 16), panel.grid = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  facet_grid(. ~ CHROM, scales = "free_x", space = "free_x", switch = "x") +
  scale_x_discrete(labels = NULL) +
  geom_point(data = subset(data, color == "lfmm"), aes(x = POS, y = ABSRDA1), color = "#fee090", size = 0.5)

# Save the plot
ggsave("manhattan_plot.jpeg", p, width = 10, height = 3, dpi = 600)  # Save the plot as a JPEG
