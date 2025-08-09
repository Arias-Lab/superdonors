#!/usr/bin/env Rscript
# ------------------------------------------
# Date:        2025-08-05
# Description: This script perform summary and Euler plot analysis
# ------------------------------------------

#execute preferably after E00_dataset_summary.R

#output files
#counts for the euler diagram
#001_euler_diagram_counts.txt

#proportional euler diagram with eulerr package
#002_proportional_euler_diagram_homozygous.pdf

#input files
filtered_data<-read.table(file="anonymized_dataset_this_study.txt", header=TRUE)

#all output files will go there
#create directory only if it does not exist
#YYYYMMDD change to relevant value

if (!dir.exists("YYYYMMDD_output_files")) {
  dir.create("YYYYMMDD_output_files")
}

#required librarires
#proportional euler diagram with Eulerr
if (!requireNamespace("eulerr", quietly = TRUE)) {
  install.packages("eulerr")
}

#libraries
library(eulerr)

# Calculate counts for each condition
euler_counts <- c(
  # HLAA1 == HLAA2 AND HLAB1 != HLAB2 AND HLAC1 != HLAC2
  "HLA-A" = sum(filtered_data$HLAA1 == filtered_data$HLAA2 & 
                             filtered_data$HLAB1 != filtered_data$HLAB2 & 
                             filtered_data$HLAC1 != filtered_data$HLAC2),
  # HLAA1 != HLAA2 AND HLAB1 == HLAB2 AND HLAC1 != HLAC2
  "HLA-B" = sum(filtered_data$HLAA1 != filtered_data$HLAA2 & 
                       filtered_data$HLAB1 == filtered_data$HLAB2 & 
                       filtered_data$HLAC1 != filtered_data$HLAC2),
  # HLAA1 != HLAA2 AND HLAB1 != HLAB2 AND HLAC1 == HLAC2
  "HLA-C" = sum(filtered_data$HLAA1 != filtered_data$HLAA2 & 
                       filtered_data$HLAB1 != filtered_data$HLAB2 & 
                       filtered_data$HLAC1 == filtered_data$HLAC2),
  # HLAA1 == HLAA2 AND HLAB1 == HLAB2 AND HLAC1 != HLAC2
  "HLA-A&HLA-B" = sum(filtered_data$HLAA1 == filtered_data$HLAA2 & 
                       filtered_data$HLAB1 == filtered_data$HLAB2 & 
                       filtered_data$HLAC1 != filtered_data$HLAC2),
  # HLAA1 == HLAA2 AND HLAB1 != HLAB2 AND HLAC1 == HLAC2
  "HLA-A&HLA-C" = sum(filtered_data$HLAA1 == filtered_data$HLAA2 & 
                       filtered_data$HLAB1 != filtered_data$HLAB2 & 
                       filtered_data$HLAC1 == filtered_data$HLAC2),
  # HLAA1 != HLAA2 AND HLAB1 == HLAB2 AND HLAC1 == HLAC2
  "HLA-B&HLA-C" = sum(filtered_data$HLAA1 != filtered_data$HLAA2 & 
                       filtered_data$HLAB1 == filtered_data$HLAB2 & 
                       filtered_data$HLAC1 == filtered_data$HLAC2),
  # HLAA1 == HLAA2 AND HLAB1 == HLAB2 AND HLAC1 == HLAC2
  "HLA-A&HLA-B&HLA-C" = sum(filtered_data$HLAA1 == filtered_data$HLAA2 & 
                                 filtered_data$HLAB1 == filtered_data$HLAB2 & 
                                 filtered_data$HLAC1 == filtered_data$HLAC2)
)

# Convert to data frame
euler_df <- data.frame(
  Category = names(euler_counts),
  Count = as.integer(euler_counts)
)

#export summary data
write.table(euler_df, file = "YYYYMMDD_output_files/001_euler_diagram_counts.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

euler_plot <- euler(euler_counts)

# Save the plot as a PDF
pdf("YYYYMMDD_output_files/002_proportional_euler_diagram_homozygous.pdf", width = 8, height = 8)
plot(euler_plot,
     fills = c("blue", "green", "red"),
     alpha = 0.5,
     labels = list(cex = 1.5),
     quantities = list(cex = 1.5)
)
dev.off()

rm(euler_counts, euler_plot, euler_df, filtered_data)
