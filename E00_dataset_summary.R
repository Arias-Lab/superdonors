#!/usr/bin/env Rscript
# ------------------------------------------
# Date:        2025-08-05
# Description: Script is to summarize HLA structures in the cohort
# ------------------------------------------

#output files
#number of homozygous HLAA, -B or -C
#001_number_of_homozygous_HLAA.txt
#002_number_of_homozygous_HLAB.txt
#003_number_of_homozygous_HLAC.txt

#number of samples with at least one homozygous
#004_number_of_homozygous_samples.txt

#contain unique rows where at least one allele is homozygous
#005_at_least_one_homozygous.txt

#summary of the allele count and frequency in the population
#006_HLAA_population_summary.txt
#007_HLAB_population_summary.txt
#008_HLAC_population_summary.txt

#summary of the allele count and frequency in the population in one table
#009_HLA_population_summary.txt

#core dataset imported
filtered_data<-read.table(file="anonymized_dataset_this_study.txt", header=TRUE)

#all output files will go there
#create directory only if it does not exist
#YYYYMMDD change to relevant value

if (!dir.exists("YYYYMMDD_output_files")) {
  dir.create("YYYYMMDD_output_files")
}

#main

#to summarize the number of rows where HLAA1 == to HLAA2, meaning they are homozygous
#to summarize the number of rows where HLAB1 == to HLAB2, meaning they are homozygous
#to summarize the number of rows where HLAC1 == to HLAC2, meaning they are homozygous
homoA <- sum(filtered_data$HLAA1 == filtered_data$HLAA2)
homoB <- sum(filtered_data$HLAB1 == filtered_data$HLAB2)
homoC <- sum(filtered_data$HLAC1 == filtered_data$HLAC2)

write.table(homoA, file = "YYYYMMDD_output_files/001_number_of_homozygous_HLAA.txt", sep = "\t", row.names = FALSE, col.names = c("HLA-A_homozygous"), quote = FALSE)
write.table(homoB, file = "YYYYMMDD_output_files/002_number_of_homozygous_HLAB.txt", sep = "\t", row.names = FALSE, col.names = c("HLA-B_homozygous"), quote = FALSE)
write.table(homoC, file = "YYYYMMDD_output_files/003_number_of_homozygous_HLAC.txt", sep = "\t", row.names = FALSE, col.names = c("HLA-C_homozygous"), quote = FALSE)

# to identify unique rows where at leaset one allele is homozygous
# Filter rows where at least one condition is met
unique_filtered_rows <- filtered_data[
  filtered_data$HLAA1 == filtered_data$HLAA2 |
  filtered_data$HLAB1 == filtered_data$HLAB2 |
  filtered_data$HLAC1 == filtered_data$HLAC2, ]
# Remove duplicate rows
unique_filtered_rows <- unique(unique_filtered_rows)

write.table(unique_filtered_rows, file = "YYYYMMDD_output_files/004_at_least_one_homozygous.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

number_of_homozygous_samples<-nrow(unique_filtered_rows)
write.table(number_of_homozygous_samples, file = "YYYYMMDD_output_files/005_number_of_homozygous_samples.txt", sep = "\t", row.names = FALSE, col.names = c("homozygous_present"), quote = FALSE)

#creates summary table with count and frequencies for HLA-A -B and -C
#define function for HLA-A -B or -C frequency analysis
summarize_allele_frequencies <- function(data, cols, colnames_out, out_file) {
  #Combine the specified columns into a single vector
  combined_values <- c(data[[cols[1]]], data[[cols[2]]])
  #Calculate the count of each unique string
  frequency_summary <- table(combined_values)
  #Convert the table into a data frame
  frequency_df <- as.data.frame(frequency_summary)
  colnames(frequency_df) <- colnames_out
  #total length (sum of lengths of both columns)
  total_length <- nrow(data) * 2
  frequency_df$freq <- round(frequency_df[[colnames_out[2]]] / total_length, 4)
  # Sort the frequency data frame in descending order of count
  frequency_df <- frequency_df[order(-frequency_df[[colnames_out[2]]]), ]
  # Export the sorted result to a text file
  write.table(frequency_df, file = out_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#perform analysis for HLA-A
summarize_allele_frequencies(data = filtered_data, cols = c("HLAA1", "HLAA2"), colnames_out = c("HLA-A", "allele_count"), out_file = "YYYYMMDD_output_files/006_HLAA_population_summary.txt")

#perform analysis for HLA-B
summarize_allele_frequencies(data = filtered_data, cols = c("HLAB1", "HLAB2"), colnames_out = c("HLA-B", "allele_count"), out_file = "YYYYMMDD_output_files/007_HLAB_population_summary.txt")

#perform analysis for HLA-C
summarize_allele_frequencies(data = filtered_data, cols = c("HLAC1", "HLAC2"), colnames_out = c("HLA-C", "allele_count"), out_file = "YYYYMMDD_output_files/008_HLAC_population_summary.txt")

#merge the population summmary into a single data frame
# Load the data frames
options(scipen = 999)
HLAA <- read.table("YYYYMMDD_output_files/006_HLAA_population_summary.txt", header = TRUE, sep = "\t")
HLAB <- read.table("YYYYMMDD_output_files/007_HLAB_population_summary.txt", header = TRUE, sep = "\t")
HLAC <- read.table("YYYYMMDD_output_files/008_HLAC_population_summary.txt", header = TRUE, sep = "\t")

# Find the maximum number of rows among the three data frames
max_rows <- max(nrow(HLAA), nrow(HLAB), nrow(HLAC))
# Pad each data frame with NA to match the maximum number of rows
HLAA_padded <- rbind(HLAA, data.frame(HLA.A = rep(NA, max_rows - nrow(HLAA)),
                                      allele_count = rep(NA, max_rows - nrow(HLAA)),
                                      freq = rep(NA, max_rows - nrow(HLAA))))
HLAB_padded <- rbind(HLAB, data.frame(HLA.B = rep(NA, max_rows - nrow(HLAB)),
                                      allele_count = rep(NA, max_rows - nrow(HLAB)),
                                      freq = rep(NA, max_rows - nrow(HLAB))))
HLAC_padded <- rbind(HLAC, data.frame(HLA.C = rep(NA, max_rows - nrow(HLAC)),
                                      allele_count = rep(NA, max_rows - nrow(HLAC)),
                                      freq = rep(NA, max_rows - nrow(HLAC))))
# Merge the data frames side by side
merged_side_by_side <- cbind(HLAA_padded, HLAB_padded, HLAC_padded)
colnames(merged_side_by_side) <- c("HLA-A", "count_A", "freq_A", 
                                   "HLA-B", "count_B", "freq_B", 
                                   "HLA-C", "count_C", "freq_C")
# Explicitly format the numeric columns to avoid scientific notation
merged_side_by_side$count_A <- format(merged_side_by_side$count_A, scientific = FALSE)
merged_side_by_side$freq_A <- format(merged_side_by_side$freq_A, scientific = FALSE, nsmall = 4)
merged_side_by_side$count_B <- format(merged_side_by_side$count_B, scientific = FALSE)
merged_side_by_side$freq_B <- format(merged_side_by_side$freq_B, scientific = FALSE, nsmall = 4)
merged_side_by_side$count_C <- format(merged_side_by_side$count_C, scientific = FALSE)
merged_side_by_side$freq_C <- format(merged_side_by_side$freq_C, scientific = FALSE, nsmall = 4)
colnames(merged_side_by_side) <- c("HLA-A", "count", "freq", "HLA-B", "count", "freq", "HLA-C", "count", "freq")
#Export the merged data frame to a text file

write.table(merged_side_by_side, file = "YYYYMMDD_output_files/009_HLA_population_summary.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

rm(HLAA, HLAB, HLAC, max_rows, HLAA_padded, HLAB_padded, HLAC_padded, merged_side_by_side, homoA, homoB, homoC, unique_filtered_rows, number_of_homozygous_samples, summarize_allele_frequencies, filtered_data)   
