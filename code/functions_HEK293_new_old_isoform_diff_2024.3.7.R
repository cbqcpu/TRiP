library(readr)
library(dplyr)
library(stringr)

process_normalize_counts <- function(renamelist_path, spikein_counts_path, counts_path) {
  # Load the spike-in count data
  spikein_data <- read_csv(spikein_counts_path)
  # Calculate the spike-in totals
  spikein_totals <- colSums(spikein_data[,-1])
  
  # Read the CSV file to get file names
  file_names_df <- read_csv(renamelist_path)
  
  # Initialize an empty data frame to collect all samples data
  all_samples_df <- NULL
  
  # Iterate over each file name and process
  for (file_name in file_names_df$new_name) {
    # Construct paths to the count and attribute tables
    isoforms_counts_path <- file.path(counts_path, file_name, "isoforms.count_table")
    isoforms_attr_path <- file.path(counts_path, file_name, "isoforms.attr_table")
    samples_path <- file.path(counts_path, file_name, "samples.table")
    
    # Read the data from the paths
    isoforms_counts <- read_delim(isoforms_counts_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    isoforms_attr <- read_delim(isoforms_attr_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    samples_norm <- read_delim(samples_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
    # Add gene information and adjust counts by internal scale
    isoforms_counts$gene_id <- isoforms_attr$gene_id
    isoforms_counts$gene_name <- isoforms_attr$gene_short_name
    isoforms_counts$length <- isoforms_attr$length
    isoforms_counts[,2] <- round(isoforms_counts[,2] * samples_norm$internal_scale[1])
    isoforms_counts[,3] <- round(isoforms_counts[,3] * samples_norm$internal_scale[2])
    
    # Merge with all_samples_df
    if (is.null(all_samples_df)) {
      all_samples_df <- isoforms_counts
    } else {
      all_samples_df <- full_join(all_samples_df, isoforms_counts, 
                                  by = c("tracking_id", "gene_id", "gene_name", "length"))
    }
  }
  
  # Normalize the counts using the spike-in counts
  sample_names <- unique(gsub("_(new|old)_0", "", names(spikein_totals)))
  for (sample_name in sample_names) {
    new_col <- paste0(sample_name, "_new_0")
    old_col <- paste0(sample_name, "_old_0")
    spike_in_count <- spikein_totals[sample_name]
    if (!is.na(spike_in_count) && spike_in_count > 0) {
      all_samples_df[[new_col]] <- all_samples_df[[new_col]] / spike_in_count * 1000000
      all_samples_df[[old_col]] <- all_samples_df[[old_col]] / spike_in_count * 1000000
    } else {
      warning(paste("No corresponding spike-in total found or is zero for sample", sample_name))
    }
  }
  
  all_samples_df <- all_samples_df %>%
    dplyr::select(tracking_id, gene_id, gene_name, length, everything())
  # Return the processed dataframe
  return(all_samples_df)
}

# Example usage:
# all_samples_df_30min <- process_normalize_counts('/path/to/renamelist_NC30min.csv', '/path/to/NC30min_yeast_all_counts.csv', './mapping_NC30min/cuffnorm_counts')
# all_samples_df_1h <- process_normalize_counts('/path/to/renamelist_NC1h.csv', '/path/to/NC1h_yeast_all_counts.csv', './mapping_NC1h/cuffnorm_counts')
# ... and so on for other samples


merge_reps <- function(data, sample_name, replicates) {
  # Create an empty data frame to store merged results
  merged_data <- data.frame(tracking_id = data$tracking_id)
  
  # Iterate over each fraction and category to calculate the mean for each
  for (fraction in 1:4) {
    # Initialize columns for each category and fraction
    new_cols <- c()
    old_cols <- c()
    
    # Identify columns corresponding to the current fraction for each replicate, for 'new' and 'old'
    for (replicate in replicates) {
      new_cols <- c(new_cols, grep(paste0(sample_name, "-", replicate, "-", fraction, "_new_0"), colnames(data), value = TRUE))
      old_cols <- c(old_cols, grep(paste0(sample_name, "-", replicate, "-", fraction, "_old_0"), colnames(data), value = TRUE))
    }
    
    # Calculate the mean across the columns identified for 'new' and 'old' separately
    merged_data[[paste0(sample_name, "_fraction", fraction, "_new_avg")]] <- rowMeans(data[, new_cols], na.rm = TRUE)
    merged_data[[paste0(sample_name, "_fraction", fraction, "_old_avg")]] <- rowMeans(data[, old_cols], na.rm = TRUE)
  }
  
  return(merged_data)
}


calculate_ratios <- function(df, sample_name) {
  # Identify columns for new and old counts for each fraction
  fractions <- 1:4
  ratio_columns <- c()
  
  for (fraction in fractions) {
    new_col_name <- paste0("_fraction", fraction, "_new_avg")
    old_col_name <- paste0("_fraction", fraction, "_old_avg")
    
    new_col <- grep(new_col_name, names(df), value = TRUE)
    old_col <- grep(old_col_name, names(df), value = TRUE)
    
    # Ensure numerical computation is correct: new / (new + old)
    # This ratio should always yield a number between 0 and 1 if computed correctly.
    ratio_name <- paste0(sample_name, "_fraction", fraction, "_ratio")
    df[[ratio_name]] <- with(df, df[[new_col]] / (df[[new_col]] + df[[old_col]]))
    
    # Replace any potential NaN values (which could occur if both new and old are 0) with 0
    df[[ratio_name]][is.na(df[[ratio_name]])] <- 0
    
    ratio_columns <- c(ratio_columns, ratio_name)
  }
  
  # Create the output dataframe containing the tracking_id and the calculated ratios
  output_df <- df[c("tracking_id", ratio_columns)]
  
  return(output_df)
}







