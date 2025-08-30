library(readr)
library(stringr)
library(tidyverse)
library(cowplot)
library(corrplot)
library(matrixStats)
read_counts <- function(dir_path){
  # Get file names from the directory
  file_names <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE)
  counts <- list()
  # Loop through each file and read it into a data frame
  for (f in file_names) {
    # Read the TSV file into a data frame
    df <- read_delim(f, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 2)
    # Add the data frame to the list
    files <- str_split(f, "/")
    fi <- files[[1]][length(files[[1]])]
    sample_name <- str_split(fi, "\\.")[[1]][1]
    print(sample_name)
    counts[[sample_name]] <- df
  }
  return(counts)
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

compare_reps_frac <- function(data, sample_name, reps) {
  # Get the data for all fractions for each replicate
  data_frac <- list()
  for (i in 1:4) {
    temp_data_frac <- list()
    for (j in reps) {
      temp_data_frac[[j]] <- data[paste(sample_name, j, i, sep = "-")][[1]]$ReadCount
    }
    data_frac[[i]] <- do.call(cbind, temp_data_frac)
  }

  # Calculate the pairwise correlations of each fraction across replicates
  corr_mat <- lapply(data_frac, function(x) cor(x))
  corr_mat <- lapply(corr_mat, function(x) {rownames(x) <- colnames(x) <- c("rep. 1", "rep. 2", "rep. 3"); return(x)})
  # Plot the correlation matrices in a grid with fraction names as plot titles
  par(mfrow = c(1, 4), mar = c(2, 2, 1, 1))
  for (i in 1:4) {
    frac_name <- paste0("Fraction ", i)
    corrplot(corr_mat[[i]], method = 'number', tl.col = "black", main = frac_name, col = COL2('PiYG'), mar=c(0,0,5,0), cl.pos="n")
  }
}

merge_reps_with_cv <- function(data, condition, time, reps, spike_in_counts){
  n_reps <- length(reps)
  
  # separate fractions
  data_frac <- list()
  data_frac_TC <- list()
  data_frac_new <- list()
  spike_in_frac <- list()
  TC_rate_merge_sum <- 0
  
  for (i in 1:4) {
    temp_data_frac <- list()
    temp_data_frac_TC <- list()
    temp_data_frac_new <- list()
    temp_spike_in_frac <- list()
    
    for (j in reps) {
      temp_data_frac[[j]] <- data[paste(condition, j, i, sep = "-")][[1]]$ReadCount
      temp_data_frac_TC[[j]] <- data[paste(condition, j, i, sep = "-")][[1]]$ConversionRate
      temp_data_frac_new[[j]] <- data[paste(condition, j, i, sep = "-")][[1]]$TcReadCount
      temp_spike_in_frac[[j]] <- sum(spike_in_counts[[paste0(condition, "-", j, "-", i)]])
    }
    
    data_frac[[i]] <- do.call(cbind, temp_data_frac)
    data_frac_TC[[i]] <- do.call(cbind, temp_data_frac_TC)
    data_frac_new[[i]] <- do.call(cbind, temp_data_frac_new)
    spike_in_frac[[i]] <- unlist(temp_spike_in_frac)
  }
  
  # Normalize the raw counts by the summed spike-in counts for each replicate
  norm_counts_cv <- lapply(seq_along(data_frac), function(i) data_frac[[i]] / spike_in_frac[[i]])
  
  # Compute the mean and standard deviation for CV
  mean_counts <- lapply(norm_counts_cv, rowMeans)
  sd_counts <- lapply(norm_counts_cv, matrixStats::rowSds)
  cv_counts <- mapply(function(mean, sd) sd / mean, mean_counts, sd_counts, SIMPLIFY = FALSE)
  # Merge the replicates in the spike-in counts
  spike_in_counts_merge <- data.frame(
    frac1 = rowSums(spike_in_counts[, paste0(condition, "-", reps, "-1")]),
    frac2 = rowSums(spike_in_counts[, paste0(condition, "-", reps, "-2")]),
    frac3 = rowSums(spike_in_counts[, paste0(condition, "-", reps, "-3")]),
    frac4 = rowSums(spike_in_counts[, paste0(condition, "-", reps, "-4")])
  )
  # Merge the replicates by summing up the counts for each gene
  raw_counts <- lapply(data_frac, rowSums)
  raw_counts_new <- lapply(data_frac_new, rowSums)
  spike_in_counts_merge <- colSums(spike_in_counts_merge)
  norm_counts <- mapply(function(x, y) x / y, raw_counts, spike_in_counts_merge, SIMPLIFY = FALSE)
  gene_length <- data[paste(condition, reps[1], "1", sep = "-")][[1]]$Length
  #data_tpm <- tpm(do.call(cbind, norm_counts), gene_length) TPM is not necessary
  data_norm_counts <- do.call(cbind, norm_counts)*sum(spike_in_counts_merge)/4
  colnames(data_norm_counts) <- paste0("frac", 1:4)
  
  # Merge the replicates by mean T to C mutation rate
  TC_rate_merge <- lapply(data_frac_TC, function(x) rowSums(x)/n_reps)
  gene_id <- data[paste(condition, reps[1], "1", sep = "-")][[1]][,1:4]
  TC_rate_merge <- do.call(cbind, TC_rate_merge)
  TC_rate_merge_sum <- TC_rate_merge_sum + sum(TC_rate_merge)
  TC_rate_merge <- TC_rate_merge/sum(TC_rate_merge)
  colnames(TC_rate_merge) <- paste0("frac", 1:4)
  
  raw_counts_new_norm <- mapply(function(x, y) x / y, raw_counts_new, spike_in_counts_merge, SIMPLIFY = FALSE)
  data_raw_counts_new_norm <- do.call(cbind, raw_counts_new_norm)*sum(spike_in_counts_merge)/4
  colnames(data_raw_counts_new_norm) <- paste0("frac", 1:4)
  new_percentage <- data_raw_counts_new_norm/data_norm_counts
  
  return(list(gene_id = gene_id, TC_rate_merge = TC_rate_merge*TC_rate_merge_sum, data_norm_counts = data_norm_counts, data_raw_counts_new_norm = data_raw_counts_new_norm, new_percentage = new_percentage, cv = cv_counts))
}


filter_high_cv <- function(dataset, high_cv_threshold) {
  # Initialize a list to store the genes_to_keep for each fraction
  all_genes_to_keep <- list()
  
  for(i in 1:4) {
    # Calculate genes_to_keep for each fraction
    genes_to_keep <- which(dataset$cv[[i]] <= high_cv_threshold)
    all_genes_to_keep[[i]] <- genes_to_keep
  }
  
  # Find the intersection of all_genes_to_keep across the four fractions
  intersection_genes <- Reduce(intersect, all_genes_to_keep)
  
  return(intersection_genes)
}

filter_dataset_by_genes <- function(dataset, common_genes) {
  # Filter gene_id and convert to a data frame
  new_gene_id <- data.frame(lapply(dataset$gene_id, function(gene_id) gene_id[common_genes]))
  # Filter the lists in the dataset using common_genes
  dataset$gene_id <- new_gene_id
  dataset$data_norm_counts <- dataset$data_norm_counts[common_genes, ]
  dataset$data_raw_counts_new_norm <- dataset$data_raw_counts_new_norm[common_genes, ]
  dataset$new_percentage <- dataset$new_percentage[common_genes, ]
  dataset$TC_rate_merge <- dataset$TC_rate_merge[common_genes, ]
  dataset$cv <- lapply(dataset$cv, function(cv) cv[common_genes])
  
  return(dataset)
}


plot_violin <- function(merge_df, time) {
  # subset data to the given time
  time_cols <- grep(paste0("^", time, "_frac"), colnames(merge_df))
  time_data <- merge_df[, time_cols]
  
  # filter out rows where all values are zero
  time_data_filtered <- time_data[rowSums(time_data) > 0, ]
  
  # reshape data for plotting
  time_data_long <- time_data_filtered %>%
    gather(key = "fraction", value = value) %>%
    mutate(fraction = gsub(paste0("^", time, "_frac"), "", fraction)) %>%
    mutate(fraction = factor(fraction, levels = c("1", "2", "3", "4")))
  
  # create violin plot
  p <- ggplot(time_data_long, aes(x = fraction, y = log10(value))) +
    geom_violin(aes(fill = fraction), alpha = 0.6) +
    labs(x = "Fractions", y = "T to C rate (log 10)", 
         title = paste0(time, " 4sU convertion rate"), 
         fill = "Fractions") +
    scale_x_discrete(labels = c("Frac. 1", "Frac. 2", "Frac. 3", "Frac. 4")) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "white") +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "#2495f1", linewidth = 1) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  p <- p + scale_fill_brewer(palette = "Set1")
  return(p)
}

plot_violin_new_ratio <- function(merge_df, time) {
  # subset data to the given time
  time_cols <- grep(paste0("^", time, "_frac"), colnames(merge_df))
  time_data <- merge_df[, time_cols]
  
  # filter out rows where all values are zero
  time_data_filtered <- time_data[rowSums(time_data) > 0, ]
  
  # reshape data for plotting
  time_data_long <- time_data_filtered %>%
    gather(key = "fraction", value = value) %>%
    mutate(fraction = gsub(paste0("^", time, "_frac"), "", fraction)) %>%
    mutate(fraction = factor(fraction, levels = c("1", "2", "3", "4")))
  
  # create violin plot
  p <- ggplot(time_data_long, aes(x = fraction, y = log10(value))) +
    geom_violin(aes(fill = fraction), alpha = 0.6) +
    labs(x = "Fractions", y = "new mRNA ratio (log 10)", 
         title = paste0(time, " new mRNA ratio"), 
         fill = "Fractions") +
    scale_x_discrete(labels = c("Frac. 1", "Frac. 2", "Frac. 3", "Frac. 4")) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", fill = "white") +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "#2495f1", linewidth = 1) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 14, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  p <- p + scale_fill_brewer(palette = "Set1")
  return(p)
}

calculate_polysome_density <- function(merge_df, time) {
  # initialize a list to store the results for each time point
  density_list <- list()

  for (t in time) {
    # subset data to the given time
    time_cols <- grep(paste0("^", t, "_frac"), colnames(merge_df))
    time_data <- merge_df[, time_cols]

    # filter out rows where all values are zero
    #time_data_filtered <- time_data[rowSums(time_data) > 0, ]
    time_data_filtered <- time_data
    normalize <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)))
    }
    
    # Apply the function to each row
    #df_normalized <- as.data.frame(t(apply(time_data_filtered, 1, normalize)))
    df_normalized <- time_data_filtered
    # calculate polysome density
    time_data_density <- data.frame(rowSums(t(apply(df_normalized, 1,function(row) row * c(0.01, 1, 4, 9)))))
    colnames(time_data_density) <- t

    # add the result to the list
    density_list[[t]] <- time_data_density
  }

  # combine the results into a single data frame
  density_df <- do.call(cbind, density_list)
  #density_df <- data.frame(gene_id = merge_df$gene_id, density_df)
  return(density_df)
}

calculate_polysome_load <- function(new_df, old_df, times) {
  # Initialize data frames with the correct number of rows
  n_rows <- nrow(new_df) # Assuming new_df and old_df have the same number of rows
  polysome_load_new_df <- data.frame(matrix(ncol = length(times), nrow = n_rows))
  polysome_load_old_df <- data.frame(matrix(ncol = length(times), nrow = n_rows))
  
  # Normalize function
  normalize <- function(x) {
    return(x / sum(x))
  }
  
  # Polysome load calculation function
  polysome_load_calc <- function(row) {
    return(row * c(0.01, 1, 4, 9))
  }
  
  # Loop through each time point
  for (i in seq_along(times)) {
    t <- times[i]
    
    # Select columns based on the time parameter
    new_time_cols <- grep(paste0("^", t, "_"), colnames(new_df), value = TRUE)
    old_time_cols <- grep(paste0("^", t, "_"), colnames(old_df), value = TRUE)
    
    new_time_df <- new_df[, new_time_cols]
    old_time_df <- old_df[, old_time_cols]
    
    # Rename columns to add prefixes
    colnames(new_time_df) <- paste0("new_", colnames(new_time_df))
    colnames(old_time_df) <- paste0("old_", colnames(old_time_df))
    
    # Combine the dataframes
    combined_df <- cbind(new_time_df, old_time_df)
    
    # Normalize the combined dataframe
    combined_df_normalized <- as.data.frame(t(apply(combined_df, 1, normalize)))
    
    # Split the normalized data back into new and old dataframes
    new_cols <- grep("^new_", colnames(combined_df_normalized), value = TRUE)
    old_cols <- grep("^old_", colnames(combined_df_normalized), value = TRUE)
    
    new_df_normalized <- combined_df_normalized[, new_cols]
    old_df_normalized <- combined_df_normalized[, old_cols]
    
    # Calculate polysome load for new and old dataframes
    polysome_load_new <- apply(new_df_normalized, 1, polysome_load_calc)
    polysome_load_old <- apply(old_df_normalized, 1, polysome_load_calc)
    
    # Store the results
    polysome_load_new_df[, i] <- rowSums(t(polysome_load_new))
    polysome_load_old_df[, i] <- rowSums(t(polysome_load_old))
    colnames(polysome_load_new_df)[i] <- paste0("", t)
    colnames(polysome_load_old_df)[i] <- paste0("", t)
  }
  
  # Return a list containing the polysome load for new and old dataframes
  return(list(new_polysome_load = polysome_load_new_df, old_polysome_load = polysome_load_old_df))
}


calculate_polysome_load_old <- function(merge_df, time) {
  # initialize a list to store the results for each time point
  density_list <- list()
  
  for (t in time) {
    # subset data to the given time
    time_cols <- grep(paste0("^", t, "_frac"), colnames(merge_df))
    time_data <- merge_df[, time_cols]
    
    # filter out rows where all values are zero
    #time_data_filtered <- time_data[rowSums(time_data) > 0, ]
    time_data_filtered <- time_data
    normalize <- function(x) {
      return (x / sum(x))
    }
    
    # Apply the function to each row
    df_normalized <- as.data.frame(t(apply(time_data_filtered, 1, normalize)))
    # calculate polysome density
    polysome_load <- data.frame(rowSums(t(apply(df_normalized, 1,function(row) row * c(0.01, 1, 4, 9)))))
    colnames(polysome_load) <- t
    
    # add the result to the list
    density_list[[t]] <- polysome_load
  }
  
  # combine the results into a single data frame
  density_df <- do.call(cbind, density_list)
  #density_df <- data.frame(gene_id = merge_df$gene_id, density_df)
  return(density_df)
}



specific_gene_plot_fraction <- function(data_frame, geneid) {
  # Specific gene:
  df_long <- data_frame[data_frame$gene_id==geneid,] %>%
    pivot_longer(cols = -gene_id, 
                 names_to = c("time", "fraction"), 
                 names_pattern = "NC_(.*)_frac(\\d)")
  
  # Convert time and fraction to factors for ordering in the plot
  df_long$time <- factor(df_long$time, levels = c("30min", "1h", "2h"))
  df_long$fraction <- factor(df_long$fraction, levels = c("1", "2", "3", "4"))
  
  # Normalize the value by subtracting the mean and dividing by the standard deviation
  df_long <- df_long %>%
    group_by(time, fraction) %>%
    mutate(value = (value - min(value)) / (max(value)-min(value)))
  # Define custom colors
  custom_colors <- c("30min" = "#AED6F1", "1h" = "#5DADE2", "2h" = "#2874A6")
  # Create the plot
  ggplot(df_long, aes(x = fraction, y = value, group = time, color = time)) +
    geom_line(size = 1) +
    scale_color_manual(values = custom_colors) +
    theme_classic()
}


get_transcript_components <- function(ensembl_mart, transcript_list) {
  
  results_df <- data.frame(transcript_id = character(), 
                           UTR5 = character(), 
                           coding = character(),
                           UTR3 = character(), 
                           stringsAsFactors = FALSE)
  
  batch_size <- 50
  batches <- split(transcript_list, ceiling(seq_along(transcript_list)/batch_size))
  
  for (batch in batches) {
    
    utr5_sequences <- getSequence(id = batch, 
                                  type = "ensembl_transcript_id", 
                                  seqType = "5utr", 
                                  mart = ensembl_mart)
    
    exon_intron_sequences <- getSequence(id = batch,
                                  type = "ensembl_transcript_id", 
                                  seqType = "coding",
                                  mart = ensembl_mart)
    
    utr3_sequences <- getSequence(id = batch, 
                                  type = "ensembl_transcript_id", 
                                  seqType = "3utr", 
                                  mart = ensembl_mart)
    
    for (transcript_id in batch) {
      results_df <- results_df %>% add_row(
        transcript_id = transcript_id, 
        UTR5 = ifelse(transcript_id %in% utr5_sequences$ensembl_transcript_id, utr5_sequences$`5utr`[utr5_sequences$ensembl_transcript_id == transcript_id], NA), 
        coding = ifelse(transcript_id %in% exon_intron_sequences$ensembl_transcript_id, exon_intron_sequences$coding[exon_intron_sequences$ensembl_transcript_id == transcript_id], NA), 
        UTR3 = ifelse(transcript_id %in% utr3_sequences$ensembl_transcript_id, utr3_sequences$`3utr`[utr3_sequences$ensembl_transcript_id == transcript_id], NA)
      )
    }
  }
  
  return(results_df)
}


get_gene_sequences <- function(ensembl_mart, gene_list) {
  
  results_df <- data.frame(gene_id = character(), 
                           gene_sequence = character(), 
                           stringsAsFactors = FALSE)
  
  batch_size <- 50
  batches <- split(gene_list, ceiling(seq_along(gene_list)/batch_size))
  
  for (batch in batches) {
    
    gene_sequences <- getSequence(id = batch, 
                                  type = "ensembl_gene_id", 
                                  seqType = "gene_exon_intron", # or choose the relevant seqType
                                  mart = ensembl_mart)
    
    for (gene_id in batch) {
      results_df <- results_df %>% add_row(
        gene_id = gene_id, 
        gene_sequence = ifelse(gene_id %in% gene_sequences$ensembl_gene_id, gene_sequences$gene_exon_intron[gene_sequences$ensembl_gene_id == gene_id], NA)
      )
    }
  }
  
  return(results_df)
}



