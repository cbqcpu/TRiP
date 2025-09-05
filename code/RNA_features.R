library(readxl)
library(dplyr)
library(ggplot2)

compare_boxplot <- function(df1, df2, feature, log = FALSE, name1 = "df1", name2 = "df2") {
  # Remove missing values from data frames
  df1 <- na.omit(df1[,c("gene_id",feature)])
  df2 <- na.omit(df2[,c("gene_id",feature)])
  
  # Concatenate the data frames vertically
  combined_df <- rbind(df1, df2)
  
  # Create a grouping variable to distinguish between the two data frames
  combined_df$Group <- c(rep(name1, nrow(df1)), rep(name2, nrow(df2)))
  
  # Apply logarithm if log is TRUE
  if (log) {
    combined_df[[feature]] <- log10(combined_df[[feature]])
    combined_df <- combined_df[is.finite(combined_df[[feature]]), ]
  }
  
  # Calculate t-test
  t_test_result <- t.test(combined_df[combined_df$Group == name1, feature],
                          combined_df[combined_df$Group == name2, feature])
  
  # Create plot
  p <- ggplot(combined_df, aes(x = Group, y = .data[[feature]], fill = Group)) +
    geom_boxplot() +
    labs(x = paste("Number of genes\n", name1, ":", nrow(df1), ",", name2, ":", nrow(df2)), y = feature, title = feature) +
    scale_fill_manual(values = c("orange", "steelblue"),
                      labels = c(name1, name2)) +
    theme_classic()
  
  # Add p-value to the plot
  p <- p + annotate("text", x=1.5, y=max(combined_df[[feature]]) * 1.1, label=paste("p =", format.pval(t_test_result$p.value)), color = "black")
  
  # Adjust plot scale
  p <- p + coord_cartesian(ylim = c(min(combined_df[[feature]]), max(combined_df[[feature]]) * 1.2))
  
  return(p)
}



plot_m6A_levels <- function(df_list, compress_ratios) {
  colors <- c("orange", "steelblue")  # Adjust this to the number of dataframes you might input
  
  compress_vector <- function(vector, compress_ratio) {
    bins <- cut(seq_along(vector), breaks = compress_ratio, labels = FALSE)
    sapply(split(vector, bins), mean, na.rm = TRUE)
  }
  
  normalize <- function(x, min_val, max_val) {
    return ((x - min_val) / (max_val - min_val))
  }
  
  data_to_plot <- data.frame()
  all_data <- list()
  
  for (i in 1:length(df_list)) {
    df <- df_list[[i]]
    # Convert string to numeric vector
    df$m6A_level_5utr <- lapply(df$m6A_level_5utr, function(x) as.numeric(strsplit(x, ",")[[1]]))
    df$m6A_level_3utr <- lapply(df$m6A_level_3utr, function(x) as.numeric(strsplit(x, ",")[[1]]))
    df$m6A_level_cds <- lapply(df$m6A_level_cds, function(x) as.numeric(strsplit(x, ",")[[1]]))
    
    # Summing the vectors
    sum_5utr <- Reduce('+', df$m6A_level_5utr)
    sum_3utr <- Reduce('+', df$m6A_level_3utr)
    sum_cds <- Reduce('+', df$m6A_level_cds)
    
    # Normalizing the sum by the number of rows in the dataframe
    sum_5utr <- sum_5utr / nrow(df)
    sum_3utr <- sum_3utr / nrow(df)
    sum_cds <- sum_cds / nrow(df)
    
    # Compressing the vectors
    sum_5utr <- compress_vector(sum_5utr, compress_ratios[1])
    sum_cds <- compress_vector(sum_cds, compress_ratios[2])
    sum_3utr <- compress_vector(sum_3utr, compress_ratios[3])
    
    # Concatenating the vectors
    concat_vector <- c(sum_5utr, sum_cds, sum_3utr)
    
    # Storing the sums
    all_data[[i]] <- concat_vector
  }
  
  # Find the overall min and max
  overall_min <- min(unlist(all_data))
  overall_max <- max(unlist(all_data))
  
  # Normalize using the overall min and max
  all_data <- lapply(all_data, normalize, min_val = overall_min, max_val = overall_max)
  
  for (i in 1:length(df_list)) {
    concat_vector <- all_data[[i]]
    
    # Create the x axis values (1 to length of the concatenated vector)
    x_values <- 1:length(concat_vector)
    
    # Create a dataframe for plotting
    df_to_plot <- data.frame(Number = x_values, Value = concat_vector, DataSet = rep(names(df_list)[i], length(concat_vector)))
    
    # Bind the data
    data_to_plot <- rbind(data_to_plot, df_to_plot)
  }
  
  # Plotting with ggplot
  ggplot(data_to_plot, aes(x = Number, y = Value, color = DataSet)) +
    geom_line(size = 2, alpha = 0.5) +
    scale_color_manual(values = colors) +
    theme_classic() +
    labs(x = "Number", y = "Value") +
    theme(legend.title = element_blank())
}

plot_psi_levels <- function(df_list, compress_ratios) {
  colors <- c("orange", "steelblue")  # Adjust this to the number of dataframes you might input
  
  compress_vector <- function(vector, compress_ratio) {
    bins <- cut(seq_along(vector), breaks = compress_ratio, labels = FALSE)
    sapply(split(vector, bins), mean, na.rm = TRUE)
  }
  
  normalize <- function(x, min_val, max_val) {
    return ((x - min_val) / (max_val - min_val))
  }
  
  data_to_plot <- data.frame()
  all_data <- list()
  
  for (i in 1:length(df_list)) {
    df <- df_list[[i]]
    # Convert string to numeric vector
    df$psi_level_5utr <- lapply(df$psi_level_5utr, function(x) as.numeric(strsplit(x, ",")[[1]]))
    df$psi_level_3utr <- lapply(df$psi_level_3utr, function(x) as.numeric(strsplit(x, ",")[[1]]))
    df$psi_level_cds <- lapply(df$psi_level_cds, function(x) as.numeric(strsplit(x, ",")[[1]]))
    
    # Summing the vectors
    sum_5utr <- Reduce('+', df$psi_level_5utr)
    sum_3utr <- Reduce('+', df$psi_level_3utr)
    sum_cds <- Reduce('+', df$psi_level_cds)
    
    # Normalizing the sum by the number of rows in the dataframe
    sum_5utr <- sum_5utr / nrow(df)
    sum_3utr <- sum_3utr / nrow(df)
    sum_cds <- sum_cds / nrow(df)
    
    # Compressing the vectors
    sum_5utr <- compress_vector(sum_5utr, compress_ratios[1])
    sum_cds <- compress_vector(sum_cds, compress_ratios[2])
    sum_3utr <- compress_vector(sum_3utr, compress_ratios[3])
    
    # Concatenating the vectors
    concat_vector <- c(sum_5utr, sum_cds, sum_3utr)
    
    # Storing the sums
    all_data[[i]] <- concat_vector
  }
  
  # Find the overall min and max
  overall_min <- min(unlist(all_data))
  overall_max <- max(unlist(all_data))
  
  # Normalize using the overall min and max
  all_data <- lapply(all_data, normalize, min_val = overall_min, max_val = overall_max)
  
  for (i in 1:length(df_list)) {
    concat_vector <- all_data[[i]]
    
    # Create the x axis values (1 to length of the concatenated vector)
    x_values <- 1:length(concat_vector)
    
    # Create a dataframe for plotting
    df_to_plot <- data.frame(Number = x_values, Value = concat_vector, DataSet = rep(names(df_list)[i], length(concat_vector)))
    
    # Bind the data
    data_to_plot <- rbind(data_to_plot, df_to_plot)
  }
  
  # Plotting with ggplot
  ggplot(data_to_plot, aes(x = Number, y = Value, color = DataSet)) +
    geom_line(size = 2, alpha = 0.5) +
    scale_color_manual(values = colors) +
    theme_classic() +
    labs(x = "Number", y = "Value") +
    theme(legend.title = element_blank())
}

cal_polyA_length <- function(TED_seq_data) {
#This function use TED-Seq data to map Ensembl transcript ID.
  # Load the necessary library
  library(biomaRt)
  
  # Select the human database
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  polyA_length_TED_Seq <- TED_seq_data
  # Define RefSeq IDs
  refSeqIDs <- polyA_length_TED_Seq$refFlatID
  
  # Use the getBM function to get the corresponding Ensembl IDs
  results <- getBM(
    attributes = c('refseq_mrna', 'ensembl_gene_id'), 
    filters = 'refseq_mrna',
    values = refSeqIDs, 
    mart = mart
  )
  # Merge the results dataframe with your original dataframe
  merged_df <- merge(polyA_length_TED_Seq, results, by.x = "refFlatID", by.y = "refseq_mrna", all.x = TRUE)
  merged_df$ensembl_gene_id[is.na(merged_df$ensembl_gene_id)] <- 'Not available'
  merged_df_new <- data.frame(gene_id = merged_df$ensembl_gene_id, polya_length = merged_df$TED1.DMSO)
  
  # Filter out rows where the gene_id is 'Not available'
  filtered_df <- merged_df_new %>%
    filter(gene_id != 'Not available')
  
  # Convert the polya_length column to numeric
  filtered_df$polya_length <- as.numeric(filtered_df$polya_length)
  
  # Group by gene_id and calculate the average polya_length for each group
  final_df <- filtered_df %>%
    group_by(gene_id) %>%
    summarise(polya_length = mean(polya_length, na.rm = TRUE))
  
  # If you want this as a regular dataframe, rather than a tibble, you can convert it with:
  final_df <- as.data.frame(final_df)
  return(final_df)
}


compare_boxplot4 <- function(df1, df2, df3, df4, feature, log = FALSE, name1 = "df1", name2 = "df2", name3 = "df3", name4 = "df4") {
  
  # Remove missing values from data frames
  df1 <- aggregate(. ~ gene_id, na.omit(df1[,c("gene_id",feature)]), mean)
  df2 <- aggregate(. ~ gene_id, na.omit(df2[,c("gene_id",feature)]), mean)
  df3 <- aggregate(. ~ gene_id, na.omit(df3[,c("gene_id",feature)]), mean)
  df4 <- aggregate(. ~ gene_id, na.omit(df4[,c("gene_id",feature)]), mean)
  
  # Concatenate the data frames vertically
  combined_df <- rbind(df1, df2, df3, df4)
  
  # Create a grouping variable to distinguish between the data frames
  combined_df$Group <- c(rep(name1, nrow(df1)), rep(name2, nrow(df2)), rep(name3, nrow(df3)), rep(name4, nrow(df4)))
  
  # Apply logarithm if log is TRUE
  if (log) {
    combined_df[[feature]] <- log10(combined_df[[feature]])
    combined_df <- combined_df[is.finite(combined_df[[feature]]), ]
  }
  
  # Create plot
  p <- ggplot(combined_df, aes(x = Group, y = .data[[feature]], fill = Group)) +
    geom_boxplot() +
    labs(x = paste("Number of genes\n", name1, ":", nrow(df1), ",", name2, ":", nrow(df2), ",", name3, ":", nrow(df3), ",", name4, ":", nrow(df4)), y = feature, title = feature) +
    scale_fill_manual(values = c("#FF7F50", "#DE3163", "#6495ED", "#40E0D0"),
                      labels = c(name1, name2, name3, name4)) +
    theme_classic()
  
  # Add p-value to the plot - removed as it's not applicable for more than 2 groups
  # p <- p + annotate("text", x=2.5, y=max(combined_df[[feature]]) * 1.1, label=paste("p =", format.pval(t_test_result$p.value)), color = "black")
  
  # Adjust plot scale
  p <- p + coord_cartesian(ylim = c(min(combined_df[[feature]]), max(combined_df[[feature]]) * 1.2))
  
  return(p)
}

