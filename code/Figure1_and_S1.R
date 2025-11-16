#Figure1
library(readr)
library(Mfuzz)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)

source("code/polysome_seq_functions.R")
load("data/HEK293T_preprocessed.RData")
load("data/RNA_features_gene_level_20240418.RData")

base_theme <- theme_bw() + theme(
  legend.position = "none",  # Remove legend
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
  panel.grid.major = element_blank(),  # Remove major grid lines
  panel.grid.minor = element_blank(),  # Remove minor grid lines
  panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
)

merge_df_counts_select <- merge_df_counts[rowMeans(merge_df_counts[,-1])>1000,]

merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts$gene_id %in% merge_df_counts_select$gene_id,]

merge_df_counts_old_select <- data.frame(gene_id = merge_df_counts_select$gene_id, merge_df_counts_select[,-1] - merge_df_counts_new_select[,-1])

MRLs <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))
MRLs2 <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_select, c("NC_30min", "NC_1h", "NC_2h"))

MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]
MRL_total <- MRLs2[2]

MRL_diff <- MRL_new-MRL_old
data_combined <- data.frame(
  NC_2h = c(MRL_old$NC_2h, MRL_new$NC_2h),
  Group = rep(c("MRL_old", "MRL_new"), 
              times = c(length(MRL_old$NC_2h), length(MRL_new$NC_2h)))
)

Fig1C_hist_MRL <- ggplot(data_combined, aes(x = NC_2h, fill = Group)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, position = "identity", bins = 30) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("#89a6d6", "#c0e4ef")) +  # User-friendly colors
  labs(title = "Density Histogram of New and Old mRNA MRL in 2 h",
       x = "MRL",
       y = "Density") +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )
#write.csv(data_combined, file = "Fig1C.csv", quote = F, row.names = F)
data <- data.frame(
  NC_2h_old = MRL_old$NC_2h,
  NC_2h_new = MRL_new$NC_2h
)


set.seed(1234)  # You can use any integer as the seed

MRL_df <- cbind(MRL_new, MRL_old)
row_max <- apply(MRL_df, 1, max)
MRL_new_normalized <- MRL_new / row_max
MRL_old_normalized <- MRL_old / row_max
colnames(MRL_new_normalized) <- c("new_30_min", "new_1h", "new_2h")
colnames(MRL_old_normalized) <- c("old_30_min", "old_1h", "old_2h")
data_matrix <- as.matrix(cbind(MRL_new_normalized, MRL_old_normalized))
rownames(data_matrix) <- merge_df_counts_select$gene_id

# Create ExpressionSet object from your data matrix
eset <- new('ExpressionSet', exprs=data_matrix)

# Determine the optimal fuzzification parameter m
m <- mestimate(eset)

# After visually inspecting the plot, choose the optimal number of clusters
optimal_c <- 2  # Set to 3 clusters

# Run mfuzz with the optimal parameters
cl <- mfuzz(eset, c=optimal_c, m=3)

# Extract membership and assign clusters
membership <- cl$membership
assigned_clusters <- apply(membership, 1, which.max)

# Sort the assigned_clusters and reorder the data accordingly
sorted_indices <- order(assigned_clusters)
sorted_data <- data_matrix[sorted_indices, ]
sorted_clusters <- assigned_clusters[sorted_indices]


# Define custom colors for clusters 1 and 2
cluster_colors <- c("1" = "#8887b9",  # Color for Cluster 1 (e.g., Tomato Red)
                    "2" = "#95c08b")  # Color for Cluster 2 (e.g., Steel Blue)

# Convert your clusters to a factor to make sure colors are assigned correctly
sorted_clusters <- as.factor(sorted_clusters)  # Ensures clusters 1 and 2 are treated as factors

# Update the row annotation based on the sorted clusters with custom colors
ha <- rowAnnotation(Cluster = sorted_clusters,
                    col = list(Cluster = cluster_colors),  # Assign custom colors to clusters
                    annotation_legend_param = list(Cluster = list(title = "Cluster")))

#save(sorted_indices, sorted_data, sorted_clusters, file = "data/figure1_heatmap_ordered_cluster.RData")

# Create the heatmap
Fig1D_phm <- Heatmap(sorted_data, 
                     name = "Normalized MRL",
                     left_annotation = ha,  # Add row annotation for cluster colors
                     col = colorRamp2(c(0, 0.5, 1),  
                                      c("white", "orange", "darkred")),  # Heatmap colors
                     show_row_names = FALSE,
                     cluster_columns = FALSE,  # Do not cluster columns
                     cluster_rows = FALSE)     # Do not cluster rows, use the sorted order

#write.csv(sorted_data, file = "Fig1D.csv", row.names = F, quote = F)
# Compare new and old mRNA loading to polysome dynamics with mRNA halflife









# test if not normalize by row, what the two cluster like

MRL_new_no_normalized <- MRL_new
MRL_old_no_normalized <- MRL_old
colnames(MRL_new_no_normalized) <- c("new_30_min", "new_1h", "new_2h")
colnames(MRL_old_no_normalized) <- c("old_30_min", "old_1h", "old_2h")
data_matrix_no <- as.matrix(cbind(MRL_new_no_normalized, MRL_old_no_normalized))
rownames(data_matrix_no) <- merge_df_counts_select$gene_id

data_range <- range(data_matrix_no, na.rm = TRUE)
phm_no <- Heatmap(data_matrix_no, 
                  name = "Normalized MRL",
                  left_annotation = ha,  # Row annotation for cluster colors
                  col = colorRamp2(c(data_range[1], 
                                     mean(data_range), 
                                     data_range[2]),  
                                   c("white", "orange", "darkred")),  # Heatmap colors
                  show_row_names = FALSE,
                  cluster_columns = FALSE,  # Do not cluster columns
                  cluster_rows = FALSE)     # Do not cluster rows

head(data_matrix_no)
head(sorted_clusters)
data_matrix_no <- as.data.frame(data_matrix_no) %>%
  rownames_to_column(var = "Gene")

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)  # For grid arrangement

# Convert the matrix to a data frame and add row names as a column

# Convert 'sorted_clusters' to a data frame for merging
cluster_assignments <- data.frame(Gene = names(sorted_clusters), Cluster = as.integer(sorted_clusters))

# Merge the data with cluster assignments
merged_data <- merge(data_matrix_no, cluster_assignments, by = "Gene")

# Reshape data to long format for plotting individual gene lines
long_data <- merged_data %>%
  pivot_longer(cols = starts_with("new_") | starts_with("old_"), names_to = c("Type", "Time"),
               names_pattern = "(new|old)_(.+)", values_to = "Expression") %>%
  mutate(Time = factor(Time, levels = c("30_min", "1h", "2h")))  # Set the order of Time points

# Calculate mean values for 'new' and 'old' mRNA across clusters and time points
cluster_means <- long_data %>%
  group_by(Cluster, Type, Time) %>%
  summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = 'drop')

# Function to create a plot for a specific cluster with 10 randomly selected genes
create_cluster_plot <- function(cluster_num, n_genes = 10) {
  # Filter data for the specified cluster and sample 10 genes
  cluster_data <- filter(long_data, Cluster == cluster_num)
  sampled_genes <- sample(unique(cluster_data$Gene), n_genes)
  sampled_data <- filter(cluster_data, Gene %in% sampled_genes)
  
  # Base ggplot object for the cluster
  p <- ggplot() +
    labs(title = paste("MRL for Cluster", cluster_num), x = "Time", y = "MRL") +
    scale_color_manual(values = c("new" = "dodgerblue", "old" = "orange")) +
    scale_x_discrete(expand = c(0, 0)) +  # Remove padding on x-axis to align time points with start/end
    theme_bw() +  # Use a minimal theme
    theme(
      legend.position = "none",  # Remove legend
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
    )
  
  # Add individual gene lines for the sampled genes with grouping by Gene and Type
  p <- p + geom_line(data = sampled_data, aes(x = Time, y = Expression, color = Type, group = interaction(Gene, Type)), 
                     size = 0.5, alpha = 0.1)
  
  # Add mean line and points for the specified cluster
  cluster_mean_data <- filter(cluster_means, Cluster == cluster_num)
  p <- p + geom_line(data = cluster_mean_data, aes(x = Time, y = MeanExpression, color = Type, group = Type), size = 1.2) +
    geom_point(data = cluster_mean_data, aes(x = Time, y = MeanExpression, color = Type), size = 2)
  
  return(p)
}



# Convert the matrix to a data frame and add row names as a column

# Convert 'sorted_clusters' to a data frame for merging
cluster_assignments <- data.frame(Gene = names(sorted_clusters), Cluster = as.integer(sorted_clusters))


# Reshape data to long format for plotting individual gene lines
# Combine Type and Time into a single variable for x-axis labeling
long_data <- merged_data %>%
  pivot_longer(cols = starts_with("new_") | starts_with("old_"), names_to = c("Type", "Time"),
               names_pattern = "(new|old)_(.+)", values_to = "Expression") %>%
  mutate(TimeType = factor(paste(Time, Type), 
                           levels = c("30_min new", "1h new", "2h new", "", "30_min old", "1h old", "2h old")),
         Cluster = factor(Cluster, levels = c(1, 2), labels = c("Cluster 1", "Cluster 2")))

# Calculate mean values for 'new' and 'old' mRNA across clusters and combined time points
cluster_means <- long_data %>%
  group_by(Cluster, TimeType) %>%
  summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = 'drop')

# Split the data into new, transition (dashed), and old segments
new_segment <- cluster_means %>% filter(TimeType %in% c("30_min new", "1h new", "2h new"))
transition_segment <- cluster_means %>% filter(TimeType %in% c("2h new", "30_min old"))
old_segment <- cluster_means %>% filter(TimeType %in% c("30_min old", "1h old", "2h old"))

# Sample 10 genes per cluster to plot individual gene lines
set.seed(123)  # For reproducibility
sampled_genes_cluster1 <- sample(unique(filter(long_data, Cluster == "Cluster 1")$Gene), 100)
sampled_genes_cluster2 <- sample(unique(filter(long_data, Cluster == "Cluster 2")$Gene), 100)

sampled_data <- long_data %>%
  filter((Cluster == "Cluster 1" & Gene %in% sampled_genes_cluster1) | (Cluster == "Cluster 2" & Gene %in% sampled_genes_cluster2))

# Create the combined plot for both clusters with a blank between segments
p <- ggplot() +
  labs(title = "MRL for Cluster 1 and Cluster 2", x = "Time", y = "MRL") +
  scale_color_manual(values = c("Cluster 1" = "#8887b9", "Cluster 2" = "#95c08b")) +
  scale_x_discrete(drop = FALSE) +  # Retain the blank level in the x-axis
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# Add individual gene lines for both clusters
p <- p + geom_line(data = sampled_data, aes(x = TimeType, y = Expression, color = Cluster, group = interaction(Gene, Cluster)),
                   size = 0.5, alpha = 0.1)

# Draw the three parts of the mean lines for each cluster
# 1. New segment (solid)
p <- p + geom_line(data = new_segment, aes(x = TimeType, y = MeanExpression, color = Cluster, group = Cluster),
                   size = 1.2, linetype = "solid")

# 2. Transition segment from 2h new to 30_min old (dashed)
p <- p + geom_line(data = transition_segment, aes(x = TimeType, y = MeanExpression, color = Cluster, group = Cluster),
                   size = 1.2, linetype = "dashed")

# 3. Old segment (solid)
p <- p + geom_line(data = old_segment, aes(x = TimeType, y = MeanExpression, color = Cluster, group = Cluster),
                   size = 1.2, linetype = "solid")

# Add points for each mean point
fig1E_p <- p + geom_point(data = cluster_means, aes(x = TimeType, y = MeanExpression, color = Cluster), size = 2)

# Display the plot


#pdf(file = "figure1/results/MRL_kinetics.pdf",width = 4, height = 3)
print(fig1E_p)
#dev.off()

#write.csv(sampled_data, file = "fig1E.csv", quote = F, row.names = F)


library(biomaRt)
library(dplyr) # Used for easier data joining at the end
genes_no_version1 <- sub("\\..*", "", sampled_genes_cluster1)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

all_transcripts_info1 <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'cds_length'),
  filters    = 'ensembl_gene_id',
  values     = genes_no_version1,
  mart       = ensembl
)

# 4. PROCESS DATA TO FIND THE LONGEST CDS PER GENE
# ------------------------------------------------
longest_cds_info1 <- all_transcripts_info1 %>%
  filter(!is.na(cds_length)) %>%  # Remove transcripts without a CDS length
  group_by(ensembl_gene_id, external_gene_name) %>% # Group by gene
  summarise(
    longest_cds_length = max(cds_length), # Find the max CDS length in each group
    .groups = 'drop' # Ungroup the data frame after summarizing
  )


# Optional: Join the results back to your original list
# This helps see which genes might not have returned a result.
original_genes_df1 <- data.frame(ensembl_gene_id_version = sampled_genes_cluster1)
original_genes_df1$ensembl_gene_id <- sub("\\..*", "", original_genes_df1$ensembl_gene_id_version)

# Perform a left join to keep all original genes
cds_length_results1 <- left_join(original_genes_df1, longest_cds_info1, by = "ensembl_gene_id")

library(biomaRt)
library(dplyr) # Used for easier data joining at the end
genes_no_version2 <- sub("\\..*", "", sampled_genes_cluster2)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

all_transcripts_info2 <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'cds_length'),
  filters    = 'ensembl_gene_id',
  values     = genes_no_version2,
  mart       = ensembl
)

# 4. PROCESS DATA TO FIND THE LONGEST CDS PER GENE
# ------------------------------------------------
longest_cds_info2 <- all_transcripts_info2 %>%
  filter(!is.na(cds_length)) %>%  # Remove transcripts without a CDS length
  group_by(ensembl_gene_id, external_gene_name) %>% # Group by gene
  summarise(
    longest_cds_length = max(cds_length), # Find the max CDS length in each group
    .groups = 'drop' # Ungroup the data frame after summarizing
  )


# Optional: Join the results back to your original list
# This helps see which genes might not have returned a result.
original_genes_df2 <- data.frame(ensembl_gene_id_version = sampled_genes_cluster2)
original_genes_df2$ensembl_gene_id <- sub("\\..*", "", original_genes_df2$ensembl_gene_id_version)

# Perform a left join to keep all original genes
cds_length_results2 <- left_join(original_genes_df2, longest_cds_info2, by = "ensembl_gene_id")


MRL_new_with_gene_id <- MRL_new
MRL_new_with_gene_id$gene_id <- merge_df_counts_select$gene_id

cds_length_cluster1 <- data.frame(gene_id = cds_length_results1$ensembl_gene_id_version, cds_length = cds_length_results1$longest_cds_length)
cds_length_cluster1_MRL_new <- merge(cds_length_cluster1, MRL_new_with_gene_id, by = "gene_id")

cds_length_cluster2 <- data.frame(gene_id = cds_length_results2$ensembl_gene_id_version, cds_length = cds_length_results2$longest_cds_length)
cds_length_cluster2_MRL_new <- merge(cds_length_cluster2, MRL_new_with_gene_id, by = "gene_id")

list_of_dfs <- list(
  cluster1 = cds_length_cluster1_MRL_new,
  cluster2 = cds_length_cluster2_MRL_new
)

# Use bind_rows() with the .id argument to create the new column
cds_length_clusters_combined_df <- bind_rows(list_of_dfs, .id = "clusters")
cds_length_clusters_combined_filter <- cds_length_clusters_combined_df %>%
  filter(!is.na(cds_length))


cds_length_clusters_combined_filter$density_2h <- cds_length_clusters_combined_filter$NC_2h/cds_length_clusters_combined_filter$cds_length


# Load required libraries
library(ggplot2)
library(dplyr)

# Define your custom theme
base_theme <- theme_bw() + theme(
  legend.position = "none", # Remove legend
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 14, color = "black"), # Change axis numbers to black and larger
  panel.grid.major = element_blank(), # Remove major grid lines
  panel.grid.minor = element_blank(), # Remove minor grid lines
  panel.border = element_rect(color = "black", linewidth = 1) # Change frame color to black
)

# Perform t-test to get p-value
ttest_result <- t.test(density_2h ~ clusters, 
                       data = cds_length_clusters_combined_filter,
                       subset = clusters %in% c("cluster1", "cluster2"))

# Format p-value for display (scientific notation for very small values)
p_value <- ttest_result$p.value
if(p_value < 0.001) {
  p_label <- sprintf("p = %.2e", p_value)
} else {
  p_label <- sprintf("p = %.3f", p_value)
}

# Calculate y position for p-value label (10% above the maximum value)
y_max <- max(cds_length_clusters_combined_filter$density_2h[cds_length_clusters_combined_filter$clusters %in% c("cluster1", "cluster2")], na.rm = TRUE)
y_position <- y_max * 1.1

# Create the box plot with exact p-value
# Create the box plot with outliers removed
p <- ggplot(cds_length_clusters_combined_filter %>%
              filter(clusters %in% c("cluster1", "cluster2")),
            aes(x = clusters, y = density_2h, fill = clusters)) +
  geom_boxplot(alpha = 0.7,
               outlier.shape = NA,  # Remove outliers (set to NA)
               width = 0.6) +
  geom_jitter(width = 0.25,
              alpha = 0.7,
              size = 1.5,
              color = "black") +
  annotate("text",
           x = 1.5,
           y = y_position,
           label = p_label,
           size = 5,
           fontface = "bold") +
  annotate("segment",
           x = 1, xend = 2,
           y = y_max * 1.02, yend = y_max * 1.02,
           color = "black", linewidth = 0.5) +
  labs(title = "Ribosome density of clusters",
       x = "Clusters",
       y = "Ribosome density 2 h") +
  base_theme + # Apply your custom theme
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15))) +
  scale_fill_manual(values = c("cluster1" = "#8484b3", "cluster2" = "#92ba88")) +
  scale_x_discrete(labels = c("cluster1" = "Cluster 1", "cluster2" = "Cluster 2"))


# Display the plot
#pdf(file = "figure1/Ribosome_density_cluster1_2.pdf", width = 3, height = 4)
print(p)
#dev.off()

library(ggplot2)
library(dplyr)

# Prepare data for cumulative distribution - proper ECDF calculation
cumulative_data <- cds_length_clusters_combined_filter %>%
  filter(clusters %in% c("cluster1", "cluster2")) %>%
  arrange(density_2h) %>%
  group_by(clusters) %>%
  mutate(cumulative_freq = (row_number() - 0.5) / n())  # Proper ECDF calculation

# Calculate mean difference for annotation
mean_cluster1 <- mean(cumulative_data$density_2h[cumulative_data$clusters == "cluster1"])
mean_cluster2 <- mean(cumulative_data$density_2h[cumulative_data$clusters == "cluster2"])
mean_diff <- round(mean_cluster2 - mean_cluster1, 3)

# Create the cumulative distribution plot
Fig1G_p <- ggplot(cumulative_data, aes(x = density_2h, y = cumulative_freq, color = clusters)) +
  geom_step(direction = "hv", size = 1.2, alpha = 0.8) +  # Step function with horizontal-vertical direction
  # Statistical annotation in top-left corner
  annotate("text", x = Inf, y = Inf, 
           label = paste0("Mean Diff: ", mean_diff, "\n", "p-value: ", p_label),
           hjust = 1.1, vjust = 1.1, size = 4, fontface = "bold", 
           color = "black") +
  labs(title = "Cumulative Distribution of Ribosome Density by Cluster",
       x = "Ribosome density 2 h",
       y = "Cumulative Probability") +
  base_theme + # Apply your custom theme
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)),
        legend.position = "bottom") +
  scale_color_manual(values = c("cluster1" = "#8484b3", "cluster2" = "#92ba88"),
                     labels = c("cluster1" = "Cluster 1", "cluster2" = "Cluster 2")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)))  # Freed x-axis with minimal padding

# Display the plot
#pdf(file = "figure1/Ribosome_density_cluster1_2_cumulative_plot.pdf", width = 5, height = 4)
print(Fig1G_p)
#dev.off()

#write.csv(cumulative_data, file = "Fig1G.csv", quote = F, row.names = F)


library(ggplot2)
library(gridExtra)

# Provided time points in hours
time_points <- c(0.5, 1, 2)  # Corresponds to 30min, 1h, 2h

# Calculate finite differences (derivative approximation)
# Apply to each row and get the derivative w.r.t. time for new and old data
derivatives_new <- t(apply(MRL_new, 1, function(row) {
  diff(row) / diff(time_points)
}))

derivatives_old <- t(apply(MRL_old, 1, function(row) {
  diff(row) / diff(time_points)
}))

# Convert the result to a data frame and add column names
derivatives_df <- as.data.frame(derivatives_new)
colnames(derivatives_df) <- c("derivative_30min_1h", "derivative_1h_2h")

# Calculate mean derivative for new and old data
mean_derivative_new <- rowMeans(derivatives_new)  # Mean of 30min-1h and 1h-2h for new data
mean_derivative_old <- rowMeans(derivatives_old)  # Mean of 30min-1h and 1h-2h for old data

# Add the mean derivatives for new and old to the dataframe
derivatives_df$mean_derivative_new <- mean_derivative_new
derivatives_df$mean_derivative_old <- mean_derivative_old

# Merge with RNA half-life data by gene_id
derivatives_df$gene_id <- merge_df_counts_select$gene_id
derivatives_df_half_life <- merge(derivatives_df, RNA_features_gene_level, by = "gene_id")

# Filter out rows where RNA_half_life is less than or equal to zero, and remove NA values
filtered_data <- derivatives_df_half_life[complete.cases(derivatives_df_half_life$mean_derivative_new, 
                                                         derivatives_df_half_life$mean_derivative_old,
                                                         derivatives_df_half_life$RNA_half_life) &
                                            derivatives_df_half_life$RNA_half_life > 0, ]

# Compute correlation for mean_derivative_new
cor_value_new <- cor(filtered_data$mean_derivative_new, log2(filtered_data$RNA_half_life))

# Compute correlation for mean_derivative_old
cor_value_old <- cor(filtered_data$mean_derivative_old, log2(filtered_data$RNA_half_life))

# Compute correlation between mean_derivative_new and mean_derivative_old
cor_value_comparison <- cor(filtered_data$mean_derivative_new, filtered_data$mean_derivative_old)


##### change to density plot
library(ggplot2)
library(FNN)  # For k-nearest neighbors

# Calculate correlation and p-value
cor_results <- cor.test(filtered_data$mean_derivative_new, log2(filtered_data$RNA_half_life), method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10  # Number of neighbors to consider
neighbors <- get.knn(data = filtered_data[, c("mean_derivative_new", "RNA_half_life")], k = k)
filtered_data$density <- 1 / rowMeans(neighbors$nn.dist)  # Density = inverse of average distance to neighbors
filtered_data$density <- filtered_data$density / max(filtered_data$density)  # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Plot with density-colored points and regression line
Fig1J_plot_new <- ggplot(filtered_data, aes(x = mean_derivative_new, y = log2(RNA_half_life), color = density)) +
  geom_point(size = 2, alpha = 0.7) +  # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +  # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") +  # Custom color bar
  labs(
    x = "Mean Derivative (New)",
    y = "log2(RNA half-life)",
    title = "Mean Derivative (New) vs log2(RNA Half-life)"
  ) +
  annotate("text", x = 1, y = 1.4, label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +  # Correlation annotation
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

Fig1J_plot_new

NC_2h_new_sum_df <- merge_df_counts_new_select[, c("gene_id", "NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")]
NC_2h_total_sum_df <- merge_df_counts_select[, c("gene_id", "NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")]

# Add the summed column
NC_2h_new_sum_df$NC_2h_counts_new <- rowSums(
  NC_2h_new_sum_df[, c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")],
  na.rm = TRUE
)

NC_2h_total_sum_df$NC_2h_counts_total <- rowSums(
  NC_2h_total_sum_df[, c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")],
  na.rm = TRUE
)

# Keep only gene_id and the new summed column
NC_2h_sum_df <- NC_2h_new_sum_df[, c("gene_id", "NC_2h_counts_new")]
NC_2h_sum_df$NC_2h_counts_total <- NC_2h_total_sum_df[, c("NC_2h_counts_total")]
NC_2h_sum_df$new_total_ratio <- NC_2h_sum_df$NC_2h_counts_new / NC_2h_sum_df$NC_2h_counts_total

filtered_data <- merge(filtered_data, NC_2h_sum_df, by = "gene_id")

##### change to density plot
library(ggplot2)
library(FNN)  # For k-nearest neighbors

# Calculate correlation and p-value
cor_results <- cor.test(filtered_data$mean_derivative_new, filtered_data$new_total_ratio, method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10  # Number of neighbors to consider
neighbors <- get.knn(data = filtered_data[, c("mean_derivative_new", "new_total_ratio")], k = k)
filtered_data$density <- 1 / rowMeans(neighbors$nn.dist)  # Density = inverse of average distance to neighbors
filtered_data$density <- filtered_data$density / max(filtered_data$density)  # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Plot with density-colored points and regression line
Fig1J_plot_new_ratio <- ggplot(filtered_data, aes(x = mean_derivative_new, y = new_total_ratio, color = density)) +
  geom_point(size = 2, alpha = 0.7) +  # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +  # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") +  # Custom color bar
  labs(
    x = "Mean Derivative (New)",
    y = "new_total_ratio",
    title = "Mean Derivative (New) vs new_total_ratio"
  ) +
  annotate("text", x = 1, y = 1.4, label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +  # Correlation annotation
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

Fig1J_plot_new_ratio


#write.csv(filtered_data, file = "Fig1J.csv", quote = F, row.names = F)


# Load required libraries
library(ggplot2)
library(RANN)        # For nearest neighbor calculation
library(circlize)    # For colorRamp2 function

# Set the number of nearest neighbors to use for density estimation
k <- 10

# Calculate the k-nearest neighbors for each point
nn_distances <- RANN::nn2(filtered_data[, c("mean_derivative_new", "mean_derivative_old")], k = k)$nn.dists

# Compute the mean distance to the k-nearest neighbors for each point
# Density is inversely proportional to distance; closer neighbors imply higher density
filtered_data$density <- 1 / rowMeans(nn_distances)

# Normalize density to the range [0, 1] for better control over color mapping
filtered_data$density_scaled <- scales::rescale(filtered_data$density, to = c(0, 1))

# Create custom color scale using colorRamp2
col <- colorRamp2(c(0, 0.3, 1), c("#4ba2dd", "white", "#f04c86"))

# Plot with density-based color using the custom color scale
Fig1F_comparison_plot <- ggplot(filtered_data, aes(x = mean_derivative_new, y = mean_derivative_old)) +
  geom_point(aes(color = density_scaled), size = 3, alpha = 0.7) +  # Map scaled density to color
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +  # Regression line
  labs(x = "Mean Derivative (New)", y = "Mean Derivative (Old)", 
       title = "Comparison of New and Old Mean Derivatives") +  # Labels
  annotate("text", x = Inf, y = Inf, label = paste("Corr =", round(cor_value_comparison, 2)),
           vjust = 2, hjust = 2, size = 5, color = "black") +  # Correlation label
  scale_color_gradientn(colors = col(seq(0, 1, length.out = 100))) +  # Apply custom color map
  theme_bw() +  # Minimal theme
  theme(
    legend.position = "right",  # Show legend for density
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Frame color to black
  )

#write.csv(filtered_data, file = "Fig1F.csv", quote = F, row.names = F)
##### change to density plot
library(ggplot2)
library(FNN)  # For k-nearest neighbors

# Calculate correlation and p-value
cor_results <- cor.test(filtered_data$mean_derivative_new, log2(filtered_data$RNA_half_life), method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10  # Number of neighbors to consider
neighbors <- get.knn(data = filtered_data[, c("mean_derivative_new", "RNA_half_life")], k = k)
filtered_data$density <- 1 / rowMeans(neighbors$nn.dist)  # Density = inverse of average distance to neighbors
filtered_data$density <- filtered_data$density / max(filtered_data$density)  # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Plot with density-colored points and regression line
Fig1I_plot_new <- ggplot(filtered_data, aes(x = mean_derivative_new, y = log2(RNA_half_life), color = density)) +
  geom_point(size = 2, alpha = 0.7) +  # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +  # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") +  # Custom color bar
  labs(
    x = "Mean Derivative (New)",
    y = "log2(RNA half-life)",
    title = "Mean Derivative (New) vs log2(RNA Half-life)"
  ) +
  annotate("text", x = 1, y = 1.4, label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +  # Correlation annotation
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

#write.csv(filtered_data, file = "Fig1I.csv", quote = F, row.names = F)

# Calculate mean values for NC_30min, NC_1h, and NC_2h fractions
mean_values_30min <- data.frame(
  fraction = c("Fraction 1", "Fraction 2", "Fraction 3", "Fraction 4"),
  mean = sapply(merge_df_new_percentage[, c("NC_30min_frac1", "NC_30min_frac2", "NC_30min_frac3", "NC_30min_frac4")], mean, na.rm = TRUE)
)

mean_values_1h <- data.frame(
  fraction = c("Fraction 1", "Fraction 2", "Fraction 3", "Fraction 4"),
  mean = sapply(merge_df_new_percentage[, c("NC_1h_frac1", "NC_1h_frac2", "NC_1h_frac3", "NC_1h_frac4")], mean, na.rm = TRUE)
)

mean_values_2h <- data.frame(
  fraction = c("Fraction 1", "Fraction 2", "Fraction 3", "Fraction 4"),
  mean = sapply(merge_df_new_percentage[, c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")], mean, na.rm = TRUE)
)

# Find the maximum mean value across all datasets
max_value <- max(c(mean_values_30min$mean, mean_values_1h$mean, mean_values_2h$mean), na.rm = TRUE)


# Define the base theme for all plots
fill_colors <- scale_fill_brewer(palette = "Paired")
# Generate plots
FigS1B_p1 <- ggplot(data = mean_values_30min, aes(x = fraction, y = mean, fill = fraction)) +
  geom_bar(stat = "identity") +
  labs(x = "Fraction", y = "Mean Value", title = "NC_30min") +
  ylim(0, max_value) +
  base_theme + fill_colors

FigS1B_p2 <- ggplot(data = mean_values_1h, aes(x = fraction, y = mean, fill = fraction)) +
  geom_bar(stat = "identity") +
  labs(x = "Fraction", title = "NC_1h") +
  ylim(0, max_value) +
  base_theme + fill_colors

FigS1B_p3 <- ggplot(data = mean_values_2h, aes(x = fraction, y = mean, fill = fraction)) +
  geom_bar(stat = "identity") +
  labs(x = "Fraction", title = "NC_2h") +
  ylim(0, max_value) +
  base_theme + fill_colors

#write.csv(mean_values_30min, file = "S1B_p1.csv", quote = F, row.names = F)
#write.csv(mean_values_1h, file = "S1B_p2.csv", quote = F, row.names = F)
#write.csv(mean_values_2h, file = "S1B_p3.csv", quote = F, row.names = F)



# Script to create mismatch plots from three groups of four mutation rate CSV files
# Groups: NC30min, NC1h, NC2h, each with frac[1-4]_overallrates.csv
# Each group generates one figure with four plots in a 1x4 grid


library(ggplot2)
library(gridExtra)

# Define three groups of files and sample names
groups <- list(
  list(
    fileNames = c(
      "data/T_to_C_mut/mutation_rate/NC30min_frac1_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC30min_frac2_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC30min_frac3_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC30min_frac4_overallrates.csv"
    ),
    sampleNames = c("NC30min_frac1", "NC30min_frac2", "NC30min_frac3", "NC30min_frac4"),
    outputFile = "data/T_to_C_mut/results/NC30min_combined.pdf",
    outputData = "data/T_to_C_mut/results/NC30min_combined_plotdata.csv"
  ),
  list(
    fileNames = c(
      "data/T_to_C_mut/mutation_rate/NC1h_frac1_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC1h_frac2_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC1h_frac3_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC1h_frac4_overallrates.csv"
    ),
    sampleNames = c("NC1h_frac1", "NC1h_frac2", "NC1h_frac3", "NC1h_frac4"),
    outputFile = "data/T_to_C_mut/results/NC1h_combined.pdf",
    outputData = "data/T_to_C_mut/results/NC1h_combined_plotdata.csv"
  ),
  list(
    fileNames = c(
      "data/T_to_C_mut/mutation_rate/NC2h_frac1_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC2h_frac2_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC2h_frac3_overallrates.csv",
      "data/T_to_C_mut/mutation_rate/NC2h_frac4_overallrates.csv"
    ),
    sampleNames = c("NC2h_frac1", "NC2h_frac2", "NC2h_frac3", "NC2h_frac4"),
    outputFile = "data/T_to_C_mut/results/NC2h_combined.pdf",
    outputData = "data/T_to_C_mut/results/NC2h_combined_plotdata.csv"
  )
)

# Function to create data frame for a single file
create_data_frame <- function(fileName, sampleName) {
  curTab <- read.table(fileName, stringsAsFactors = FALSE)
  
  # Normalize rates to fractions
  curTab[, c("A", "C", "G", "T")] <- curTab[, c("A", "C", "G", "T")] / rowSums(curTab[, c("A", "C", "G", "T")])
  curTab[, c("a", "c", "g", "t")] <- curTab[, c("a", "c", "g", "t")] / rowSums(curTab[, c("a", "c", "g", "t")])
  
  data.frame(
    Mutation_Label = c("AT", "AT", "AC", "AC", "AG", "AG",
                       "TA", "TA", "TC", "TC", "TG", "TG",
                       "CA", "CA", "CT", "CT", "CG", "CG",
                       "GA", "GA", "GT", "GT", "GC", "GC"),
    Strand = rep(c("+", "-"), 12),
    Mutation_Rate = c(
      curTab["A", "T"], curTab["A", "t"], curTab["A", "C"], curTab["A", "c"], curTab["A", "G"], curTab["A", "g"],
      curTab["T", "A"], curTab["T", "a"], curTab["T", "C"], curTab["T", "c"], curTab["T", "G"], curTab["T", "g"],
      curTab["C", "A"], curTab["C", "a"], curTab["C", "T"], curTab["C", "t"], curTab["C", "G"], curTab["C", "g"],
      curTab["G", "A"], curTab["G", "a"], curTab["G", "T"], curTab["G", "t"], curTab["G", "C"], curTab["G", "c"]
    ),
    File = sampleName
  )
}

# Process each group
for (group in groups) {
  plotList <- list()
  combined_data <- data.frame()  # <- to store data for all samples in this group
  
  for (i in 1:length(group$fileNames)) {
    all_data <- create_data_frame(group$fileNames[i], group$sampleNames[i])
    combined_data <- rbind(combined_data, all_data)  # collect data
    
    curPlot <- ggplot(all_data, aes(x = Mutation_Label, y = Mutation_Rate, fill = Strand)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("+" = "#6baed6", "-" = "#f768a1"), labels = c("+", "-")) +
      scale_y_continuous(limits = c(0, 0.04)) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text = element_blank()
      ) +
      labs(x = "Mutation", y = "Mutation Rate", title = group$sampleNames[i]) +
      facet_grid(. ~ File, scales = "free_x")
    
    plotList[[i]] <- curPlot
  }
  
  # Save plot data for this figure
  #write.csv(combined_data, group$outputData, row.names = FALSE)
  # save the data as original data.
  # Save PDF with plots
  pdf(group$outputFile, width = 15, height = 3.5)
  do.call(grid.arrange, c(plotList, ncol = 4, nrow = 1))
  dev.off()
}



library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

# Function to calculate T2C_per_new_mRNA from directory
calculate_t2c_from_dir <- function(dir_path) {
  files <- list.files(path = dir_path, pattern = "*.tsv", full.names = TRUE)
  
  calculate_t2c <- function(file) {
    df <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 2)
    t2c <- df$ConversionsOnTs / df$TcReadCount
    return(t2c)
  }
  
  results <- map(files, calculate_t2c)
  col_names <- basename(files) %>% 
    sub("\\.merge\\.fq_slamdunk_mapped_filtered_tcount\\.tsv", "", .)
  
  result_df <- data.frame(results)
  colnames(result_df) <- col_names
  return(result_df)
}

# Function to clean combined dataframe by removing rows with NA or 0
clean_combined_dataframe <- function(df_list) {
  # Combine all dataframes
  combined_df <- bind_cols(df_list)
  
  # Remove rows where any value is NA or 0
  cleaned_df <- combined_df %>% 
    filter_all(all_vars(!is.na(.) & . != 0))
  
  # Split back into individual dataframes
  result <- list()
  start_col <- 1
  for (df in df_list) {
    n_cols <- ncol(df)
    result[[length(result) + 1]] <- cleaned_df[, start_col:(start_col + n_cols - 1)]
    start_col <- start_col + n_cols
  }
  return(result)
}

# Function to calculate mean and create pie chart, saving as PDF
plot_pie_to_pdf <- function(df, df_name) {
  # Calculate mean of all columns
  mean_values <- rowMeans(df, na.rm = TRUE)
  
  # Bin the mean values
  bins <- cut(mean_values, 
              breaks = c(-Inf, 0.5, 1, 2, Inf),
              labels = c("0-0.5", "0.5-1", "1-2", ">2"),
              include.lowest = TRUE)
  
  # Create a dataframe for plotting
  bin_counts <- as.data.frame(table(bins))
  colnames(bin_counts) <- c("Bin", "Count")
  
  # Create pie chart
  p <- ggplot(bin_counts, aes(x = "", y = Count, fill = Bin)) +
    geom_bar(stat = "identity", width = 0.4) +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(title = paste("Pie Chart of Mean T2C_per_new_mRNA for", df_name)) +
    theme(legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))
  
  # Save as PDF
  pdf_file <- paste0(df_name, "_pie_chart.pdf")
  pdf(pdf_file, width = 4, height = 4)
  print(p)
  dev.off()
  
  return(mean_values)
}

# Load dataframes
NC_30min <- calculate_t2c_from_dir("data/T_to_C_mut/TC_counts/NC_30min/")
NC_1h <- calculate_t2c_from_dir("data/T_to_C_mut/TC_counts/NC_1h/")
NC_2h <- calculate_t2c_from_dir("data/T_to_C_mut/TC_counts/NC_2h/")

# Combine and clean dataframes
df_list <- list(NC_30min = NC_30min, NC_1h = NC_1h, NC_2h = NC_2h)
cleaned_dfs <- clean_combined_dataframe(df_list)

# Assign cleaned dataframes back
NC_30min <- cleaned_dfs[[1]]
NC_1h <- cleaned_dfs[[2]]
NC_2h <- cleaned_dfs[[3]]

# Verify row counts (should be identical)
cat("Row counts after cleaning:\n")
cat("NC_30min:", nrow(NC_30min), "\n")
cat("NC_1h:", nrow(NC_1h), "\n")
cat("NC_2h:", nrow(NC_2h), "\n")

# Generate pie charts and save as PDFs, also collect means
# FigS1D
mean_30min <- plot_pie_to_pdf(NC_30min, "NC_30min")
mean_1h <- plot_pie_to_pdf(NC_1h, "NC_1h")
mean_2h <- plot_pie_to_pdf(NC_2h, "NC_2h")

# Combine means into a single dataframe
mean_df <- data.frame(NC_30min = mean_30min, NC_1h = mean_1h, NC_2h = mean_2h)

# Save the mean dataframe to CSV
#write.csv(mean_df, "FigS1D.csv", row.names = FALSE)



library(ggplot2)
library(readr)

# Load and process 4sU and DMSO datasets without normalization
data_4sU_1 <- read_csv("data/polysome_profilling_UV/2025.8.13-2h_4sU-1.csv", skip = 32)
data_4sU_1_cut <- data_4sU_1[data_4sU_1$`Distance(mm)` > 8 & data_4sU_1$`Distance(mm)` < 65,]
data_4sU_1_cut$Replicate <- "Replicate 1"
data_4sU_1_cut$Condition <- "4sU"

data_4sU_2 <- read_csv("data/polysome_profilling_UV/2025.8.13-2h_4sU-2.csv", skip = 32)
data_4sU_2_cut <- data_4sU_2[data_4sU_2$`Distance(mm)` > 8 & data_4sU_2$`Distance(mm)` < 65,]
data_4sU_2_cut$Replicate <- "Replicate 2"
data_4sU_2_cut$Condition <- "4sU"

data_4sU_3 <- read_csv("data/polysome_profilling_UV/2025.8.13-2h_4sU-3.csv", skip = 32)
data_4sU_3_cut <- data_4sU_3[data_4sU_3$`Distance(mm)` > 8 & data_4sU_3$`Distance(mm)` < 65,]
data_4sU_3_cut$Replicate <- "Replicate 3"
data_4sU_3_cut$Condition <- "4sU"

data_DMSO_1 <- read_csv("data/polysome_profilling_UV/2025.8.13-2h_DMSO-1.csv", skip = 32)
data_DMSO_1_cut <- data_DMSO_1[data_DMSO_1$`Distance(mm)` > 8 & data_DMSO_1$`Distance(mm)` < 65,]
data_DMSO_1_cut$Replicate <- "Replicate 1"
data_DMSO_1_cut$Condition <- "DMSO"

data_DMSO_2 <- read_csv("data/polysome_profilling_UV/2025.8.13-2h_DMSO-2.csv", skip = 32)
data_DMSO_2_cut <- data_DMSO_2[data_DMSO_2$`Distance(mm)` > 8 & data_DMSO_2$`Distance(mm)` < 65,]
data_DMSO_2_cut$Replicate <- "Replicate 2"
data_DMSO_2_cut$Condition <- "DMSO"

data_DMSO_3 <- read_csv("data/polysome_profilling_UV/2025.8.13-2h_DMSO-3.csv", skip = 32)
data_DMSO_3_cut <- data_DMSO_3[data_DMSO_3$`Distance(mm)` > 8 & data_DMSO_3$`Distance(mm)` < 65,]
data_DMSO_3_cut$Replicate <- "Replicate 3"
data_DMSO_3_cut$Condition <- "DMSO"


data_Puro_1 <- read_csv("data/polysome_profilling_UV/20230321_WT_puro_2h_4.csv", skip = 32)
data_Puro_1_cut <- data_Puro_1[data_Puro_1$`Distance(mm)` > 8 & data_Puro_1$`Distance(mm)` < 65,]
data_Puro_1_cut$Replicate <- "Replicate 1"
data_Puro_1_cut$Condition <- "Puro"

data_Puro_2 <- read_csv("data/polysome_profilling_UV/20230321_WT_puro_2h_2.csv", skip = 32)
data_Puro_2_cut <- data_Puro_2[data_Puro_2$`Distance(mm)` > 8 & data_Puro_2$`Distance(mm)` < 65,]
data_Puro_2_cut$Replicate <- "Replicate 2"
data_Puro_2_cut$Condition <- "Puro"

data_Puro_3 <- read_csv("data/polysome_profilling_UV/20230321_WT_puro_2h_3.csv", skip = 32)
data_Puro_3_cut <- data_Puro_3[data_Puro_3$`Distance(mm)` > 8 & data_Puro_3$`Distance(mm)` < 65,]
data_Puro_3_cut$Replicate <- "Replicate 3"
data_Puro_3_cut$Condition <- "Puro"


# Combine all datasets into one with an additional column to differentiate between replicates and conditions
combined_data <- rbind(data_4sU_1_cut, data_4sU_2_cut, data_4sU_3_cut, 
                       data_DMSO_1_cut, data_DMSO_2_cut, data_DMSO_3_cut, data_Puro_1_cut, data_Puro_2_cut, data_Puro_3_cut)

# Plot
FigS1C_p <- ggplot(combined_data, aes(x = `Distance(mm)`, y = Absorbance, color = Condition, group = interaction(Condition, Replicate))) +
  geom_line(size = 0.5) +  # Adjust the line width here
  scale_color_manual(values = c("4sU" = "#ec73a8", "DMSO" = "#5db480", "Puro" = "#000000")) +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 12, colour = "black")
  ) +
  labs(title = "Polysome profilling of 500uM 4sU vs DMSO 2h treatment", x = "Distance (mm)", y = "UV Absorbance a.u.")

FigS1C_p

polysome_profilling_data <- bind_rows(data_4sU_1_cut, data_4sU_2_cut, data_4sU_3_cut, 
                                      data_DMSO_1_cut, data_DMSO_2_cut, data_DMSO_3_cut,
                                      data_Puro_1_cut, data_Puro_2_cut, data_Puro_3_cut)

#write.csv(polysome_profilling_data, file = "FigS1C.csv", quote = F, row.names = F)


load("data/data_for_rep_test.RData")

# Assuming your dataset is named NC_2h_rep1_counts
merged_NC_2h_rep1_counts <- NC_2h_rep1_counts %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum))

merged_NC_2h_rep2_counts <- NC_2h_rep2_counts %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum))


merged_NC_2h_rep1_counts_new <- NC_2h_rep1_counts_new %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum))

merged_NC_2h_rep2_counts_new <- NC_2h_rep2_counts_new %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum))


merged_NC_2h_rep1_counts_select <- merged_NC_2h_rep1_counts %>%
  filter(gene_id %in% merge_df_counts_select$gene_id)

merged_NC_2h_rep2_counts_select <- merged_NC_2h_rep2_counts %>%
  filter(gene_id %in% merge_df_counts_select$gene_id)


merged_NC_2h_rep1_counts_select_new <- merged_NC_2h_rep1_counts_new %>%
  filter(gene_id %in% merge_df_counts_select$gene_id)

merged_NC_2h_rep2_counts_select_new <- merged_NC_2h_rep2_counts_new %>%
  filter(gene_id %in% merge_df_counts_select$gene_id)


merged_NC_2h_rep1_counts_select_old <- data.frame(gene_id = merged_NC_2h_rep1_counts_select$gene_id, merged_NC_2h_rep1_counts_select[,-1] - merged_NC_2h_rep1_counts_select_new[,-1])
merged_NC_2h_rep2_counts_select_old <- data.frame(gene_id = merged_NC_2h_rep2_counts_select$gene_id, merged_NC_2h_rep2_counts_select[,-1] - merged_NC_2h_rep2_counts_select_new[,-1])

MRLs_rep1 <- calculate_polysome_load(merged_NC_2h_rep1_counts_select_new, merged_NC_2h_rep1_counts_select_old, c("NC_2h"))
MRLs_rep2 <- calculate_polysome_load(merged_NC_2h_rep2_counts_select_new, merged_NC_2h_rep2_counts_select_old, c("NC_2h"))

library(ggplot2)

# Remove missing values before calculating correlations
# Also ensure there are no constant values in the columns

# Prepare data for correlation
data_1_2 <- data.frame(x = MRLs_rep1$new_polysome_load$NC_2h, y = MRLs_rep2$new_polysome_load$NC_2h)

# Filter out rows with any NA
data_1_2 <- na.omit(data_1_2)

# Compute correlations, checking for constant vectors
corr_1_2 <- ifelse(sd(data_1_2$x) == 0 | sd(data_1_2$y) == 0, NA, cor(data_1_2$x, data_1_2$y))

# Plot 1: MRLs_rep1 vs MRLs_rep2
FigS1E_p_rep_new <- ggplot(data = data_1_2, aes(x = x, y = y)) +
  geom_point(color = '#4ba2dd', alpha = 0.1) +
  geom_smooth(method = 'lm', se = FALSE, color = '#08306b') +
  labs(title = "new mRNA MRL_rep1 vs MRLs_rep2", 
       x = "NC 2h rep1", 
       y = "NC 2h rep2") +
  annotate("text", x = Inf, y = Inf, label = paste("Corr:", round(corr_1_2, 2)), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black")+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

#write.csv(data_1_2, file = "FigS1E_new.csv", quote = F, row.names = F)


# Prepare data for correlation
data_1_2 <- data.frame(x = MRLs_rep1$old_polysome_load$NC_2h, y = MRLs_rep2$old_polysome_load$NC_2h)

# Filter out rows with any NA
data_1_2 <- na.omit(data_1_2)

# Compute correlations, checking for constant vectors
corr_1_2 <- ifelse(sd(data_1_2$x) == 0 | sd(data_1_2$y) == 0, NA, cor(data_1_2$x, data_1_2$y))

# Plot 1: MRLs_rep1 vs MRLs_rep2
FigS1E_p_rep_old <- ggplot(data = data_1_2, aes(x = x, y = y)) +
  geom_point(color = '#4ba2dd', alpha = 0.1) +
  geom_smooth(method = 'lm', se = FALSE, color = '#08306b') +
  labs(title = "old mRNA MRL_rep1 vs MRLs_rep2", 
       x = "NC 2h rep1", 
       y = "NC 2h rep2") +
  annotate("text", x = Inf, y = Inf, label = paste("Corr:", round(corr_1_2, 2)), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black")+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

FigS1E_p_rep_new

FigS1E_p_rep_old

#write.csv(data_1_2, file = "FigS1E_old.csv", quote = F, row.names = F)


merge_df_counts_new_select
merge_df_counts_old_select

# Load necessary libraries
library(org.Hs.eg.db)  # Human genome annotation package
library(AnnotationDbi) # For mapIds function

# Function to map Ensembl IDs to Gene Symbols and add them as a new column
map_ensembl_to_gene_symbol <- function(df, id_column_name = "gene_id") {
  # Extract Ensembl IDs from the specified column
  ensembl_ids <- df[[id_column_name]]
  
  # Remove version numbers (everything after the dot in Ensembl IDs)
  ensembl_ids <- sub("\\..*", "", ensembl_ids)
  
  # Map Ensembl IDs to Gene Symbols using org.Hs.eg.db
  gene_symbols <- mapIds(org.Hs.eg.db, 
                         keys = ensembl_ids, 
                         column = "SYMBOL", 
                         keytype = "ENSEMBL", 
                         multiVals = "first")
  
  # Add the gene symbol as a new column in the dataframe
  df$gene_name <- gene_symbols
  
  # Return the modified dataframe
  return(df)
}

# Apply the function to your dataframes
merge_df_counts_new_select <- map_ensembl_to_gene_symbol(merge_df_counts_new_select, "gene_id")
merge_df_counts_old_select <- map_ensembl_to_gene_symbol(merge_df_counts_old_select, "gene_id")

# If you want to remove entries without valid gene symbols (optional)
merge_df_counts_new_select <- merge_df_counts_new_select[!is.na(merge_df_counts_new_select$gene_name), ]
merge_df_counts_old_select <- merge_df_counts_old_select[!is.na(merge_df_counts_old_select$gene_name), ]

library(ggplot2)
library(gridExtra)
library(reshape2)

# Function to normalize and plot the gene expression data for both new and old data
plot_gene_expression_combined <- function(gene_name, df_new, df_old) {
  
  # Subset the dataframe for the specific gene in both new and old datasets
  gene_data_new <- df_new[df_new$gene_name == gene_name, ]
  gene_data_old <- df_old[df_old$gene_name == gene_name, ]
  
  # Check if the gene is found in both datasets
  if (nrow(gene_data_new) == 0 || nrow(gene_data_old) == 0) {
    stop("Gene not found in one or both dataframes")
  }
  
  # Combine the 8 fractions for each time point across new and old
  combined_NC_30min <- c(as.numeric(gene_data_new[, c("NC_30min_frac1", "NC_30min_frac2", "NC_30min_frac3", "NC_30min_frac4")]),
                         as.numeric(gene_data_old[, c("NC_30min_frac1", "NC_30min_frac2", "NC_30min_frac3", "NC_30min_frac4")]))
  
  combined_NC_1h <- c(as.numeric(gene_data_new[, c("NC_1h_frac1", "NC_1h_frac2", "NC_1h_frac3", "NC_1h_frac4")]),
                      as.numeric(gene_data_old[, c("NC_1h_frac1", "NC_1h_frac2", "NC_1h_frac3", "NC_1h_frac4")]))
  
  combined_NC_2h <- c(as.numeric(gene_data_new[, c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")]),
                      as.numeric(gene_data_old[, c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")]))
  
  # Normalize each time point across new and old combined data
  max_NC_30min <- max(combined_NC_30min, na.rm = TRUE)
  max_NC_1h <- max(combined_NC_1h, na.rm = TRUE)
  max_NC_2h <- max(combined_NC_2h, na.rm = TRUE)
  
  # Create data frames for both new and old, normalizing by the combined max values
  expression_data_new <- data.frame(
    Fraction = factor(c("frac1", "frac2", "frac3", "frac4")),
    NC_30min = as.numeric(gene_data_new[, c("NC_30min_frac1", "NC_30min_frac2", "NC_30min_frac3", "NC_30min_frac4")]) / max_NC_30min,
    NC_1h = as.numeric(gene_data_new[, c("NC_1h_frac1", "NC_1h_frac2", "NC_1h_frac3", "NC_1h_frac4")]) / max_NC_1h,
    NC_2h = as.numeric(gene_data_new[, c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")]) / max_NC_2h
  )
  
  expression_data_old <- data.frame(
    Fraction = factor(c("frac1", "frac2", "frac3", "frac4")),
    NC_30min = as.numeric(gene_data_old[, c("NC_30min_frac1", "NC_30min_frac2", "NC_30min_frac3", "NC_30min_frac4")]) / max_NC_30min,
    NC_1h = as.numeric(gene_data_old[, c("NC_1h_frac1", "NC_1h_frac2", "NC_1h_frac3", "NC_1h_frac4")]) / max_NC_1h,
    NC_2h = as.numeric(gene_data_old[, c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4")]) / max_NC_2h
  )
  
  # Melt the data for ggplot2 compatibility
  expression_data_new_melt <- melt(expression_data_new, id.vars = "Fraction")
  expression_data_old_melt <- melt(expression_data_old, id.vars = "Fraction")
  
  # Create line plots for new and old datasets
  plot_new <- ggplot(expression_data_new_melt, aes(x = Fraction, y = value, color = variable, group = variable)) +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    scale_color_manual(values = c("NC_30min" = "#6baed6", "NC_1h" = "#2171b5", "NC_2h" = "#f768a1")) +  # Custom colors for each line
    labs(title = paste("Expression of", gene_name, "(New)"),
         x = "Fraction",
         y = "Normalized Expression",
         color = "Condition") +
    theme_bw() +  # Use a minimal theme
    theme(
      legend.position = "none",  # Remove legend
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
    )
  
  plot_old <- ggplot(expression_data_old_melt, aes(x = Fraction, y = value, color = variable, group = variable)) +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    scale_color_manual(values = c("NC_30min" = "#6baed6", "NC_1h" = "#2171b5", "NC_2h" = "#f768a1")) +  # Custom colors for each line
    labs(title = paste("Expression of", gene_name, "(Old)"),
         x = "Fraction",
         y = "Normalized Expression",
         color = "Condition") +
    theme_bw() +  # Use a minimal theme
    theme(
      legend.position = "none",  # Remove legend
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
    )
  
  # Return both plots as a list
  return(list(plot_new = plot_new, plot_old = plot_old))
}

# Generate plots for different genes
plots_FKBP4 <- plot_gene_expression_combined("FKBP4", merge_df_counts_new_select, merge_df_counts_old_select)
plots_JUN <- plot_gene_expression_combined("JUN", merge_df_counts_new_select, merge_df_counts_old_select)
plots_SOX4 <- plot_gene_expression_combined("SOX4", merge_df_counts_new_select, merge_df_counts_old_select)
plots_ACADS <- plot_gene_expression_combined("ACADS", merge_df_counts_new_select, merge_df_counts_old_select)
plots_PSAP <- plot_gene_expression_combined("PSAP", merge_df_counts_new_select, merge_df_counts_old_select)
plots_NFKB <- plot_gene_expression_combined("NFKBIA", merge_df_counts_new_select, merge_df_counts_old_select)


#FigS1F

# Display the plots in a grid
#pdf(file = "../figure1_output/polysome_distribution_FKBP.pdf", width = 8, height = 4)
grid.arrange(plots_FKBP4$plot_new, plots_FKBP4$plot_old, ncol = 2)
#dev.off()
#pdf(file = "../figure1_output/polysome_distribution_JUN.pdf", width = 8, height = 4)
grid.arrange(plots_JUN$plot_new, plots_JUN$plot_old, ncol = 2)
#dev.off()
#pdf(file = "../figure1_output/polysome_distribution_SOX4.pdf", width = 8, height = 4)
grid.arrange(plots_SOX4$plot_new, plots_SOX4$plot_old, ncol = 2)
#dev.off()
#pdf(file = "../figure1_output/polysome_distribution_ACADS.pdf", width = 8, height = 4)
grid.arrange(plots_ACADS$plot_new, plots_ACADS$plot_old, ncol = 2)
#dev.off()
#pdf(file = "../figure1_output/polysome_distribution_PSAP.pdf", width = 8, height = 4)
grid.arrange(plots_PSAP$plot_new, plots_PSAP$plot_old, ncol = 2)
#dev.off()
#pdf(file = "../figure1_output/polysome_distribution_NFKB.pdf", width = 8, height = 4)
grid.arrange(plots_NFKB$plot_new, plots_NFKB$plot_old, ncol = 2)
#dev.off()


selected_genes <- c("FKBP4", "JUN", "SOX4", "ACADS", "PSAP", "NFKBIA")

# Subset the dataframe for those genes
selected_rows_merge_df_counts_new_select <- merge_df_counts_new_select[merge_df_counts_new_select$gene_name %in% selected_genes, ]
selected_rows_merge_df_counts_old_select <- merge_df_counts_old_select[merge_df_counts_old_select$gene_name %in% selected_genes, ]

#write.csv(selected_rows_merge_df_counts_new_select, file = "FigS1F_new.csv", quote = F, row.names = F)
#write.csv(selected_rows_merge_df_counts_old_select, file = "FigS1F_old.csv", quote = F, row.names = F)










##################### add ribosome density

library(dplyr)
library(biomaRt)

# Connect to the Ensembl database
ensembl_mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get longest CDS length from Ensembl for the genes in your dataframe
longest_cds_info <- getBM(
  attributes = c('ensembl_gene_id', 'cds_length'),
  filters    = 'ensembl_gene_id',
  values     = sub("\\..*", "", filtered_data$gene_id), # Remove version for query
  mart       = ensembl_mart
) %>%
  filter(!is.na(cds_length)) %>%
  group_by(ensembl_gene_id) %>%
  summarise(cds_length = max(cds_length), .groups = 'drop')

#save(longest_cds_info, file = "data/figure1_longest_cds_info.RData")


# Add the new 'cds_length' column to your original dataframe
filtered_data_cds_length <- filtered_data %>%
  mutate(ensembl_gene_id = sub("\\..*", "", gene_id)) %>%
  left_join(longest_cds_info, by = "ensembl_gene_id") %>%
  dplyr::select(-ensembl_gene_id) # Clean up the temporary join key

names(filtered_data_cds_length)[names(filtered_data_cds_length) == "cds_length.y"] <- "cds_length"

MRL_new_gene_id <- MRL_new
MRL_new_gene_id_total <- MRL_total$old_polysome_load
colnames(MRL_new_gene_id_total) <- c("NC_30min_total", "NC_1h_total", "NC_2h_total")
MRL_new_gene_id <- cbind(MRL_new_gene_id, MRL_new_gene_id_total)

MRL_new_gene_id$gene_id <- merge_df_counts_select$gene_id

filtered_data_cds_length_MRL <- merge(filtered_data_cds_length, MRL_new_gene_id, by = "gene_id")
filtered_data_cds_length_MRL$ribosome_density_2h <- filtered_data_cds_length_MRL$NC_2h/filtered_data_cds_length_MRL$cds_length
filtered_data_cds_length_MRL$ribosome_density_2h_total <- filtered_data_cds_length_MRL$NC_2h_total/filtered_data_cds_length_MRL$cds_length

plot(log2(filtered_data_cds_length_MRL$RNA_half_life), filtered_data_cds_length_MRL$ribosome_density_2h)

library(ggplot2)
library(FNN) # For k-nearest neighbors

# Remove rows with NAs for the variables of interest
data_clean <- filtered_data_cds_length_MRL[complete.cases(filtered_data_cds_length_MRL[, c("RNA_half_life", "ribosome_density_2h")]), ]

# Calculate correlation and p-value
cor_results <- cor.test(data_clean$RNA_half_life, 
                        data_clean$ribosome_density_2h, 
                        method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10 # Number of neighbors to consider
neighbors <- get.knn(data = data_clean[, c("RNA_half_life", "ribosome_density_2h")], k = k)
data_clean$density <- 1 / rowMeans(neighbors$nn.dist) # Density = inverse of average distance to neighbors
data_clean$density <- data_clean$density / max(data_clean$density) # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Plot with density-colored points and regression line
FigS2G_ribosome_density_plot <- ggplot(data_clean, 
                                aes(x = log2(RNA_half_life), 
                                    y = ribosome_density_2h, 
                                    color = density)) +
  geom_point(size = 2, alpha = 0.7) + # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) + # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") + # Custom color bar
  labs(
    x = "log2(RNA half-life)",
    y = "Ribosome Density (2h)",
    title = "log2(RNA Half-life) vs Ribosome Density (2h)"
  ) +
  annotate("text", x = min(log2(data_clean$RNA_half_life)), 
           y = max(data_clean$ribosome_density_2h), 
           label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 0, vjust = 1, size = 5, color = "black", fontface = "bold") + # Correlation annotation
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# Display the plot
print(FigS2G_ribosome_density_plot)

write.csv(data_clean, file = "FigS2G.csv", quote = F, row.names = F)

library(ggplot2)
library(FNN) # For k-nearest neighbors

# Remove rows with NAs for the variables of interest
data_clean <- filtered_data_cds_length_MRL[complete.cases(filtered_data_cds_length_MRL[, c("RNA_half_life", "ribosome_density_2h_total")]), ]

# Calculate correlation and p-value
cor_results <- cor.test(data_clean$RNA_half_life, 
                        data_clean$ribosome_density_2h_total, 
                        method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10 # Number of neighbors to consider
neighbors <- get.knn(data = data_clean[, c("RNA_half_life", "ribosome_density_2h_total")], k = k)
data_clean$density <- 1 / rowMeans(neighbors$nn.dist) # Density = inverse of average distance to neighbors
data_clean$density <- data_clean$density / max(data_clean$density) # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Plot with density-colored points and regression line
ribosome_density_plot_total_2h <- ggplot(data_clean, 
                                         aes(x = log2(RNA_half_life), 
                                             y = ribosome_density_2h_total, 
                                             color = density)) +
  geom_point(size = 2, alpha = 0.7) + # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) + # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") + # Custom color bar
  labs(
    x = "log2(RNA half-life)",
    y = "Ribosome Density (2h total)",
    title = "log2(RNA Half-life) vs Ribosome Density (2h total)"
  ) +
  annotate("text", x = min(log2(data_clean$RNA_half_life)), 
           y = max(data_clean$ribosome_density_2h_total), 
           label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 0, vjust = 1, size = 5, color = "black", fontface = "bold") + # Correlation annotation
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# Display the plot
pdf(file = "figure1/Ribosome_density_total_2h_vs_RNA_half_life.pdf", width = 4.5, height = 3)
print(ribosome_density_plot_total_2h)
dev.off()


load("figure1E_cluster.RData")
head(sampled_genes_cluster1)
head(sampled_genes_cluster2)

library(ggplot2)
library(ggrepel) # For better text labels

# Remove rows with NAs for the variables of interest
data_clean <- filtered_data_cds_length_MRL[complete.cases(filtered_data_cds_length_MRL[, c("ribosome_density_2h_total", "ribosome_density_2h")]), ]

# Define the two clusters of interest
cluster1_genes <- sampled_genes_cluster1
cluster2_genes <- sampled_genes_cluster2

# Create a simple cluster identifier
data_clean$cluster_type <- "Other"

# Assign cluster 1 to genes in sampled_genes_cluster1
data_clean$cluster_type[data_clean$gene_id %in% cluster1_genes] <- "Cluster 1"

# Assign cluster 2 to genes in sampled_genes_cluster2
data_clean$cluster_type[data_clean$gene_id %in% cluster2_genes] <- "Cluster 2"

# Convert to factor for proper plotting
data_clean$cluster_type <- factor(data_clean$cluster_type,
                                  levels = c("Other", "Cluster 1", "Cluster 2"))

# Calculate correlation and p-value for all data (for the overall annotation)
cor_results_all <- cor.test(data_clean$ribosome_density_2h_total,
                            data_clean$ribosome_density_2h,
                            method = "pearson")
cor_value_all <- round(cor_results_all$estimate, 3)
p_value_all <- signif(cor_results_all$p.value, 3)

# Calculate correlation and p-value for Cluster 1 only
cluster1_data <- data_clean[data_clean$cluster_type == "Cluster 1", ]
if(nrow(cluster1_data) > 1) {
  cor_results_cluster1 <- cor.test(cluster1_data$ribosome_density_2h_total,
                                   cluster1_data$ribosome_density_2h,
                                   method = "pearson")
  cor_value_cluster1 <- round(cor_results_cluster1$estimate, 3)
  p_value_cluster1 <- signif(cor_results_cluster1$p.value, 3)
} else {
  cor_value_cluster1 <- NA
  p_value_cluster1 <- NA
}

# Calculate correlation and p-value for Cluster 2 only
cluster2_data <- data_clean[data_clean$cluster_type == "Cluster 2", ]
if(nrow(cluster2_data) > 1) {
  cor_results_cluster2 <- cor.test(cluster2_data$ribosome_density_2h_total,
                                   cluster2_data$ribosome_density_2h,
                                   method = "pearson")
  cor_value_cluster2 <- round(cor_results_cluster2$estimate, 3)
  p_value_cluster2 <- signif(cor_results_cluster2$p.value, 3)
} else {
  cor_value_cluster2 <- NA
  p_value_cluster2 <- NA
}

# Determine limits for same scale
min_val <- min(c(data_clean$ribosome_density_2h_total, data_clean$ribosome_density_2h), na.rm = TRUE)
max_val <- max(c(data_clean$ribosome_density_2h_total, data_clean$ribosome_density_2h), na.rm = TRUE)

# Plot for Cluster 1
FigS1G_ribosome_density_new_vs_total_cluster <- ggplot(data_clean,
                                                aes(x = ribosome_density_2h, # Swapped: now y becomes x
                                                    y = ribosome_density_2h_total)) + # Swapped: now x becomes y
  # First layer: All other genes (background) in grey
  geom_point(data = subset(data_clean, cluster_type == "Other"),
             color = "grey70", size = 1.5, alpha = 0.6) +
  # Second layer: Cluster 1 genes on top
  geom_point(data = subset(data_clean, cluster_type == "Cluster 1"),
             color = "#8887b9", size = 2.5, alpha = 0.8) +
  # Third layer: Cluster 2 genes on top (also grey for this plot)
  geom_point(data = subset(data_clean, cluster_type == "Cluster 2"),
             color = "#95c08b", size = 2.5, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 1) + # Diagonal line
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) + # Regression line
  scale_x_continuous(limits = c(0, 0.02)) +
  scale_y_continuous(limits = c(0, 0.02)) +
  coord_fixed(ratio = 1) + # Ensure same scale
  labs(
    x = "MRD new 2 h", # Swapped label
    y = "MRD total 2 h", # Swapped label
    title = "MRD on steady-state vs. MRD on new mRNA (Cluster 1)"
  ) +
  # Overall correlation annotation
  annotate("text", x = 0.005,
           y = 0.018,
           label = paste0("Overall: r = ", cor_value_all, "\nP = ", p_value_all),
           hjust = 0, vjust = 1, size = 4, color = "black", fontface = "bold") +
  # Cluster 1 specific correlation annotation
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# Display the plot for Cluster 1
#pdf(file = "figure1/Ribosome_density_new_2h_vs_total_2h_cluster_2025.9.29.pdf", width = 5, height = 4)
print(FigS1G_ribosome_density_new_vs_total_cluster)
#dev.off()

data_clean_output <- data_clean[,c("ribosome_density_2h_total", "ribosome_density_2h")]

write.csv(data_clean_output, file = "FigS1G.csv", quote = F, row.names = F)

