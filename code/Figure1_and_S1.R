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

#MRL_total <- calculate_polysome_load(merge_df_counts_select, c("NC_30min", "NC_1h", "NC_2h"))
MRLs <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))
MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]
MRL_total <- MRL_new + MRL_old

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


# Compare new and old mRNA loading to polysome dynamics with mRNA halflife

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



