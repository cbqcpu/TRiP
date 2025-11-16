source("code/polysome_seq_functions.R")
load("data/HEK293T_preprocessed.RData")
load("data/HEK293T_pulse_chase2.RData")
load("data/figure1_heatmap_ordered_cluster.RData")

# Pulse Only experiment
merge_df_counts_select <- merge_df_counts[rowMeans(merge_df_counts[, -1]) > 1000, ]
merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts$gene_id %in% merge_df_counts_select$gene_id, ]
merge_df_counts_old_select <- data.frame(gene_id = merge_df_counts_select$gene_id, merge_df_counts_select[, -1] - merge_df_counts_new_select[, -1])

# Pulse-Chase experiment
merge_df_counts_select_pulse_chase <- merge_df_counts_pulse_chase[rowMeans(merge_df_counts[, -1]) > 1000, ]
merge_df_counts_new_select_pulse_chase <- merge_df_counts_new_pulse_chase[merge_df_counts$gene_id %in% merge_df_counts_select$gene_id, ]
merge_df_counts_old_select_pulse_chase <- data.frame(gene_id = merge_df_counts_select_pulse_chase$gene_id, merge_df_counts_select_pulse_chase[, -1] - merge_df_counts_new_select_pulse_chase[, -1])


MRLs <- calculate_polysome_load(cbind(merge_df_counts_new_select, merge_df_counts_new_select_pulse_chase), cbind(merge_df_counts_old_select, merge_df_counts_old_select_pulse_chase),
                                c("NC_30min", "NC_1h", "NC_2h", "pulse2h_chase1h", "pulse2h_chase2h", "pulse2h_chase8h")
)

MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]


library(readr)
set.seed(1234) # You can use any integer as the seed


# Calculate the max value for each row (gene) using ONLY the new mRNA data
row_max <- apply(MRL_new, 1, max)


MRL_new_normalized <- MRL_new / row_max
# The following line is no longer strictly necessary for the heatmap but is kept for completeness
MRL_old_normalized <- MRL_old / row_max

colnames(MRL_new_normalized) <- c("new_30_min", "new_1h", "new_2h", "new_pulse2h_chase1h", "new_pulse2h_chase2h", "new_pulse2h_chase8h")
colnames(MRL_old_normalized) <- c("old_30_min", "old_1h", "old_2h", "old_pulse2h_chase1h", "old_pulse2h_chase2h", "old_pulse2h_chase8h")

# This line remains the same from the last step, ensuring we only plot new mRNA
data_matrix <- as.matrix(MRL_new_normalized)
rownames(data_matrix) <- merge_df_counts_select$gene_id


library(Mfuzz)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)

# Sort the assigned_clusters and reorder the data accordingly
sorted_data <- data_matrix[sorted_indices, ]

# Define custom colors for clusters 1 and 2
cluster_colors <- c("1" = "#8887b9", # Color for Cluster 1
                    "2" = "#95c08b") # Color for Cluster 2

# Convert your clusters to a factor to make sure colors are assigned correctly
sorted_clusters <- as.factor(sorted_clusters)

# Update the row annotation based on the sorted clusters with custom colors
ha <- rowAnnotation(Cluster = sorted_clusters,
                    col = list(Cluster = cluster_colors),
                    annotation_legend_param = list(Cluster = list(title = "Cluster"))
)

# Create the heatmap
FigS2B_phm <- Heatmap(sorted_data,
                   name = "Normalized MRL",
                   left_annotation = ha,
                   col = colorRamp2(c(0, 0.5, 1),
                                    c("white", "orange", "darkred")
                   ),
                   show_row_names = FALSE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE
)

FigS2B_phm

write.csv(sorted_data, file = "FigS2B.csv", quote = F, row.names = F)


# test if not normalize by row, what the two cluster like

MRL_new_no_normalized <- MRL_new
MRL_old_no_normalized <- MRL_old
colnames(MRL_new_no_normalized) <- c("new_30_min", "new_1h", "new_2h", "new_pulse2h_chase1h", "new_pulse2h_chase2h", "new_pulse2h_chase8h")
colnames(MRL_old_no_normalized) <- c("old_30_min", "old_1h", "old_2h", "old_pulse2h_chase1h", "old_pulse2h_chase2h", "old_pulse2h_chase8h")
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


# Load required libraries
# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Convert 'sorted_clusters' to a data frame for merging
cluster_assignments <- data.frame(Gene = names(sorted_clusters), Cluster = as.integer(sorted_clusters))

# Merge the data with cluster assignments
merged_data <- merge(data_matrix_no, cluster_assignments, by = "Gene")

# Reshape data to long format for plotting individual gene lines
long_data <- merged_data %>%
  pivot_longer(cols = starts_with("new_"),  # <- keep only "new"
               names_to = c("Type", "Time"),
               names_pattern = "(new)_(.+)", values_to = "Expression") %>%
  mutate(TimeType = factor(paste(Time, Type),
                           levels = c("30_min new", "1h new", "2h new",
                                      "pulse2h_chase1h new", "pulse2h_chase2h new", "pulse2h_chase8h new")),
         Cluster = factor(Cluster, levels = c(1, 2), labels = c("Cluster 1", "Cluster 2")))

# Calculate mean values for 'new' mRNA across clusters and time points
cluster_means <- long_data %>%
  group_by(Cluster, TimeType) %>%
  summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = 'drop')

# Sample 100 genes per cluster to plot individual gene lines
set.seed(123)  # For reproducibility
sampled_genes_cluster1 <- sample(unique(filter(long_data, Cluster == "Cluster 1")$Gene), 100)
sampled_genes_cluster2 <- sample(unique(filter(long_data, Cluster == "Cluster 2")$Gene), 100)

sampled_data <- long_data %>%
  filter((Cluster == "Cluster 1" & Gene %in% sampled_genes_cluster1) |
           (Cluster == "Cluster 2" & Gene %in% sampled_genes_cluster2))

# Create the combined plot for both clusters
Fig_p <- ggplot() +
  labs(title = "MRL for Cluster 1 and Cluster 2", x = "Time", y = "MRL") +
  scale_color_manual(values = c("Cluster 2" = "#95c08b", "Cluster 1" = "#8887b9")) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Add individual gene lines (only "new")
Fig_p <- Fig_p + geom_line(data = sampled_data,
                           aes(x = TimeType, y = Expression, color = Cluster, group = interaction(Gene, Cluster)),
                           size = 0.5, alpha = 0.1)

# Add mean line (only "new", solid)
Fig_p <- Fig_p + geom_line(data = cluster_means,
                           aes(x = TimeType, y = MeanExpression, color = Cluster, group = Cluster),
                           size = 1.2, linetype = "solid")

# Add mean points (only "new")
Fig_p <- Fig_p + geom_point(data = cluster_means,
                            aes(x = TimeType, y = MeanExpression, color = Cluster), size = 2)
FigS2C <- Fig_p
# Display the plot
print(FigS2C)

write.csv(sorted_data, file = "FigS2C.csv", quote = F, row.names = F)


time_points_pulse <- c(0.5, 1, 2)  # Corresponds to 30min, 1h, 2h

# Calculate finite differences (derivative approximation)
# Apply to each row and get the derivative w.r.t. time for new and old data
derivatives_new_pulse <- t(apply(MRL_new[,1:3], 1, function(row) {
  diff(row) / diff(time_points_pulse)
}))

derivatives_old_pulse <- t(apply(MRL_old[,1:3], 1, function(row) {
  diff(row) / diff(time_points_pulse)
}))

time_points_chase <- c(2, 3, 4, 10)
derivatives_new_chase <- t(apply(MRL_new[,3:6], 1, function(row) {
  diff(row) / diff(time_points_chase)
}))

# Convert the result to a data frame and add column names
derivatives_df_new_pulse <- as.data.frame(derivatives_new_pulse)
colnames(derivatives_df_new_pulse) <- c("derivative_30min_1h", "derivative_1h_2h")

derivatives_df_old_pulse <- as.data.frame(derivatives_old_pulse)
colnames(derivatives_df_old_pulse) <- c("derivative_30min_1h", "derivative_1h_2h")

derivatives_df_new_chase <- as.data.frame(derivatives_new_chase)
colnames(derivatives_df_new_chase) <- c("derivative_2h_chase1h", "derivative_chase1h_chase2h", "derivative_chase2h_chase8h")

# Calculate mean derivative for new and old data
mean_derivative_new_pulse <- rowMeans(derivatives_df_new_pulse)  # Mean of 30min-1h and 1h-2h for new data
mean_derivative_old_pulse <- rowMeans(derivatives_df_old_pulse)  # Mean of 30min-1h and 1h-2h for old data
mean_derivative_new_chase <- rowMeans(derivatives_df_new_chase)  # Mean of 30min-1h and 1h-2h for new data

# Add the mean derivatives for new and old to the dataframe
derivatives_mean <- data.frame(new_pulse = mean_derivative_new_pulse, old_pulse = mean_derivative_old_pulse, new_chase = mean_derivative_new_chase)

# Merge with RNA half-life data by gene_id
derivatives_mean$gene_id <- merge_df_counts_select$gene_id

plot(derivatives_mean$new_pulse, derivatives_mean$new_chase)
plot(derivatives_mean$old_pulse, derivatives_mean$new_chase)

library(ggplot2)
library(FNN) # For k-nearest neighbors

# --- Calculations for Plot 1 ---

# 1. Calculate correlation and p-value
cor_results1 <- cor.test(derivatives_mean$new_pulse, derivatives_mean$new_chase, method = "pearson")
cor_value1 <- round(cor_results1$estimate, 3)
p_value1 <- signif(cor_results1$p.value, 3)

# 2. Calculate density using k-nearest neighbors
k <- 10 # Number of neighbors
neighbors1 <- get.knn(data = derivatives_mean[, c("new_pulse", "new_chase")], k = k)
derivatives_mean$density1 <- 1 / rowMeans(neighbors1$nn.dist) # Density = inverse of avg distance
derivatives_mean$density1 <- derivatives_mean$density1 / max(derivatives_mean$density1) # Normalize density

# Custom color palette from your template
custom_colors <- c("#54a5de", "white", "#f15389")

# --- Generate Plot 1 ---

FigS2E_plot1 <- ggplot(derivatives_mean, aes(x = new_pulse, y = new_chase, color = density1)) +
  geom_point(size = 2, alpha = 0.7) + # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, se = FALSE) + # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") + # Custom color bar
  labs(
    x = "Loading rate (Pulse)",
    y = "Unloading rate (Chase)",
    title = "Loading Rate Pulse vs Chase"
  ) +
  # Add correlation annotation. You may need to adjust x and y for best placement.
  annotate("text", x = Inf, y = -Inf, label = paste0("r = ", cor_value1, "\nP = ", p_value1),
           hjust = 1.1, vjust = -0.5, size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# To display the plot, run:
# print(plot1)


library(ggplot2)
library(FNN) # For k-nearest neighbors

# --- Calculations for Plot 2 ---

# 1. Calculate correlation and p-value
cor_results2 <- cor.test(derivatives_mean$old_pulse, derivatives_mean$new_chase, method = "pearson")
cor_value2 <- round(cor_results2$estimate, 3)
p_value2 <- signif(cor_results2$p.value, 3)

# 2. Calculate density using k-nearest neighbors
k <- 10 # Number of neighbors
neighbors2 <- get.knn(data = derivatives_mean[, c("old_pulse", "new_chase")], k = k)
derivatives_mean$density2 <- 1 / rowMeans(neighbors2$nn.dist) # Density = inverse of avg distance
derivatives_mean$density2 <- derivatives_mean$density2 / max(derivatives_mean$density2) # Normalize density

# Custom color palette from your template
custom_colors <- c("#54a5de", "white", "#f15389")

# --- Generate Plot 2 ---

FigS2F_plot2 <- ggplot(derivatives_mean, aes(x = old_pulse, y = new_chase, color = density2)) +
  geom_point(size = 2, alpha = 0.7) + # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, se = FALSE) + # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") + # Custom color bar
  labs(
    x = "Unloading rate (Pulse)",
    y = "Unloading rate (Chase)",
    title = "Old Pulse vs. New Chase"
  ) +
  # Add correlation annotation. You may need to adjust x and y for best placement.
  annotate("text", x = Inf, y = -Inf, label = paste0("r = ", cor_value2, "\nP = ", p_value2),
           hjust = 1.1, vjust = -0.5, size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# To display the plot, run:
# print(plot2)

FigS2E_plot1


FigS2F_plot2


write.csv(derivatives_mean, file = "FigS2E.csv", quote = F, row.names = F)
write.csv(derivatives_mean, file = "FigS2F.csv", quote = F, row.names = F)


#FigS2G is in Figure1_and_S1.R


