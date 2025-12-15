library(readr)
load("data/HEK293T_preprocessed.RData")
merge_df_counts_293 <- merge_df_counts[rowMeans(merge_df_counts[,-1])>1000,]

load("data/THP1_preprocessed.RData")
merge_df_counts_THP1 <- merge_df_counts[rowMeans(merge_df_counts[,-1])>1000,]

load("data/MRLs_293.RData")
MRL_new_293 <- MRLs_293[[1]]
MRL_old_293 <- MRLs_293[[2]]

load("data/MRLs_THP1.RData")
MRL_new_THP1 <- MRLs_THP1[[1]]
MRL_old_THP1 <- MRLs_THP1[[2]]

library(pracma)
time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new_293 <- apply(MRL_new_293, 1, calculate_auc)
auc_MRL_old_293 <- apply(MRL_old_293, 1, calculate_auc)
auc_diff_293 <- auc_MRL_new_293 - auc_MRL_old_293

auc_MRL_new_THP1 <- apply(MRL_new_THP1, 1, calculate_auc)
auc_MRL_old_THP1 <- apply(MRL_old_THP1, 1, calculate_auc)
auc_diff_THP1 <- auc_MRL_new_THP1 - auc_MRL_old_THP1

auc_diff_293_data <- data.frame(gene_id = merge_df_counts_293$gene_id, auc_diff_293 = auc_diff_293)
auc_diff_THP1_data <- data.frame(gene_id = merge_df_counts_THP1$gene_id, auc_diff_THP1 = auc_diff_THP1)

auc_diff_293_THP1 <- merge(auc_diff_293_data, auc_diff_THP1_data, by = "gene_id")


library(ggplot2)
library(dplyr)

# Fit a linear model
model_auc <- lm(auc_diff_THP1 ~ auc_diff_293, data = auc_diff_293_THP1)
summary(model_auc)

# Calculate residuals and standard deviation of residuals
auc_diff_293_THP1$residuals <- resid(model_auc)
std_dev <- sd(auc_diff_293_THP1$residuals)

# Identify higher and lower outliers
auc_diff_293_THP1 <- auc_diff_293_THP1 %>%
  mutate(
    rank_residuals = rank(-abs(residuals)),  # Rank by absolute values of residuals
    highlight_upper = rank(-residuals) <= 379,  # Top 50 positive residuals
    highlight_lower = rank(residuals) <= 379    # Top 50 negative residuals
  )

# Create the scatter plot
cor_test_Fig4C <- cor.test(auc_diff_293_THP1$auc_diff_293, auc_diff_293_THP1$auc_diff_THP1)
Fig4C_plot_auc <- ggplot(auc_diff_293_THP1, aes(x = auc_diff_293, y = auc_diff_THP1)) +
  geom_point(aes(color = ifelse(highlight_upper, "#e95966", 
                                ifelse(highlight_lower, "#00bad5", "grey")), alpha = 0.5)) +
  geom_smooth(method = "lm", se = FALSE) +  # Add regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "AUC Diff 293", y = "AUC Diff THP1", title = "AUC Diff 293 vs AUC Diff THP1 with Outliers Highlighted") +
  scale_color_identity() +  # Use color as-is without legend
  coord_cartesian(xlim = c(-10, 5), ylim = c(-10, 5)) +  # Set axis limits
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Add black border
  )

write.csv(auc_diff_293_THP1, file = "Fig4C.csv", quote = F, row.names = F)

load("data/cliffs_delta_293.RData")
load("data/cliffs_deltas_df_THP1.RData")

cliffs_deltas_df_293$Feature <- rownames(cliffs_deltas_df_293)
colnames(cliffs_deltas_df_293) <- c("Feature", "HEK293")
colnames(cliffs_deltas_df_THP1) <- c("Feature", "THP1")

# Merge the two data frames based on the Feature column
merged_df <- merge(cliffs_deltas_df_293, cliffs_deltas_df_THP1, by = "Feature")

library(ggplot2)
library(ggrepel)  # For better label placement
# Create a new column to mark points starting with "m6A"
merged_df$highlight <- grepl("^m6A", merged_df$Feature)

# Create the scatter plot
Fig4D_plot_features <- ggplot(merged_df, aes(x = HEK293, y = THP1, label = Feature)) +
  geom_point(aes(color = highlight), size = 3, alpha = 0.8) +  # Color points conditionally
  geom_text_repel(
    aes(color = ifelse(highlight, "red", "black")),  # Label color conditional on highlight
    size = 3, 
    max.overlaps = 15, 
    segment.color = "grey70"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Diagonal line (y = x)
  labs(
    x = "Cliff's Delta (HEK293)",
    y = "Cliff's Delta (THP1)",
    title = "Comparison of Cliff's Delta for HEK293 and THP1"
  ) +
  scale_color_manual(values = c("FALSE" = "#4ba2dd", "TRUE" = "red"), guide = "none") +  # Set colors
  coord_cartesian(xlim = c(-0.55, 0.4), ylim = c(-0.55, 0.4)) +  # Set x and y axis limits
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

write.csv(merged_df, file = "Fig4D.csv", quote = F, row.names = F)



library(dplyr)
library(VennDiagram)

# Get the top and bottom 1000 genes for auc_diff_293_data
top_1000_293 <- auc_diff_293_data %>%
  arrange(desc(auc_diff_293)) %>%  # Sort descending for top genes
  slice_head(n = 1000) %>%
  pull(gene_id)

bottom_1000_293 <- auc_diff_293_data %>%
  arrange(auc_diff_293) %>%  # Sort ascending for bottom genes
  slice_head(n = 1000) %>%
  pull(gene_id)

# Get the top and bottom 1000 genes for auc_diff_THP1_data
top_1000_THP1 <- auc_diff_THP1_data %>%
  arrange(desc(auc_diff_THP1)) %>%  # Sort descending for top genes
  slice_head(n = 1000) %>%
  pull(gene_id)

bottom_1000_THP1 <- auc_diff_THP1_data %>%
  arrange(auc_diff_THP1) %>%  # Sort ascending for bottom genes
  slice_head(n = 1000) %>%
  pull(gene_id)

# Calculate overlaps
overlap_top <- length(intersect(top_1000_293, top_1000_THP1))
overlap_bottom <- length(intersect(bottom_1000_293, bottom_1000_THP1))


# Create Venn diagrams
# For top 1000 genes
dev.off()
Fig4B_venn_top <- draw.pairwise.venn(
  area1 = length(top_1000_293),
  area2 = length(top_1000_THP1),
  cross.area = overlap_top,
  category = c("HEK293", "THP1"),
  fill = c("#4ba2dd", "#e95966"),
  alpha = 0.5,
  cat.pos = c(-30, 30),
  cat.dist = c(0.05, 0.05),
  main = "Venn Diagram for Top 1000 Genes"
)


# For bottom 1000 genes
dev.off()
Fig4B_venn_bottom <- draw.pairwise.venn(
  area1 = length(bottom_1000_293),
  area2 = length(bottom_1000_THP1),
  cross.area = overlap_bottom,
  category = c("HEK293", "THP1"),
  fill = c("#4ba2dd", "#e95966"),
  alpha = 0.5,
  cat.pos = c(-30, 30),
  cat.dist = c(0.05, 0.05),
  main = "Venn Diagram for Bottom 1000 Genes"
)

dev.off()

write.csv(top_1000_293, file = "Fig4B_293T_top1000.csv", quote = F, row.names = F)
write.csv(bottom_1000_293, file = "Fig4B_293T_bottom1000.csv", quote = F, row.names = F)
write.csv(top_1000_THP1, file = "Fig4B_THP1_top1000.csv", quote = F, row.names = F)
write.csv(bottom_1000_THP1, file = "Fig4B_THP1_bottom1000.csv", quote = F, row.names = F)



#### calculate the derivative to test the difference of loading and unloading rate
derivatives_new_293 <- t(apply(MRL_new_293, 1, function(row) {
  diff(row) / diff(time_points)
}))

derivatives_old_293 <- t(apply(MRL_old_293, 1, function(row) {
  diff(row) / diff(time_points)
}))

# Convert the result to a data frame and add column names
derivatives_df_293 <- as.data.frame(derivatives_new_293)
colnames(derivatives_df_293) <- c("derivative_30min_1h", "derivative_1h_2h")

# Calculate mean derivative for new and old data
mean_derivative_new_293 <- rowMeans(derivatives_new_293)  # Mean of 30min-1h and 1h-2h for new data
mean_derivative_old_293 <- rowMeans(derivatives_old_293)  # Mean of 30min-1h and 1h-2h for old data

# Add the mean derivatives for new and old to the dataframe
derivatives_df_293$mean_derivative_new <- mean_derivative_new_293
derivatives_df_293$mean_derivative_old <- mean_derivative_old_293

derivatives_df_293$gene_id <- merge_df_counts_293$gene_id

derivatives_new_THP1 <- t(apply(MRL_new_THP1, 1, function(row) {
  diff(row) / diff(time_points)
}))

derivatives_old_THP1 <- t(apply(MRL_old_THP1, 1, function(row) {
  diff(row) / diff(time_points)
}))

# Convert the result to a data frame and add column names
derivatives_df_THP1 <- as.data.frame(derivatives_new_THP1)
colnames(derivatives_df_THP1) <- c("derivative_30min_1h", "derivative_1h_2h")

# Calculate mean derivative for new and old data
mean_derivative_new_THP1 <- rowMeans(derivatives_new_THP1)  # Mean of 30min-1h and 1h-2h for new data
mean_derivative_old_THP1 <- rowMeans(derivatives_old_THP1)  # Mean of 30min-1h and 1h-2h for old data

# Add the mean derivatives for new and old to the dataframe
derivatives_df_THP1$mean_derivative_new <- mean_derivative_new_THP1
derivatives_df_THP1$mean_derivative_old <- mean_derivative_old_THP1

derivatives_df_THP1$gene_id <- merge_df_counts_THP1$gene_id


library(ggplot2)
library(gridExtra)

# Provided time points in hours
time_points <- c(0.5, 1, 2)  # Corresponds to 30min, 1h, 2h

# Calculate finite differences (derivative approximation)
# Apply to each row and get the derivative w.r.t. time for new and old data
derivatives_new <- t(apply(MRL_new_THP1, 1, function(row) {
  diff(row) / diff(time_points)
}))

derivatives_old <- t(apply(MRL_old_THP1, 1, function(row) {
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
derivatives_df$gene_id <- merge_df_counts_THP1$gene_id
derivatives_df_half_life <- derivatives_df_THP1

# Filter out rows where RNA_half_life is less than or equal to zero, and remove NA values
filtered_data <- derivatives_df_half_life[complete.cases(derivatives_df_half_life$mean_derivative_new, 
                                                         derivatives_df_half_life$mean_derivative_old,
                                                         derivatives_df_half_life$RNA_half_life)&
                                            (derivatives_df_half_life$mean_derivative_new>-2.5)&
                                            (derivatives_df_half_life$mean_derivative_old>-5),]

# Compute correlation between mean_derivative_new and mean_derivative_old
cor_value_comparison <- cor(filtered_data$mean_derivative_new, filtered_data$mean_derivative_old)

# Create the comparison plot for mean_derivative_new vs mean_derivative_old with a regression line
cor_test_FigS4C <- cor.test(filtered_data$mean_derivative_new, filtered_data$mean_derivative_old)
FigS4C_comparison_plot <- ggplot(filtered_data, aes(x = mean_derivative_new, y = mean_derivative_old)) +
  geom_point(color = "#4ba2dd", size = 3, alpha = 0.1) +  # Add points with transparency
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +  # Add regression line
  labs(x = "Mean Derivative (New)", y = "Mean Derivative (Old)", 
       title = "Comparison of New and Old Mean Derivatives") +  # Labels
  annotate("text", x = Inf, y = Inf, label = paste("Corr =", round(cor_value_comparison, 2)),
           vjust = 2, hjust = 2, size = 5, color = "black") +  # Add correlation label
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

write.csv(filtered_data, file = "FigS5C.csv", quote = F, row.names = F)

# Find common genes between the two data frames
common_genes <- intersect(derivatives_df_293$gene_id, derivatives_df_THP1$gene_id)

# Filter both data frames to include only the common genes
df_293_common <- derivatives_df_293[derivatives_df_293$gene_id %in% common_genes, ]
df_THP1_common <- derivatives_df_THP1[derivatives_df_THP1$gene_id %in% common_genes, ]

# Merge the filtered data frames by 'gene_id'
common_df <- merge(df_293_common, df_THP1_common, by = "gene_id", suffixes = c("_293", "_THP1"))


# Load libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(FNN)  # For nearest neighbor search

# Calculate differences
common_df$diff_new <- common_df$mean_derivative_new_THP1 - common_df$mean_derivative_new_293
common_df$diff_old <- common_df$mean_derivative_old_THP1 - common_df$mean_derivative_old_293

# Calculate correlation and p-value
cor_results <- cor.test(common_df$diff_new, common_df$diff_old, method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10  # Number of neighbors to consider
neighbors <- get.knn(data = common_df[, c("diff_new", "diff_old")], k = k)
common_df$density <- 1 / rowMeans(neighbors$nn.dist)  # Density = inverse of average distance to neighbors
common_df$density <- common_df$density / max(common_df$density)  # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Contour plot with density-colored points and fitting line
Fig4F_loading_vs_unloading <- ggplot(common_df, aes(x = diff_new, y = diff_old, color = density)) +
  geom_point(size = 2, alpha = 0.7) +  # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +  # Fitting line
  geom_hline(yintercept = 0, color = "#555555", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "#555555", linetype = "dashed") +
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") +  # Custom color bar
  labs(x = "Difference in Mean Derivative New (THP1 - 293)", 
       y = "Difference in Mean Derivative Old (THP1 - 293)") +
  annotate("text", x = 1, y = 1.4, label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +
  theme_bw() +  # Use a minimal theme
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )


write.csv(common_df, file = "Fig4F.csv", quote = F, row.names = F)

source("code/polysome_seq_functions.R")
load("data/THP1_preprocessed.RData")
load("data/RNA_features_gene_level.RData")

library(ggplot2)
library(dplyr)
library(patchwork)
library(pracma)
library(gridExtra)
base_theme <- theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


merge_df_counts_select <- merge_df_counts[rowMeans(merge_df_counts[,-1])>1000,]

merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts$gene_id %in% merge_df_counts_select$gene_id,]

merge_df_counts_old_select <- data.frame(gene_id = merge_df_counts_select$gene_id, merge_df_counts_select[,-1] - merge_df_counts_new_select[,-1])

#MRL_total <- calculate_polysome_load(merge_df_counts_select, c("NC_30min", "NC_1h", "NC_2h"))
MRLs <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))
MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]
MRL_total <- MRL_new + MRL_old


#MRL_old <- calculate_polysome_load(merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))
#save(MRL_new, MRL_old, merge_df_counts_select, merge_df_counts_new_select, merge_df_counts_old_select, file = "data/HEK293T_counts_MRL.RData")
# Identifying rows with NA values in each dataframe
na_rows_total <- which(rowSums(is.na(MRL_total)) > 0)
na_rows_new <- which(rowSums(is.na(MRL_new)) > 0)
na_rows_old <- which(rowSums(is.na(MRL_old)) > 0)

# Combining the indices to get a unique set of indices with NA values in any dataframe
na_rows_combined <- unique(c(na_rows_total, na_rows_new, na_rows_old))

# Removing these rows from all three dataframes
if(length(na_rows_combined) != 0) {
  MRL_total <- MRL_total[-na_rows_combined, ]
  MRL_new <- MRL_new[-na_rows_combined, ]
  MRL_old <- MRL_old[-na_rows_combined, ]
  merge_df_counts_select <- merge_df_counts_select[-na_rows_combined, ]
  merge_df_counts_new_select <- merge_df_counts_new_select[-na_rows_combined, ]
  merge_df_counts_old_select <- merge_df_counts_old_select[-na_rows_combined, ]
  
}

MRL_diff <- MRL_new-MRL_old

MRL_diff_output <- MRL_diff
MRL_diff_output$gene_id <- merge_df_counts_select$gene_id


# Calculate the column sums
column_sums <- colSums(MRL_total)

# Calculate the mean of column sums
mean_column_sums <- mean(column_sums)

# Function to scale the dataframes
scale_dataframe <- function(df, column_sums, mean_sum) {
  t(t(df) / column_sums) * mean_sum
}

# Scale MRL_new and MRL_old
MRL_new_scaled <- as.data.frame(scale_dataframe(MRL_new, column_sums, mean_column_sums))
MRL_old_scaled <- as.data.frame(scale_dataframe(MRL_old, column_sums, mean_column_sums))

MRL_new <- MRL_new_scaled
MRL_old <- MRL_old_scaled

colSums(MRL_new)
colSums(MRL_old)

time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new <- apply(MRL_new, 1, calculate_auc)
auc_MRL_old <- apply(MRL_old, 1, calculate_auc)

auc_diff <- auc_MRL_new - auc_MRL_old

MRL_diff$auc_diff <- auc_diff
MRL_diff <- MRL_diff[(!is.infinite(MRL_diff$NC_30min))&(!is.infinite(MRL_diff$NC_1h))&(!is.infinite(MRL_diff$NC_2h)), ]
cor_30min <- cor(MRL_diff$auc_diff, MRL_diff$NC_30min)
cor_1h <- cor(MRL_diff$auc_diff, MRL_diff$NC_1h)
cor_2h <- cor(MRL_diff$auc_diff, MRL_diff$NC_2h)

cor_test_FigS4E_plot1 <- cor.test(MRL_diff$auc_diff, MRL_diff$NC_30min)
cor_test_FigS4E_plot2 <- cor.test(MRL_diff$auc_diff, MRL_diff$NC_1h)
cor_test_FigS4E_plot3 <- cor.test(MRL_diff$auc_diff, MRL_diff$NC_2h)

FigS5E_plot1 <- ggplot(MRL_diff, aes(x = auc_diff, y = NC_30min)) +
  stat_density_2d(geom = "path", color = "#99c3e5") +
  labs(x = "AUC Difference", y = "MRL Difference", title = "NC 30min") +
  base_theme+
  annotate("text", x = Inf, y = Inf, label = sprintf("r = %.2f", cor_30min), 
           hjust = 1.1, vjust = 2, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )


FigS5E_plot2 <- ggplot(MRL_diff, aes(x = auc_diff, y = NC_1h)) +
  stat_density_2d(geom = "path", color = "#99c3e5") +
  labs(x = "AUC Difference", y = "MRL Difference", title = "NC 1h") +
  base_theme +
  annotate("text", x = Inf, y = Inf, label = sprintf("r = %.2f", cor_1h), 
           hjust = 1.1, vjust = 2, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )


FigS5E_plot3 <- ggplot(MRL_diff, aes(x = auc_diff, y = NC_2h)) +
  stat_density_2d(geom = "path", color = "#99c3e5") +
  labs(x = "AUC Difference", y = "MRL Difference", title = "NC 2h") +
  base_theme +
  annotate("text", x = Inf, y = Inf, label = sprintf("r = %.2f", cor_2h), 
           hjust = 1.1, vjust = 2, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

write.csv(MRL_diff, file = "FigS5E.csv", quote = F, row.names = F)



base_theme <- theme_bw() + theme(
  legend.position = "none",  # Remove legend
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
  panel.grid.major = element_blank(),  # Remove major grid lines
  panel.grid.minor = element_blank(),  # Remove minor grid lines
  panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
)

source("code/polysome_seq_functions.R")
load("data/THP1_preprocessed.RData")
merge_df_counts_select <- merge_df_counts[merge_df_counts$gene_id %in% genes_logCPM[genes_logCPM$logCPM>1,"gene_id"],]
merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts$gene_id %in% genes_logCPM[genes_logCPM$logCPM>1,"gene_id"],]

merge_df_counts_old_select <- data.frame(gene_id = merge_df_counts_select$gene_id, merge_df_counts_select[,-1] - merge_df_counts_new_select[,-1])

#MRL_total <- calculate_polysome_load(merge_df_counts_select, c("NC_30min", "NC_1h", "NC_2h"))
MRLs <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))

MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]
MRL_total <- MRL_new + MRL_old
#MRL_old <- calculate_polysome_load(merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))

# Identifying rows with NA values in each dataframe
na_rows_total <- which(rowSums(is.na(MRL_total)) > 0)
na_rows_new <- which(rowSums(is.na(MRL_new)) > 0)
na_rows_old <- which(rowSums(is.na(MRL_old)) > 0)

# Combining the indices to get a unique set of indices with NA values in any dataframe
na_rows_combined <- unique(c(na_rows_total, na_rows_new, na_rows_old))

# Removing these rows from all three dataframes
if(!isempty(na_rows_combined)) {
  MRL_total <- MRL_total[-na_rows_combined, ]
  MRL_new <- MRL_new[-na_rows_combined, ]
  MRL_old <- MRL_old[-na_rows_combined, ]
  merge_df_counts_select <- merge_df_counts_select[-na_rows_combined, ]
  merge_df_counts_new_select <- merge_df_counts_new_select[-na_rows_combined, ]
  merge_df_counts_old_select <- merge_df_counts_old_select[-na_rows_combined, ]
}

# Calculate the column sums
column_sums <- colSums(MRL_total)

# Calculate the mean of column sums
mean_column_sums <- mean(column_sums)

# Function to scale the dataframes
scale_dataframe <- function(df, column_sums, mean_sum) {
  t(t(df) / column_sums) * mean_sum
}

# Scale MRL_new and MRL_old
MRL_new_scaled <- as.data.frame(scale_dataframe(MRL_new, column_sums, mean_column_sums))
MRL_old_scaled <- as.data.frame(scale_dataframe(MRL_old, column_sums, mean_column_sums))

MRL_new <- MRL_new_scaled
MRL_old <- MRL_old_scaled

colSums(MRL_new)
colSums(MRL_old)



library(pracma)
time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new <- apply(MRL_new, 1, calculate_auc)
auc_MRL_old <- apply(MRL_old, 1, calculate_auc)

auc_diff <- auc_MRL_new - auc_MRL_old

auc_diff_ranked <- rank(auc_diff, ties.method = "first")
high_coupling_indices <- which(auc_diff_ranked > length(auc_diff_ranked) - 1000)
low_coupling_indices <- which(auc_diff_ranked <= 1000)

# Load the ggplot2 package
library(ggplot2)


# Calculate the values for the top and bottom 2000 indices
high_coupling_value <- min(auc_diff[high_coupling_indices])
low_coupling_value <- max(auc_diff[low_coupling_indices])
high_coupling_genes <- merge_df_counts_select[high_coupling_indices,1]
low_coupling_genes <- merge_df_counts_select[low_coupling_indices,1]


#write.table(high_coupling_genes, file = "./high_coupling_genes_THP1.csv", quote = FALSE, sep = ',')
#write.table(low_coupling_genes, file = "./low_coupling_genes_THP1.csv", quote = FALSE, sep = ',')
# Compute the density
d <- density(auc_diff)

# Create a data frame from the density object
density_data <- data.frame(x = d$x, y = d$y)

# Split the data based on the thresholds
low_segment <- subset(density_data, x <= low_coupling_value)
mid_segment <- subset(density_data, x > low_coupling_value & x < high_coupling_value)
high_segment <- subset(density_data, x >= high_coupling_value)

# Calculate y-values for the vertical lines at the thresholds
low_density_value <- d$y[which.min(abs(d$x - low_coupling_value))]
high_density_value <- d$y[which.min(abs(d$x - high_coupling_value))]

# Create the density plot with adjusted vertical lines
FigS5B <- ggplot() +
  geom_ribbon(data = low_segment, aes(x = x, ymin = 0, ymax = y), fill = "#95c08b", alpha = 0.5) +
  geom_ribbon(data = mid_segment, aes(x = x, ymin = 0, ymax = y), fill = "#acb2b5", alpha = 0.5) +
  geom_ribbon(data = high_segment, aes(x = x, ymin = 0, ymax = y), fill = "#8887b9", alpha = 0.5) +
  geom_line(data = low_segment, aes(x = x, y = y), color = "#95c08b", size = 1.5) +
  geom_line(data = mid_segment, aes(x = x, y = y), color = "#acb2b5", size = 1.5) +
  geom_line(data = high_segment, aes(x = x, y = y), color = "#8887b9", size = 1.5) +
  geom_segment(aes(x = low_coupling_value, xend = low_coupling_value, y = 0, yend = low_density_value), 
               color="#95c08b", linetype="dashed", size=1) +
  geom_segment(aes(x = high_coupling_value, xend = high_coupling_value, y = 0, yend = high_density_value), 
               color="#8887b9", linetype="dashed", size=1) +
  labs(title="separate high and low by coupling AUC difference", x="AUC diff", y="Density") +
  theme_classic()


write.csv(density_data, file = "FigS5B.csv", quote = F, row.names = F)





# -----------------------------------------------------------------------------
# Section 1: Load Libraries
# -----------------------------------------------------------------------------
# Ensure all required packages are installed and loaded.
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(Mfuzz)
library(Biobase)
library(tibble)
library(tools)

# -----------------------------------------------------------------------------
# Section 2: Define Reusable Functions
# -----------------------------------------------------------------------------

# --- TPM Normalization Function ---
# Note: This function assumes a fixed gene length of 1kb as a placeholder.
# For accurate TPM, you should use actual gene lengths.
calculate_tpm <- function(counts) {
  gene_lengths_kb <- rep(1, nrow(counts)) # Replace with actual gene lengths if available
  rpk <- counts / gene_lengths_kb
  scaling_factors <- colSums(rpk, na.rm = TRUE)
  tpm <- t(t(rpk) / scaling_factors * 1e6)
  return(tpm)
}

# --- Fold Change and P-value Calculation Function ---
calculate_foldchange_pvalue <- function(df, cols_time1, cols_time2) {
  # Ensure columns are numeric before calculation
  df[cols_time1] <- sapply(df[cols_time1], as.numeric)
  df[cols_time2] <- sapply(df[cols_time2], as.numeric)
  
  fold_changes <- rowMeans(log2(df[cols_time1] + 1), na.rm = TRUE) - rowMeans(log2(df[cols_time2] + 1), na.rm = TRUE)
  
  p_values <- apply(df[, c(cols_time1, cols_time2)], 1, function(row) {
    # Check for sufficient data points to perform a t-test
    if (sum(!is.na(row[cols_time1])) < 2 | sum(!is.na(row[cols_time2])) < 2) {
      return(NA)
    }
    # Using a standard t-test as paired test might not be appropriate if samples are independent
    t_test_result <- t.test(as.numeric(row[cols_time1]), as.numeric(row[cols_time2]))
    return(t_test_result$p.value)
  })
  
  return(data.frame(fold_changes, p_values))
}


# --- Gene ID Truncation Function for GSEA plots ---
truncate_geneID <- function(gene_string) {
  genes <- strsplit(gene_string, "/")[[1]]
  # Keep only the first 5 genes for cleaner plot labels
  if (length(genes) > 5) {
    return(paste(genes[1:5], collapse = "/"))
  } else {
    return(gene_string)
  }
}


# -----------------------------------------------------------------------------
# Section 3: Data Loading and Pre-processing
# -----------------------------------------------------------------------------
# This section prepares the dataframes required for all three figures.

# --- Load and Normalize mRNA Data ---
# !! IMPORTANT: Update this path to your file location.
THP1_LPS_Total_RNA_seq <- read_delim("data/THP1_LPS_Total_RNA_seq.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

countData <- THP1_LPS_Total_RNA_seq[, c(8:19)]
THP1_mRNA <- cbind(THP1_LPS_Total_RNA_seq[, c(1, 7)], countData)

normalized_THP1_mRNA <- THP1_mRNA %>%
  dplyr::select(-Geneid, -gene_name) %>%
  as.matrix() %>%
  calculate_tpm() %>%
  as.data.frame()

normalized_THP1_mRNA <- cbind(THP1_mRNA[, 1:2], normalized_THP1_mRNA)
normalized_THP1_mRNA$Geneid <- sub("\\..*", "", normalized_THP1_mRNA$Geneid)
colnames(normalized_THP1_mRNA)[colnames(normalized_THP1_mRNA) == 'Geneid'] <- "gene_id"


# --- Load and Normalize Protein Data ---
# !! IMPORTANT: Update this path to your file location.
load("data/THP1_data_abundance.RData") # This loads the 'THP1_data_ensembl' object

normalized_THP1_data_ensembl <- THP1_data_ensembl %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-accession, -gene_id, -MW) %>%
  as.matrix() %>%
  calculate_tpm() %>%
  as.data.frame()

normalized_THP1_data_ensembl <- cbind(THP1_data_ensembl[, 1:3], normalized_THP1_data_ensembl)


# --- Load MRL Difference Data ---
# !! IMPORTANT: Update this path to your file location.
THP1_gene_level_MRL_diff <- read_csv("data/THP1_gene_level_MRL_diff.csv")
THP1_gene_level_MRL_diff$gene_id <- sub("\\..*", "", THP1_gene_level_MRL_diff$gene_id)
THP1_gene_level_MRL_diff <- data.frame(THP1_gene_level_MRL_diff$gene_id, THP1_gene_level_MRL_diff$NC_2h)
colnames(THP1_gene_level_MRL_diff) <- c("gene_id", "MRL_diff_2h")


# --- Merge mRNA, Protein, and MRL Data ---
THP1_data_protein_mRNA_compare <- merge(normalized_THP1_data_ensembl, normalized_THP1_mRNA, by = "gene_id")
THP1_data_protein_mRNA_compare <- merge(THP1_data_protein_mRNA_compare, THP1_gene_level_MRL_diff, by = "gene_id")
THP1_data_protein_mRNA_compare_clear <- na.omit(THP1_data_protein_mRNA_compare)


# --- Calculate Fold Changes and P-values for All Time Points ---
# This creates the 'results' dataframe needed for FigS5A

# Define column groups
protein_cols_NC <- c("pro_NC-1", "pro_NC-2", "pro_NC-3")
protein_cols_2h <- c("pro_LPS2h-1", "pro_LPS2h-2", "pro_LPS2h-3")
protein_cols_12h <- c("pro_LPS12h-1", "pro_LPS12h-2", "pro_LPS12h-3")
protein_cols_24h <- c("pro_LPS24h-1", "pro_LPS24h-2", "pro_LPS24h-3")

mRNA_cols_NC <- c("NC-1", "NC-2", "NC-3")
mRNA_cols_2h <- c("LPS2h-1", "LPS2h-2", "LPS2h-3")
mRNA_cols_12h <- c("LPS12h-1", "LPS12h-2", "LPS12h-3")
mRNA_cols_24h <- c("LPS24h-1", "LPS24h-2", "LPS24h-3")

# Perform calculations
protein_foldchange_pvalue_2h <- calculate_foldchange_pvalue(THP1_data_protein_mRNA_compare_clear, protein_cols_2h, protein_cols_NC)
mRNA_foldchange_pvalue_2h <- calculate_foldchange_pvalue(THP1_data_protein_mRNA_compare_clear, mRNA_cols_2h, mRNA_cols_NC)

protein_foldchange_pvalue_12h <- calculate_foldchange_pvalue(THP1_data_protein_mRNA_compare_clear, protein_cols_12h, protein_cols_NC)
mRNA_foldchange_pvalue_12h <- calculate_foldchange_pvalue(THP1_data_protein_mRNA_compare_clear, mRNA_cols_12h, mRNA_cols_NC)

protein_foldchange_pvalue_24h <- calculate_foldchange_pvalue(THP1_data_protein_mRNA_compare_clear, protein_cols_24h, protein_cols_NC)
mRNA_foldchange_pvalue_24h <- calculate_foldchange_pvalue(THP1_data_protein_mRNA_compare_clear, mRNA_cols_24h, mRNA_cols_NC)

# Combine all results into a single dataframe
results <- data.frame(
  gene_id = THP1_data_protein_mRNA_compare_clear$gene_id,
  gene_name = THP1_data_protein_mRNA_compare_clear$gene_name,
  protein_fold_change_2h = protein_foldchange_pvalue_2h$fold_changes,
  protein_p_value_2h = protein_foldchange_pvalue_2h$p_values,
  mRNA_fold_change_2h = mRNA_foldchange_pvalue_2h$fold_changes,
  mRNA_p_value_2h = mRNA_foldchange_pvalue_2h$p_values,
  protein_fold_change_12h = protein_foldchange_pvalue_12h$fold_changes,
  protein_p_value_12h = protein_foldchange_pvalue_12h$p_values,
  mRNA_fold_change_12h = mRNA_foldchange_pvalue_12h$fold_changes,
  mRNA_p_value_12h = mRNA_foldchange_pvalue_12h$p_values,
  protein_fold_change_24h = protein_foldchange_pvalue_24h$fold_changes,
  protein_p_value_24h = protein_foldchange_pvalue_24h$p_values,
  mRNA_fold_change_24h = mRNA_foldchange_pvalue_24h$fold_changes,
  mRNA_p_value_24h = mRNA_foldchange_pvalue_24h$p_values,
  MRL_diff = THP1_data_protein_mRNA_compare_clear$MRL_diff_2h
)

# Replace NaN with NA for cleaner plotting
results[is.na(results)] <- NA


# -----------------------------------------------------------------------------
# Section 4: Generate Figure S5A (All 6 Volcano Plots)
# -----------------------------------------------------------------------------

# --- Plot 1: mRNA changes at 2h ---
results$significance_2 <- ifelse(results$mRNA_p_value_2h < 0.05 & abs(results$mRNA_fold_change_2h) > 1,
                                 ifelse(results$mRNA_fold_change_2h > 1, "Increase", "Decrease"), "Not Significant")

FigS6A_p_mRNA_2h <- ggplot(results, aes(x = mRNA_fold_change_2h, y = -log10(mRNA_p_value_2h), color = significance_2, label = ifelse(significance_2 != "Not Significant", gene_name, ""))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#bd0026", "Decrease" = "#253494", "Not Significant" = "grey")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(title = "LPS 2 h mRNA change", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel()

# --- Plot 2: mRNA changes at 12h ---
results$significance_12 <- ifelse(results$mRNA_p_value_12h < 0.05 & abs(results$mRNA_fold_change_12h) > 1,
                                  ifelse(results$mRNA_fold_change_12h > 1, "Increase", "Decrease"), "Not Significant")

FigS6A_p_mRNA_12h <- ggplot(results, aes(x = mRNA_fold_change_12h, y = -log10(mRNA_p_value_12h), color = significance_12, label = ifelse(significance_12 != "Not Significant", gene_name, ""))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#bd0026", "Decrease" = "#253494", "Not Significant" = "grey")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(title = "LPS 12 h mRNA change", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel()

# --- Plot 3: mRNA changes at 24h ---
results$significance_24 <- ifelse(results$mRNA_p_value_24h < 0.05 & abs(results$mRNA_fold_change_24h) > 1,
                                  ifelse(results$mRNA_fold_change_24h > 1, "Increase", "Decrease"), "Not Significant")

FigS6A_p_mRNA_24h <- ggplot(results, aes(x = mRNA_fold_change_24h, y = -log10(mRNA_p_value_24h), color = significance_24, label = ifelse(significance_24 != "Not Significant", gene_name, ""))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#bd0026", "Decrease" = "#253494", "Not Significant" = "grey")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(title = "LPS 24 h mRNA change", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel()

# --- Plot 4: Protein changes at 2h ---
results$significance_2p <- ifelse(results$protein_p_value_2h < 0.05 & abs(results$protein_fold_change_2h) > 1,
                                  ifelse(results$protein_fold_change_2h > 1, "Increase", "Decrease"), "Not Significant")

FigS6A_p_protein_2h <- ggplot(results, aes(x = protein_fold_change_2h, y = -log10(protein_p_value_2h), color = significance_2p, label = ifelse(significance_2p != "Not Significant", gene_name, ""))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#bd0026", "Decrease" = "#253494", "Not Significant" = "grey")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(title = "LPS 2 h protein change", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel()

# --- Plot 5: Protein changes at 12h ---
results$significance_12p <- ifelse(results$protein_p_value_12h < 0.05 & abs(results$protein_fold_change_12h) > 1,
                                   ifelse(results$protein_fold_change_12h > 1, "Increase", "Decrease"), "Not Significant")

FigS6A_p_protein_12h <- ggplot(results, aes(x = protein_fold_change_12h, y = -log10(protein_p_value_12h), color = significance_12p, label = ifelse(significance_12p != "Not Significant", gene_name, ""))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#bd0026", "Decrease" = "#253494", "Not Significant" = "grey")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(title = "LPS 12 h protein change", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel()

# --- Plot 6: Protein changes at 24h ---
results$significance_24p <- ifelse(results$protein_p_value_24h < 0.05 & abs(results$protein_fold_change_24h) > 1,
                                   ifelse(results$protein_fold_change_24h > 1, "Increase", "Decrease"), "Not Significant")

FigS6A_p_protein_24h <- ggplot(results, aes(x = protein_fold_change_24h, y = -log10(protein_p_value_24h), color = significance_24p, label = ifelse(significance_24p != "Not Significant", gene_name, ""))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Increase" = "#bd0026", "Decrease" = "#253494", "Not Significant" = "grey")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(title = "LPS 24 h protein change", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel()

# --- Display FigS5A plots ---
print(FigS6A_p_mRNA_2h)
print(FigS6A_p_mRNA_12h)
print(FigS6A_p_mRNA_24h)
print(FigS6A_p_protein_2h)
print(FigS6A_p_protein_12h)
print(FigS6A_p_protein_24h)


write.csv(results, file = "FigS6A.csv", quote = F, row.names = F)

# -----------------------------------------------------------------------------
# Section 5: Generate Figure 4I (GSEA Plot)
# -----------------------------------------------------------------------------

# --- Prepare data for GSEA ---
# !! IMPORTANT: Update this path to your file location.
load("data/MRL_LPS_diff.RData") # This loads the 'MRL_LPS_diff' object
MRL_LPS_diff$gene_id <- sapply(strsplit(MRL_LPS_diff$gene_id, "\\."), `[`, 1)

# Merge with previously cleared data to get all necessary columns
MRL_LPS_diff <- merge(MRL_LPS_diff, THP1_data_protein_mRNA_compare_clear, by = "gene_id")

MRL_LPS_diff <- MRL_LPS_diff %>%
  mutate(MRL_diff_increase = LPS_2h - NC_2h)

# --- Create ranked gene list ---
ranked_metric <- MRL_LPS_diff$MRL_diff_increase
names(ranked_metric) <- MRL_LPS_diff$gene_id
ranked_metric <- sort(ranked_metric, decreasing = TRUE)
ranked_metric <- ranked_metric[!is.na(ranked_metric)] # Remove NAs

# --- Map ENSEMBL IDs to Gene Symbols ---
ensembl_ids <- sub("\\..*", "", names(ranked_metric))
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first"
)
valid_indices <- !is.na(gene_symbols)
ranked_metric <- ranked_metric[valid_indices]
names(ranked_metric) <- gene_symbols[valid_indices]

# --- Load Gene Set Data and Run GSEA ---
# !! IMPORTANT: Update this path to your file location.
gmt_file_c2 <- "data/GSEA_sets/c2.all.v2024.1.Hs.symbols.gmt"
gene_sets_c2 <- read.gmt(gmt_file_c2)

gsea_result_c2 <- GSEA(ranked_metric,
                       TERM2GENE = gene_sets_c2,
                       pvalueCutoff = 0.05,
                       verbose = FALSE
)
result_c2 <- as.data.frame(gsea_result_c2@result)


# --- Process GSEA results for plotting ---
top_10_positive_c2 <- result_c2 %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice(1:10)

bottom_10_negative_c2 <- result_c2 %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice(1:10)

top_bottom_10_c2 <- bind_rows(top_10_positive_c2, bottom_10_negative_c2)

# Format the c2 Gene Set Descriptions
top_bottom_10_c2$Description <- top_bottom_10_c2$Description %>%
  gsub("_", " ", .) %>%
  toTitleCase()

# Add Cluster and abs_NES columns
top_bottom_10_c2 <- top_bottom_10_c2 %>%
  mutate(Cluster = ifelse(NES > 0, "Fast", "Slow"),
         abs_NES = abs(NES))

top_bottom_10_c2 <- read_csv("data/top_bottom_10_c2_LPS_activation_mod.csv")
top_bottom_10_c2$geneID <- sapply(top_bottom_10_c2$geneID, truncate_geneID)
top_bottom_10_c2 <- top_bottom_10_c2 %>%
  arrange(Cluster, desc(abs_NES)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# --- Generate the plot for Fig4I ---
Fig4I_p2 <- ggplot(top_bottom_10_c2, aes(x = abs_NES, y = Description, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(x = 0.1, y = Description, label = Description), size = 3.5, hjust = 0) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.5),
    axis.text.x = element_text(colour = 'black', size = 10),
    axis.ticks.x = element_line(colour = 'black'),
    axis.title.x = element_text(colour = 'black', size = 12),
    legend.position = "none"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("Slow" = "#8887b9", "Fast" = "#95c08b")) +
  geom_text(
    data = top_bottom_10_c2,
    aes(x = 0.1, y = Description, label = geneID, color = Cluster),
    size = 4, fontface = 'italic', hjust = 0, vjust = 2.3
  ) +
  scale_color_manual(values = c("Slow" = "#8887b9", "Fast" = "#95c08b")) +
  scale_y_discrete(expand = c(0.1, 0)) +
  ggtitle("Curated Gene Enrichment in LPS activation")

# --- Display Fig4I plot ---
print(Fig4I_p2)

write.csv(top_bottom_10_c2, file = "Fig4I.csv", quote = F, row.names = F)

# -----------------------------------------------------------------------------
# Section 6: Generate Figure 4K (Clustered Heatmap)
# -----------------------------------------------------------------------------

# --- Prepare data for clustering ---
# Rename mRNA columns to avoid conflicts
MRL_LPS_diff <- MRL_LPS_diff %>%
  rename_with(~paste0("mRNA_", .x), .cols = all_of(c(mRNA_cols_NC, mRNA_cols_2h, mRNA_cols_12h, mRNA_cols_24h)))

# --- Calculate protein/mRNA ratios and normalize ---
pro_mRNA_ratio_df <- MRL_LPS_diff %>%
  rowwise() %>%
  mutate(
    pro_NC_mean = mean(c_across(starts_with("pro_NC")), na.rm = TRUE),
    pro_LPS2h_mean = mean(c_across(starts_with("pro_LPS2h")), na.rm = TRUE),
    pro_LPS12h_mean = mean(c_across(starts_with("pro_LPS12h")), na.rm = TRUE),
    pro_LPS24h_mean = mean(c_across(starts_with("pro_LPS24h")), na.rm = TRUE),
    mRNA_NC_mean = mean(c_across(starts_with("mRNA_NC")), na.rm = TRUE),
    mRNA_LPS2h_mean = mean(c_across(starts_with("mRNA_LPS2h")), na.rm = TRUE),
    mRNA_LPS12h_mean = mean(c_across(starts_with("mRNA_LPS12h")), na.rm = TRUE),
    mRNA_LPS24h_mean = mean(c_across(starts_with("mRNA_LPS24h")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    ratio_NC = pro_NC_mean / mRNA_NC_mean,
    ratio_LPS2h = pro_LPS2h_mean / mRNA_LPS2h_mean,
    ratio_LPS12h = pro_LPS12h_mean / mRNA_LPS12h_mean,
    ratio_LPS24h = pro_LPS24h_mean / mRNA_LPS24h_mean
  ) %>%
  dplyr::select(gene_id, ratio_NC, ratio_LPS2h, ratio_LPS12h, ratio_LPS24h)

# --- Normalize ratios row-wise ---
pro_mRNA_ratio_normalized <- pro_mRNA_ratio_df %>%
  rowwise() %>%
  mutate(
    across(starts_with("ratio_"), ~ (. - min(c_across(starts_with("ratio_")), na.rm = TRUE)) /
             (max(c_across(starts_with("ratio_")), na.rm = TRUE) - min(c_across(starts_with("ratio_")), na.rm = TRUE)))
  ) %>%
  ungroup()

# --- Prepare matrix for clustering ---
heatmap_matrix <- pro_mRNA_ratio_normalized %>%
  dplyr::select(-gene_id) %>%
  as.matrix()
rownames(heatmap_matrix) <- pro_mRNA_ratio_normalized$gene_id

heatmap_matrix_clean <- heatmap_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  group_by(gene_id) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()

# Remove any rows with NA/NaN values
valid_rows <- complete.cases(heatmap_matrix_clean) & !apply(heatmap_matrix_clean, 1, function(x) any(is.nan(x)))
heatmap_matrix_clean <- heatmap_matrix_clean[valid_rows, ]

# --- Run Mfuzz clustering ---
eset <- new('ExpressionSet', exprs = heatmap_matrix_clean)
m_value <- mestimate(eset)
cl <- mfuzz(eset, c = 2, m = m_value)

# --- Prepare data for heatmap plotting ---
cluster_assignments <- cl$cluster
pro_mRNA_ratio_normalized_clustered <- as.data.frame(heatmap_matrix_clean) %>%
  mutate(Cluster = factor(cluster_assignments)) %>%
  arrange(Cluster)


heatmap_matrix_sorted <- pro_mRNA_ratio_normalized_clustered %>%
  dplyr::select(-Cluster) %>%
  as.matrix()

# --- Define row annotation for the heatmap ---
row_annotation <- rowAnnotation(
  Cluster = pro_mRNA_ratio_normalized_clustered$Cluster,
  col = list(Cluster = c("1" = "#95c08b", "2" = "#8887b9"))
)

# --- Generate the plot for Fig4K ---
Fig4K_heatmap_protein_to_mRNA <- Heatmap(heatmap_matrix_sorted,
                                         name = "Protein/mRNA Ratio",
                                         col = colorRamp2(c(0, 0.5, 1), c("white", "orange", "darkred")),
                                         cluster_rows = FALSE,
                                         cluster_columns = FALSE,
                                         show_row_names = FALSE,
                                         row_split = pro_mRNA_ratio_normalized_clustered$Cluster,
                                         left_annotation = row_annotation,
                                         column_title = "Conditions",
                                         row_title = "Genes"
)

# --- Display Fig4K plot ---
draw(Fig4K_heatmap_protein_to_mRNA)


write.csv(heatmap_matrix_sorted, file = "Fig4K.csv", quote = F, row.names = F)


result_c2 <- read_csv("data/Fig4K_curated_gene_GSEA_table_mod.csv")



# --- Process GSEA results for plotting ---
top_10_positive_c2 <- result_c2 %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  dplyr::slice(1:10)

bottom_10_negative_c2 <- result_c2 %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  dplyr::slice(1:10)

top_bottom_10_c2 <- bind_rows(top_10_positive_c2, bottom_10_negative_c2)

# Format the c2 Gene Set Descriptions
top_bottom_10_c2$Description <- top_bottom_10_c2$Description %>%
  gsub("_", " ", .) %>%
  toTitleCase()

# Add Cluster and abs_NES columns
top_bottom_10_c2 <- top_bottom_10_c2 %>%
  mutate(Cluster = ifelse(NES > 0, "Fast", "Slow"),
         abs_NES = abs(NES))

#top_bottom_10_c2 <- read_csv("data/top_bottom_10_c2_LPS_activation_mod.csv")
top_bottom_10_c2$geneID <- sapply(top_bottom_10_c2$core_enrichment, truncate_geneID)
top_bottom_10_c2 <- top_bottom_10_c2 %>%
  arrange(Cluster, desc(abs_NES)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# --- Generate the plot for Fig4I ---
FigS6C_p2 <- ggplot(top_bottom_10_c2, aes(x = abs_NES, y = Description, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(x = 0.1, y = Description, label = Description), size = 3.5, hjust = 0) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.5),
    axis.text.x = element_text(colour = 'black', size = 10),
    axis.ticks.x = element_line(colour = 'black'),
    axis.title.x = element_text(colour = 'black', size = 12),
    legend.position = "none"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("Fast" = "#8887b9", "Slow" = "#95c08b")) +
  geom_text(
    data = top_bottom_10_c2,
    aes(x = 0.1, y = Description, label = geneID, color = Cluster),
    size = 4, fontface = 'italic', hjust = 0, vjust = 2.3
  ) +
  scale_color_manual(values = c("Fast" = "#8887b9", "Slow" = "#95c08b")) +
  scale_y_discrete(expand = c(0.1, 0)) +
  ggtitle("curated Gene Enrichment in LPS activation")


write.csv(result_c2, file = "FigS6C.csv", quote = F, row.names = F)



# -----------------------------------------------------------------------------
# Section 1: Load Libraries
# -----------------------------------------------------------------------------
# Ensure all required packages are installed and loaded.
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr) # For pivot_longer
library(Mfuzz)
library(Biobase)
library(tibble)

# -----------------------------------------------------------------------------
# Section 2: Define Reusable Functions
# -----------------------------------------------------------------------------

# --- TPM Normalization Function ---
# Note: This function assumes a fixed gene length of 1kb as a placeholder.
# For accurate TPM, you should use actual gene lengths.
calculate_tpm <- function(counts) {
  gene_lengths_kb <- rep(1, nrow(counts)) # Replace with actual gene lengths if available
  rpk <- counts / gene_lengths_kb
  scaling_factors <- colSums(rpk, na.rm = TRUE)
  tpm <- t(t(rpk) / scaling_factors * 1e6)
  return(tpm)
}

# -----------------------------------------------------------------------------
# Section 3: Data Loading and Pre-processing
# -----------------------------------------------------------------------------
# This section prepares the dataframes required for both figures.

# --- Load and Normalize mRNA Data ---
# !! IMPORTANT: Update this path to your file location.
THP1_LPS_Total_RNA_seq <- read_delim("data/THP1_LPS_Total_RNA_seq.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

countData <- THP1_LPS_Total_RNA_seq[, c(8:19)]
THP1_mRNA <- cbind(THP1_LPS_Total_RNA_seq[, c(1, 7)], countData)

normalized_THP1_mRNA <- THP1_mRNA %>%
  dplyr::select(-Geneid, -gene_name) %>%
  as.matrix() %>%
  calculate_tpm() %>%
  as.data.frame()

normalized_THP1_mRNA <- cbind(THP1_mRNA[, 1:2], normalized_THP1_mRNA)
normalized_THP1_mRNA$Geneid <- sub("\\..*", "", normalized_THP1_mRNA$Geneid)
colnames(normalized_THP1_mRNA)[colnames(normalized_THP1_mRNA) == 'Geneid'] <- "gene_id"

# --- Load and Normalize Protein Data ---
# !! IMPORTANT: Update this path to your file location.
load("data/THP1_data_abundance.RData") # This loads the 'THP1_data_ensembl' object

normalized_THP1_data_ensembl <- THP1_data_ensembl %>%
  replace(is.na(.), 0) %>%
  dplyr::select(-accession, -gene_id, -MW) %>%
  as.matrix() %>%
  calculate_tpm() %>%
  as.data.frame()

normalized_THP1_data_ensembl <- cbind(THP1_data_ensembl[, 1:3], normalized_THP1_data_ensembl)

# --- Load MRL Difference Data ---
# !! IMPORTANT: Update this path to your file location.
THP1_gene_level_MRL_diff <- read_csv("data/THP1_gene_level_MRL_diff.csv")
THP1_gene_level_MRL_diff$gene_id <- sub("\\..*", "", THP1_gene_level_MRL_diff$gene_id)
THP1_gene_level_MRL_diff <- data.frame(THP1_gene_level_MRL_diff$gene_id, THP1_gene_level_MRL_diff$NC_2h)
colnames(THP1_gene_level_MRL_diff) <- c("gene_id", "MRL_diff_2h")

# --- Merge Data for Initial Analysis ---
THP1_data_protein_mRNA_compare <- merge(normalized_THP1_data_ensembl, normalized_THP1_mRNA, by = "gene_id")
THP1_data_protein_mRNA_compare <- merge(THP1_data_protein_mRNA_compare, THP1_gene_level_MRL_diff, by = "gene_id")
THP1_data_protein_mRNA_compare_clear <- na.omit(THP1_data_protein_mRNA_compare)

# --- Load and Process MRL_LPS_diff Data ---
# !! IMPORTANT: Update this path to your file location.
load("data/MRL_LPS_diff.RData") # This loads the 'MRL_LPS_diff' object
MRL_LPS_diff$gene_id <- sapply(strsplit(MRL_LPS_diff$gene_id, "\\."), `[`, 1)
MRL_LPS_diff <- merge(MRL_LPS_diff, THP1_data_protein_mRNA_compare_clear, by = "gene_id")


# -----------------------------------------------------------------------------
# Section 4: Generate Figure 4J (Protein Ratio Boxplot)
# -----------------------------------------------------------------------------

# --- Calculate Ratios ---
MRL_LPS_diff_fig4j <- MRL_LPS_diff %>%
  mutate(
    protein_ratio = (`pro_LPS12h-1` / `pro_NC-1`)
  )

# --- Prepare Data for Plotting ---
upper_data <- data.frame(
  protein_ratio = MRL_LPS_diff_fig4j$protein_ratio[MRL_LPS_diff_fig4j$highlight_upper1 == TRUE],
  group = "Upper"
)
lower_data <- data.frame(
  protein_ratio = MRL_LPS_diff_fig4j$protein_ratio[MRL_LPS_diff_fig4j$highlight_lower1 == TRUE],
  group = "Lower"
)

plot_data_fig4j <- rbind(upper_data, lower_data)
plot_data_fig4j <- plot_data_fig4j[is.finite(plot_data_fig4j$protein_ratio), ]

# --- Dynamically set plot limits ---
q_upper <- quantile(plot_data_fig4j$protein_ratio[plot_data_fig4j$group == "Upper"], probs = c(0.25, 0.75), na.rm = TRUE)
iqr_upper <- q_upper[2] - q_upper[1]
upper_whisker_upper <- q_upper[2] + 1.5 * iqr_upper

q_lower <- quantile(plot_data_fig4j$protein_ratio[plot_data_fig4j$group == "Lower"], probs = c(0.25, 0.75), na.rm = TRUE)
iqr_lower <- q_lower[2] - q_lower[1]
upper_whisker_lower <- q_lower[2] + 1.5 * iqr_lower

y_axis_max <- max(upper_whisker_upper, upper_whisker_lower)
label_y_position <- y_axis_max * 0.9

# --- Generate the Plot for Fig4J ---





# -----------------------------------------------------------------------------
# Section 5: Generate Figure S5C (Grouped MRL Difference Boxplot)
# -----------------------------------------------------------------------------

# --- Step 5.1: Prepare Data for Clustering ---
MRL_LPS_diff_2 <- MRL_LPS_diff %>%
  rowwise() %>%
  mutate(
    pro_NC_mean = mean(c_across(starts_with("pro_NC")), na.rm = TRUE),
    pro_LPS2h_mean = mean(c_across(starts_with("pro_LPS2h")), na.rm = TRUE),
    pro_LPS12h_mean = mean(c_across(starts_with("pro_LPS12h")), na.rm = TRUE),
    pro_LPS24h_mean = mean(c_across(starts_with("pro_LPS24h")), na.rm = TRUE),
    mRNA_NC_mean = mean(c_across(starts_with("NC-")), na.rm = TRUE),
    mRNA_LPS2h_mean = mean(c_across(starts_with("LPS2h-")), na.rm = TRUE),
    mRNA_LPS12h_mean = mean(c_across(starts_with("LPS12h-")), na.rm = TRUE),
    mRNA_LPS24h_mean = mean(c_across(starts_with("LPS24h-")), na.rm = TRUE)
  ) %>%
  ungroup()

pro_mRNA_ratio_df <- MRL_LPS_diff_2 %>%
  mutate(
    ratio_NC = pro_NC_mean / mRNA_NC_mean,
    ratio_LPS2h = pro_LPS2h_mean / mRNA_LPS2h_mean,
    ratio_LPS12h = pro_LPS12h_mean / mRNA_LPS12h_mean,
    ratio_LPS24h = pro_LPS24h_mean / mRNA_LPS24h_mean
  ) %>%
  dplyr::select(gene_id, ratio_NC, ratio_LPS2h, ratio_LPS12h, ratio_LPS24h)

pro_mRNA_ratio_normalized <- pro_mRNA_ratio_df %>%
  rowwise() %>%
  mutate(
    across(starts_with("ratio_"), ~ (. - min(c_across(starts_with("ratio_")), na.rm=TRUE)) /
             (max(c_across(starts_with("ratio_")), na.rm=TRUE) - min(c_across(starts_with("ratio_")), na.rm=TRUE)))
  ) %>%
  ungroup()

# --- Step 5.2: Perform Mfuzz Clustering ---
heatmap_matrix <- pro_mRNA_ratio_normalized %>%
  dplyr::select(-gene_id) %>%
  as.matrix()
rownames(heatmap_matrix) <- pro_mRNA_ratio_normalized$gene_id

heatmap_matrix_clean <- heatmap_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  group_by(gene_id) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()

valid_rows <- complete.cases(heatmap_matrix_clean) & !apply(heatmap_matrix_clean, 1, function(x) any(is.nan(x)))
heatmap_matrix_clean <- heatmap_matrix_clean[valid_rows, ]

eset <- new('ExpressionSet', exprs = heatmap_matrix_clean)
m_value <- mestimate(eset)
cl <- mfuzz(eset, c = 2, m = m_value)

# --- Step 5.3: Assign Clusters and Prepare Data for Plotting ---
cluster_assignments <- cl$cluster
pro_mRNA_ratio_normalized_clustered <- as.data.frame(heatmap_matrix_clean) %>%
  mutate(gene_id = rownames(.), Cluster = factor(cluster_assignments))

MRL_LPS_diff_2_clustered <- MRL_LPS_diff_2 %>%
  inner_join(pro_mRNA_ratio_normalized_clustered[, c("gene_id", "Cluster")], by = "gene_id")

MRL_LPS_diff_2_long <- MRL_LPS_diff_2_clustered %>%
  dplyr::select(Cluster, NC_2h, LPS_2h, LPS_12h, LPS_24h) %>%
  pivot_longer(
    cols = c(NC_2h, LPS_2h, LPS_12h, LPS_24h),
    names_to = "Condition",
    values_to = "MRL_diff"
  ) %>%
  mutate(Condition = factor(Condition, levels = c("NC_2h", "LPS_2h", "LPS_12h", "LPS_24h")))

# --- Step 5.4: Generate the Plot for FigS5C ---
FigS6D_grouped_plot <- ggplot(MRL_LPS_diff_2_long, aes(x = Condition, y = MRL_diff, fill = Cluster)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "MRL Difference Across Conditions by Cluster",
    x = "Condition",
    y = "MRL Difference"
  ) +
  scale_fill_manual(
    values = c("1" = "#8887b9", "2" = "#95c08b"),
    name = "Cluster"
  ) +
  # Compare clusters within each condition
  stat_compare_means(aes(group = Cluster), label = "p.format", method = "t.test")

# --- Display the Plot ---
print(FigS6D_grouped_plot)

write.csv(MRL_LPS_diff_2_long, file = "FigS6D.csv", quote = F, row.names = F)




# -----------------------------------------------------------------------------
## Section 1: Setup and Environment
# -----------------------------------------------------------------------------
# This section loads necessary libraries and sources your custom R scripts.

# Load required packages for plotting and data manipulation
library(ggplot2)
library(dplyr)

# Source the script containing the custom 'calculate_polysome_load' function
# Ensure this file path is correct for your environment.
source("code/polysome_seq_functions.R")

# -----------------------------------------------------------------------------
## Section 2: Data Loading
# -----------------------------------------------------------------------------
# This section loads the preprocessed sequencing and proteomics data.

# Load the main preprocessed RData object
# Ensure this file path is correct for your environment.
load("data/THP1_preprocessed.RData")

# Source the script that loads the proteomics data (gene ID to gene name mapping)
# This script should create the 'THP1_proteomic_id_protein_name' data frame.
# Ensure this file path is correct for your environment.
source("code/THP1_Proteomic_2024.1.30.R")


# -----------------------------------------------------------------------------
## Section 3: Initial Data Filtering and Cleaning
# -----------------------------------------------------------------------------
# This section filters genes based on expression levels and handles missing values,
# which is a crucial step to ensure data consistency for downstream analysis.

# Filter for genes with logCPM > 1 to remove lowly expressed genes
merge_df_counts_select <- merge_df_counts[merge_df_counts$gene_id %in% genes_logCPM[genes_logCPM$logCPM > 1, "gene_id"], ]
merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts_new$gene_id %in% genes_logCPM[genes_logCPM$logCPM > 1, "gene_id"], ]

# Calculate 'old' mRNA counts by subtracting 'new' from 'total'
merge_df_counts_old_select <- data.frame(gene_id = merge_df_counts_select$gene_id, merge_df_counts_select[, -1] - merge_df_counts_new_select[, -1])

# Calculate MRLs for the initial time course to identify and remove genes with NA values
# This cleaning step is important because it modifies the data frames used later.
MRLs <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))
MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]

# Identify and remove rows that contain any NA values in the MRL calculations
na_rows_combined <- which(rowSums(is.na(MRL_new)) > 0 | rowSums(is.na(MRL_old)) > 0)

if (length(na_rows_combined) > 0) {
  merge_df_counts_select <- merge_df_counts_select[-na_rows_combined, ]
  merge_df_counts_new_select <- merge_df_counts_new_select[-na_rows_combined, ]
  merge_df_counts_old_select <- merge_df_counts_old_select[-na_rows_combined, ]
}


# -----------------------------------------------------------------------------
## Section 4: MRL Calculation for LPS Time Course
# -----------------------------------------------------------------------------
# This section calculates the Mean Ribosome Load (MRL) difference for the
# LPS stimulation conditions, which is the primary data for the plots.

# Calculate MRLs for both new and old mRNAs across LPS time points
MRL_LPS <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_2h", "LPS_2h", "LPS_12h", "LPS_24h"))
MRL_LPS_new <- MRL_LPS[[1]]
MRL_LPS_old <- MRL_LPS[[2]]

# Calculate the MRL difference (New MRL - Old MRL)
MRL_LPS_diff <- MRL_LPS_new - MRL_LPS_old
MRL_LPS_diff$gene_id <- merge_df_counts_select$gene_id

# Merge with proteomics data to add gene names, then remove any remaining NAs
MRL_LPS_diff <- merge(MRL_LPS_diff, THP1_proteomic_id_protein_name, "gene_id")
MRL_LPS_diff <- na.omit(MRL_LPS_diff)

# Define the set of genes to be labeled in all three plots
label_gene_set <- c("IL1B", "IFIT2", "IFIT3", "CXCL10", "JUNB", "MT2A", "CCL5", "RIG-1", "OAS2", "IFIT5", "NT5C3A", "ABCA1", "IFI44L", "SAMD9", "MARCKS")


# -----------------------------------------------------------------------------
## Section 5: Generate and Display Figure 4H Plots
# -----------------------------------------------------------------------------
# This section creates the three scatter plots comparing MRL differences.

### Plot 1: NC 2h vs. LPS 2h
model1 <- lm(LPS_2h ~ NC_2h, data = MRL_LPS_diff)
MRL_LPS_diff$residuals1 <- resid(model1)
std_dev1 <- sd(MRL_LPS_diff$residuals1)
MRL_LPS_diff$highlight_upper1 <- MRL_LPS_diff$residuals1 > (1 * std_dev1)
MRL_LPS_diff$highlight_lower1 <- MRL_LPS_diff$residuals1 < (-1 * std_dev1)
MRL_LPS_diff$label1 <- ifelse(MRL_LPS_diff$gene_name %in% label_gene_set, as.character(MRL_LPS_diff$gene_name), "")
MRL_LPS_diff$shape_factor1 <- ifelse(MRL_LPS_diff$gene_name %in% label_gene_set, "highlighted", "normal")

cor_test_Fig4H_plot1 <- cor.test(MRL_LPS_diff$NC_2h, MRL_LPS_diff$LPS_2h)
cor_test_Fig4H_plot2 <- cor.test(MRL_LPS_diff$NC_2h, MRL_LPS_diff$LPS_12h)
cor_test_Fig4H_plot3 <- cor.test(MRL_LPS_diff$NC_2h, MRL_LPS_diff$LPS_24h)

Fig4H_plot1 <- ggplot(MRL_LPS_diff, aes(x = NC_2h, y = LPS_2h)) +
  geom_point(aes(color = ifelse(gene_name %in% label_gene_set, "#54278f",
                                ifelse(highlight_upper1, "#e95966",
                                       ifelse(highlight_lower1, "#00bad5", "grey"))),
                 shape = shape_factor1,
                 alpha = 0.5, # Alpha adjusted for better point visibility
                 size = ifelse(shape_factor1 == "highlighted", 3, 1.5))) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(aes(label = label1)) + # For better label placement, consider using the ggrepel package
  scale_shape_manual(values = c("normal" = 16, "highlighted" = 17)) +
  labs(x = "MRL Difference (NC 2h)", y = "MRL Difference (LPS 2h)", title = "NC 2h vs LPS 2h") +
  scale_color_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-6, 4.5)) +
  scale_y_continuous(limits = c(-6, 4.5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )


### Plot 2: NC 2h vs. LPS 12h
model2 <- lm(LPS_12h ~ NC_2h, data = MRL_LPS_diff)
MRL_LPS_diff$residuals2 <- resid(model2)
std_dev2 <- sd(MRL_LPS_diff$residuals2)
MRL_LPS_diff$highlight_upper2 <- MRL_LPS_diff$residuals2 > (1 * std_dev2)
MRL_LPS_diff$highlight_lower2 <- MRL_LPS_diff$residuals2 < (-1 * std_dev2)
MRL_LPS_diff$label2 <- ifelse(MRL_LPS_diff$gene_name %in% label_gene_set, as.character(MRL_LPS_diff$gene_name), "")
MRL_LPS_diff$shape_factor2 <- ifelse(MRL_LPS_diff$gene_name %in% label_gene_set, "highlighted", "normal")

Fig4H_plot2 <- ggplot(MRL_LPS_diff, aes(x = NC_2h, y = LPS_12h)) +
  geom_point(aes(color = ifelse(gene_name %in% label_gene_set, "#54278f",
                                ifelse(highlight_upper2, "#e95966",
                                       ifelse(highlight_lower2, "#00bad5", "grey"))),
                 shape = shape_factor2,
                 alpha = 0.5,
                 size = ifelse(shape_factor2 == "highlighted", 3, 1.5))) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(aes(label = label2)) +
  scale_shape_manual(values = c("normal" = 16, "highlighted" = 17)) +
  labs(x = "MRL Difference (NC 2h)", y = "MRL Difference (LPS 12h)", title = "NC 2h vs LPS 12h") +
  scale_color_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-6, 4.5)) +
  scale_y_continuous(limits = c(-6, 4.5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )


### Plot 3: NC 2h vs. LPS 24h
model3 <- lm(LPS_24h ~ NC_2h, data = MRL_LPS_diff)
MRL_LPS_diff$residuals3 <- resid(model3)
std_dev3 <- sd(MRL_LPS_diff$residuals3)
MRL_LPS_diff$highlight_upper3 <- MRL_LPS_diff$residuals3 > (1 * std_dev3)
MRL_LPS_diff$highlight_lower3 <- MRL_LPS_diff$residuals3 < (-1 * std_dev3)
MRL_LPS_diff$label3 <- ifelse(MRL_LPS_diff$gene_name %in% label_gene_set, as.character(MRL_LPS_diff$gene_name), "")
MRL_LPS_diff$shape_factor3 <- ifelse(MRL_LPS_diff$gene_name %in% label_gene_set, "highlighted", "normal")

Fig4H_plot3 <- ggplot(MRL_LPS_diff, aes(x = NC_2h, y = LPS_24h)) +
  geom_point(aes(color = ifelse(gene_name %in% label_gene_set, "#54278f",
                                ifelse(highlight_upper3, "#e95966",
                                       ifelse(highlight_lower3, "#00bad5", "grey"))),
                 shape = shape_factor3,
                 alpha = 0.5,
                 size = ifelse(shape_factor3 == "highlighted", 3, 1.5))) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_text(aes(label = label3)) +
  scale_shape_manual(values = c("normal" = 16, "highlighted" = 17)) +
  labs(x = "MRL Difference (NC 2h)", y = "MRL Difference (LPS 24h)", title = "NC 2h vs LPS 24h") +
  scale_color_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-6, 4.5)) +
  scale_y_continuous(limits = c(-6, 4.5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# Display the final plots
print(Fig4H_plot1)
print(Fig4H_plot2)
print(Fig4H_plot3)

write.csv(MRL_LPS_diff, file = "Fig4H.csv", quote = F, row.names = F)

m6A_FC <- read_excel("data/18_2024_5261_MOESM8_ESM.xlsx", sheet = "M1 to Mϕ ", skip = 1)
colnames(m6A_FC) <- c(colnames(m6A_FC)[1:3],"log2_fold_change_m6A", "gene_name")
load("data/THP1_m6A_IP_seq/MRL_LPS_diff.RData")
MRL_LPS_diff_all_genes <- MRL_LPS_diff
load("data/THP1_m6A_IP_seq/MRL_LPS_diff_2h_for_motif_enrichment.RData")
MRL_LPS_diff_with_protein_select <- MRL_LPS_diff
gene_id_name_pairs <- read_csv("data/THP1_m6A_IP_seq/gene_id_name_pairs.csv")
MRL_LPS_diff_gene_name <- merge(MRL_LPS_diff_all_genes, gene_id_name_pairs, by = "gene_id")
MRL_LPS_diff_m6A <- merge(m6A_FC, MRL_LPS_diff_gene_name, by = "gene_name")

# Load ggplot2
library(ggplot2)

# Create the plot
FigS5H_p1 <- ggplot(MRL_LPS_diff_m6A, aes(x = log2_fold_change_m6A, y = LPS_2h - NC_2h)) +
  geom_point(alpha = 0.3, size = 2, color = "#377eb8") +  # Add scatter points
  labs(
    title = "Scatter Plot: log2 Fold Change vs LPS_2h - NC_2h",
    x = "log2 Fold Change m6A",
    y = "Difference (LPS_2h - NC_2h)"
  ) +
  theme_bw() +  # Use a clean, minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),  # Set axis titles to size 14
    axis.text = element_text(size = 14, color = "black"),  # Set axis labels to size 14 and black
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Set panel border color to black and width to 1
  )

write.csv(MRL_LPS_diff_m6A, file = "FigS5H.csv", quote = F, row.names = F)

# Load necessary library
library(ggplot2)
library(dplyr)

# Filter rows in MRL_LPS_diff_with_protein_select for upper or lower highlights
highlighted_genes <- MRL_LPS_diff_with_protein_select %>%
  filter(highlight_upper1 == TRUE | highlight_lower1 == TRUE) %>%
  dplyr::select(gene_id, highlight_upper1, highlight_lower1)  # Keep relevant columns

# Create a new dataframe from MRL_LPS_diff_m6A based on the highlighted genes
MRL_LPS_diff_m6A_up_down <- MRL_LPS_diff_m6A %>%
  inner_join(highlighted_genes, by = "gene_id")  # Match based on gene_id


# Filter and prepare data for Upper and Lower groups only
filtered_data <- MRL_LPS_diff_m6A_up_down %>%
  mutate(group = case_when(
    highlight_upper1 ~ "Upper",
    highlight_lower1 ~ "Lower"
  )) %>%
  filter(!is.na(group))  # Exclude rows where group is NA (Other)

# Perform a t-test
t_test_result <- t.test(
  log2_fold_change_m6A ~ group, 
  data = filtered_data, 
  var.equal = TRUE # Assume equal variance
)

# Extract the p-value
p_value <- t_test_result$p.value

# Visualization: Boxplot with t-test p-value label
FigS5H_p3 <- ggplot(filtered_data, aes(x = group, y = log2_fold_change_m6A, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("Upper" = "#e95966", "Lower" = "#00bad5")) +
  labs(
    title = "THP-1 LPS activation m6A peaks in Upper vs Lower genes",
    x = "Group in LPS 2 h activation",
    y = "log2 Fold Change m6A peaks"
  ) +
  annotate(
    "text", 
    x = 1.5, 
    y = max(filtered_data$log2_fold_change_m6A, na.rm = TRUE) + 1, 
    label = paste0("p = ", format(p_value, digits = 3)),
    size = 5
  ) +
  theme_bw() +  # Use a clean, minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),  # Adjust axis title font size
    axis.text = element_text(size = 14, color = "black"),  # Adjust axis text size and color
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Add a black border to the plot
  )

write.csv(filtered_data, file = "FigS5H.csv", quote = F, row.names = F)




# get the main isoform for each genes
library(readr)
library(dplyr)

load("data/for_isoform/THP1_isoform/identified_isoform_cleaned.RData")

main_transcript_isoform <- identified_isoform_cleaned %>%
  group_by(gene_id) %>%
  filter(NC_RNA_exp_iso == max(NC_RNA_exp_iso, na.rm = TRUE)) %>%
  dplyr::slice(1) %>%  # In case of ties, take the first one
  ungroup()


load("data/MRL_LPS_diff_2h_for_motif_enrichment.RData")
GRCh38_gene_name_id <- read_csv("data/for_isoform/THP1_isoform/GRCh38_gene_name_id.csv")


# Clean up the `gene_id` column by removing everything after the point
GRCh38_gene_name_id <- GRCh38_gene_name_id %>%
  mutate(gene_id_clean = sub("\\..*", "", gene_id))

MRL_LPS_diff <- MRL_LPS_diff %>%
  mutate(gene_id = sub("\\..*", "", gene_id)) %>%
  filter(highlight_upper1 == TRUE | highlight_lower1 == TRUE)

# Join the two dataframes based on `transcript_id` and `tracking_id` (removing version number from `tracking_id`)
main_transcript_with_gene_id <- main_transcript_isoform %>%
  mutate(transcript_id_clean = sub("\\..*", "", transcript_id))

MRL_LPS_diff_transcript <- merge(MRL_LPS_diff, main_transcript_with_gene_id, by = "gene_id")

library(biomaRt)

# Connect to the Ensembl database for human genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to retrieve sequences for a specific type
get_sequence <- function(transcript_ids, sequence_type) {
  attributes <- c("ensembl_transcript_id", sequence_type)
  
  # Ensure transcript IDs are clean (e.g., no version numbers)
  transcript_ids_clean <- gsub("\\.\\d+$", "", transcript_ids)
  
  # Retrieve sequences
  sequences <- getBM(
    attributes = attributes,
    filters = "ensembl_transcript_id",
    values = transcript_ids_clean,
    mart = ensembl
  )
  
  # Rename the sequence column for clarity
  colnames(sequences)[2] <- sequence_type
  return(sequences)
}

# Example usage:
transcript_ids <- MRL_LPS_diff_transcript$transcript_id  # Replace with your actual transcript IDs

# Retrieve 5'UTR, CDS, and 3'UTR sequences separately
utr5_sequences <- get_sequence(transcript_ids, "5utr")
cds_sequences <- get_sequence(transcript_ids, "coding")
utr3_sequences <- get_sequence(transcript_ids, "3utr")

# Combine all sequences into one dataframe
sequences_combined <- Reduce(function(x, y) merge(x, y, by = "ensembl_transcript_id", all = TRUE),
                             list(utr5_sequences, cds_sequences, utr3_sequences))

colnames(sequences_combined) <- c("transcript_id", "UTR5", "coding", "UTR3")


MRL_LPS_diff_transcript_2 <- merge(MRL_LPS_diff_transcript, sequences_combined, by = "transcript_id")

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# General motif enrichment function
enrich_motif <- function(sequences_df, part_column, motif_pattern, group_name, part_name) {
  # Filter out rows where the sequence part is unavailable or NA
  sequences_df <- sequences_df %>%
    filter(!is.na(!!sym(part_column)) & !!sym(part_column) != "Sequence unavailable")
  
  # Calculate sequence length
  sequences_df <- sequences_df %>% 
    mutate(sequence_length = nchar(!!sym(part_column)))
  
  # Count motifs in each sequence part
  sequences_df <- sequences_df %>%
    rowwise() %>%
    mutate(motif_count = sum(gregexpr(motif_pattern, !!sym(part_column), perl = TRUE)[[1]] > 0))
  
  # Calculate ratio
  sequences_df <- sequences_df %>%
    mutate(ratio = motif_count / sequence_length,
           group = group_name,
           part = part_name)
  
  # Return relevant columns
  return(sequences_df %>% 
           dplyr::select(gene_id, transcript_id, motif_count, sequence_length, ratio, group, part))
}

# Function to process a specific group and all sequence parts
process_group <- function(sequences_df, motif_pattern, group_name) {
  # Enrich for UTR5
  result_utr5 <- enrich_motif(sequences_df, "UTR5", motif_pattern, group_name, "UTR5")
  
  # Enrich for coding
  result_coding <- enrich_motif(sequences_df, "coding", motif_pattern, group_name, "CDS")
  
  # Enrich for UTR3
  result_utr3 <- enrich_motif(sequences_df, "UTR3", motif_pattern, group_name, "UTR3")
  
  # Combine results
  return(bind_rows(result_utr5, result_coding, result_utr3))
}

# Define the DRACH motif pattern
motif_pattern <- "[AGT][AG]AC[ACT]"
RRACH_pattern <- "[AG][AG]AC[ACT]"
YTHDC1 <- "GGAC"
YTHDC2 <- "TGGACT"
YTHDF1 <- "G[AG]AC"
YTHDF2 <- "[AGT]GAC[AT]"
YTHDF3 <- "[AGT][AG]AC[ACT]"


# Filter upper and lower genes
sequences_upper <- MRL_LPS_diff_transcript_2 %>% filter(highlight_upper1)
sequences_lower <- MRL_LPS_diff_transcript_2 %>% filter(highlight_lower1)


#motif_pattern
# Process both groups
result_upper <- process_group(sequences_upper, motif_pattern, "Upper")
result_lower <- process_group(sequences_lower, motif_pattern, "Lower")

# Combine results for all parts and groups
combined_results <- bind_rows(result_upper, result_lower)
combined_results <- combined_results %>%
  mutate(part = factor(part, levels = c("UTR5", "CDS", "UTR3")))
# Create subfigures for each sequence part with custom colors
FigS5G_plot_motif_pattern <- ggplot(combined_results, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot() +
  facet_wrap(~part, scales = "free_y") +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1),  # Change frame color to black
    strip.background = element_blank(),  # Remove grey background from facet headers
    strip.text = element_text(size = 14)  # Adjust facet header text size
  ) +
  ylab("Motif Ratio") +
  xlab("Group") +
  ggtitle("Motif Enrichment Across UTR5, CDS, and UTR3") +
  scale_fill_manual(values = c("Upper" = "#db5763", "Lower" = "#1eafc7")) +
  stat_compare_means(method = "t.test")


write.csv(combined_results, file = "FigS5G.csv", quote = F, row.names = F)



#RRACH_pattern
# Process both groups
result_upper <- process_group(sequences_upper, RRACH_pattern, "Upper")
result_lower <- process_group(sequences_lower, RRACH_pattern, "Lower")

# Combine results for all parts and groups
combined_results <- bind_rows(result_upper, result_lower)
combined_results <- combined_results %>%
  mutate(part = factor(part, levels = c("UTR5", "CDS", "UTR3")))
# Create subfigures for each sequence part with custom colors
FigS5G_plot_RRACH_pattern <- ggplot(combined_results, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot() +
  facet_wrap(~part, scales = "free_y") +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1),  # Change frame color to black
    strip.background = element_blank(),  # Remove grey background from facet headers
    strip.text = element_text(size = 14)  # Adjust facet header text size
  ) +
  ylab("Motif Ratio") +
  xlab("Group") +
  ggtitle("Motif Enrichment Across UTR5, CDS, and UTR3") +
  scale_fill_manual(values = c("Upper" = "#db5763", "Lower" = "#1eafc7")) +
  stat_compare_means(method = "t.test")








## --------------------------------------------------------------
## 1. Load libraries & data
## --------------------------------------------------------------
source("code/polysome_seq_functions.R")
load("data/THP1_preprocessed.RData")
load("data/RNA_features_gene_level.RData")
load("data/THP1_data_abundance.RData")   # contains THP1_data_ensembl
load("pro_mRNA_ratio_df.R")              # contains protein/mRNA ratios

library(dplyr)
library(ggplot2)
library(readr)
library(ggsignif)   # only needed for p3 significance bar

## --------------------------------------------------------------
## 2. Build the core dataframe MRL_LPS_diff
## --------------------------------------------------------------
# ---- Polysome load -------------------------------------------------
merge_df_counts_select <- merge_df_counts[rowMeans(merge_df_counts[,-1]) > 1000, ]
merge_df_counts_new_select  <- merge_df_counts_new[ merge_df_counts$gene_id %in% merge_df_counts_select$gene_id, ]
merge_df_counts_old_select  <- data.frame(gene_id = merge_df_counts_select$gene_id,
                                          merge_df_counts_select[,-1] - merge_df_counts_new_select[,-1])

MRL_LPS <- calculate_polysome_load(
  merge_df_counts_new_select,
  merge_df_counts_old_select,
  c("NC_2h","LPS_2h","LPS_12h","LPS_24h")
)
MRL_LPS_new <- MRL_LPS[[1]]
MRL_LPS_old <- MRL_LPS[[2]]

MRL_LPS_diff <- MRL_LPS_new - MRL_LPS_old
MRL_LPS_diff$gene_id <- merge_df_counts_select$gene_id

# ---- Total RNA-seq (mRNA) -----------------------------------------
THP1_LPS_Total_RNA_seq <- read_delim("data/THP1_LPS_Total_RNA_seq.tsv",
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

countData <- THP1_LPS_Total_RNA_seq[, c(8:19)]
THP1_mRNA <- cbind(THP1_LPS_Total_RNA_seq[, c(1,7)], countData)

count_columns <- c("mRNA_LPS12h-1","mRNA_LPS12h-2","mRNA_LPS12h-3",
                   "mRNA_LPS24h-1","mRNA_LPS24h-2","mRNA_LPS24h-3",
                   "mRNA_LPS2h-1","mRNA_LPS2h-2","mRNA_LPS2h-3",
                   "mRNA_NC-1","mRNA_NC-2","mRNA_NC-3")
colnames(THP1_mRNA)[3:14] <- count_columns

library_sizes <- colSums(THP1_mRNA[, count_columns])
THP1_mRNA_CPM <- THP1_mRNA
for (col in count_columns)
  THP1_mRNA_CPM[[col]] <- (THP1_mRNA[[col]] / library_sizes[col]) * 1e6

THP1_mRNA_CPM <- THP1_mRNA_CPM %>% 
  dplyr::select(Geneid, gene_name, all_of(count_columns))

THP1_mRNA_CPM$gene_id <- sub("\\..*","", THP1_mRNA_CPM$Geneid)

THP1_mRNA_avg <- THP1_mRNA_CPM %>% 
  mutate(
    mRNA_LPS12h = rowMeans(dplyr::select(., `mRNA_LPS12h-1`,`mRNA_LPS12h-2`,`mRNA_LPS12h-3`), na.rm=TRUE),
    mRNA_LPS24h = rowMeans(dplyr::select(., `mRNA_LPS24h-1`,`mRNA_LPS24h-2`,`mRNA_LPS24h-3`), na.rm=TRUE),
    mRNA_LPS2h  = rowMeans(dplyr::select(., `mRNA_LPS2h-1`, `mRNA_LPS2h-2`, `mRNA_LPS2h-3`),  na.rm=TRUE),
    mRNA_NC     = rowMeans(dplyr::select(., `mRNA_NC-1`,    `mRNA_NC-2`,    `mRNA_NC-3`),    na.rm=TRUE)
  ) %>% 
  dplyr::select(gene_id, gene_name, mRNA_LPS12h, mRNA_LPS24h, mRNA_LPS2h, mRNA_NC)

# ---- Protein data -------------------------------------------------
THP1_data_avg <- THP1_data_ensembl %>% 
  mutate(
    pro_NC     = rowMeans(dplyr::select(., `pro_NC-1`,`pro_NC-2`,`pro_NC-3`), na.rm=TRUE),
    pro_LPS2h  = rowMeans(dplyr::select(., `pro_LPS2h-1`,`pro_LPS2h-2`,`pro_LPS2h-3`), na.rm=TRUE),
    pro_LPS12h = rowMeans(dplyr::select(., `pro_LPS12h-1`,`pro_LPS12h-2`,`pro_LPS12h-3`), na.rm=TRUE),
    pro_LPS24h = rowMeans(dplyr::select(., `pro_LPS24h-1`,`pro_LPS24h-2`,`pro_LPS24h-3`), na.rm=TRUE)
  ) %>% 
  dplyr::select(gene_id, accession, pro_NC, pro_LPS2h, pro_LPS12h, pro_LPS24h)

# ---- Combine mRNA + protein ---------------------------------------
protein_mRNA_combined_data <- THP1_mRNA_avg %>% 
  inner_join(THP1_data_avg, by = "gene_id") %>% 
  dplyr::select(gene_id, gene_name, accession,
                mRNA_LPS12h, mRNA_LPS24h, mRNA_LPS2h, mRNA_NC,
                pro_NC, pro_LPS2h, pro_LPS12h, pro_LPS24h)

protein_mRNA_combined_data_with_ratios <- protein_mRNA_combined_data %>% 
  mutate(
    ratio_NC     = ifelse(mRNA_NC == 0, NA, pro_NC / mRNA_NC),
    ratio_LPS2h  = ifelse(mRNA_LPS2h == 0, NA, pro_LPS2h / mRNA_LPS2h),
    ratio_LPS12h = ifelse(mRNA_LPS12h == 0, NA, pro_LPS12h / mRNA_LPS12h),
    ratio_LPS24h = ifelse(mRNA_LPS24h == 0, NA, pro_LPS24h / mRNA_LPS24h)
  )

# ---- Merge with MRL differences -----------------------------------
MRL_LPS_diff$gene_id <- sub("\\..*","", MRL_LPS_diff$gene_id)
MRL_LPS_diff <- merge(MRL_LPS_diff, protein_mRNA_combined_data_with_ratios, by = "gene_id")
MRL_LPS_diff <- na.omit(MRL_LPS_diff)

## --------------------------------------------------------------
## 3. Add columns required for p1 & p3
## --------------------------------------------------------------
# mRNA fold-change (LPS 2h vs NC)
MRL_LPS_diff$mRNA_fc <- MRL_LPS_diff$mRNA_LPS2h / MRL_LPS_diff$mRNA_NC

# Translation-efficiency fold-change
MRL_LPS_diff$TE_fc <- MRL_LPS_diff$ratio_LPS2h / MRL_LPS_diff$ratio_NC

# dynaRDS increase (used in p3)
MRL_LPS_diff$dynaRDS_increase <- MRL_LPS_diff$LPS_2h - MRL_LPS_diff$NC_2h

# Highlight genes >1 SD above/below the NC-2h vs LPS-2h regression line
model1 <- lm(LPS_2h ~ NC_2h, data = MRL_LPS_diff)
MRL_LPS_diff$residuals1       <- resid(model1)
MRL_LPS_diff$highlight_upper1 <- MRL_LPS_diff$residuals1 > sd(MRL_LPS_diff$residuals1)
MRL_LPS_diff$highlight_lower1 <- MRL_LPS_diff$residuals1 < -sd(MRL_LPS_diff$residuals1)

# Global TE-vs-mRNA fit (needed for p3 position relative to line)
fit_global <- lm(log2(TE_fc) ~ log2(mRNA_fc), data = MRL_LPS_diff)

## --------------------------------------------------------------
## 4. Processed dataframe for the two plots
## --------------------------------------------------------------
MRL_LPS_diff_processed <- MRL_LPS_diff %>% 
  mutate(
    highlight_group = case_when(
      highlight_upper1 ~ "Upper",
      highlight_lower1 ~ "Lower",
      TRUE              ~ "Other"
    ),
    position_relative_to_line = if_else(residuals(fit_global) > 0, "Above Line", "Below Line")
  )

# Data filtered for p3 (mRNA_fc > 1, i.e. log2(mRNA_fc) > 0)
plot_data_filtered <- MRL_LPS_diff_processed %>% 
  filter(log2(mRNA_fc) > 0)

## --------------------------------------------------------------
## 5. Plot p1 – Scatter of log2(TE_fc) vs log2(mRNA_fc)
## --------------------------------------------------------------
Fig4J_p1 <- ggplot(MRL_LPS_diff_processed,
             aes(x = log2(mRNA_fc), y = log2(TE_fc))) +
  geom_point(aes(color = highlight_group, size = highlight_group), alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "log2(mRNA fold change)",
       y = "log2(TE fold change)") +
  scale_color_manual(name = "Gene Group",
                     values = c("Upper" = "#e95966",
                                "Lower" = "#00bad5",
                                "Other" = "grey")) +
  scale_size_manual(name = "Gene Group",
                    values = c("Upper" = 2, "Lower" = 2, "Other" = 1.5)) +
  coord_cartesian(xlim = c(min(log2(MRL_LPS_diff_processed$mRNA_fc), na.rm = TRUE), 2),
                  ylim = c(-2, 2)) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1))


## --------------------------------------------------------------
## 6. Plot p3 – Box-plot of dynaRDS_increase by position relative to line
## --------------------------------------------------------------
Fig4J_p3 <- ggplot(plot_data_filtered,
             aes(x = position_relative_to_line, y = dynaRDS_increase,
                 fill = position_relative_to_line)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  scale_fill_manual(name = "Gene Group",
                    values = c("Above Line" = "#e95966", "Below Line" = "#00bad5")) +
  geom_signif(comparisons = list(c("Above Line", "Below Line")),
              test = "t.test", map_signif_level = FALSE,
              y_position = max(plot_data_filtered$dynaRDS_increase, na.rm = TRUE) * 1.1,
              textsize = 5, vjust = -0.5) +
  labs(x = "Gene Group", y = "dynaRDS Increase") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1))


## --------------------------------------------------------------
## 7. Show the plots in the R session
## --------------------------------------------------------------
print(Fig4J_p1)
print(Fig4J_p3)

write.csv(MRL_LPS_diff_processed, file = "Fig4J.csv", quote = F, row.names = F)



#Fig4E

top_bottom_10_c2 <- read_csv("data/Fig4E_top_bottom_10_c2_enrichment_mod.csv")
top_bottom_10_c2 <- top_bottom_10_c2 %>%
  mutate(Cluster = ifelse(NES > 0, "Fast", "Slow")) %>%
  # Add 'abs_NES' column as the absolute value of 'NES'
  mutate(abs_NES = abs(NES))

# View the updated dataframe
print(top_bottom_10_c2)
#设置celltype展示的顺序
top_bottom_10_c2$Cluster <- factor(top_bottom_10_c2$Cluster, levels = c("Fast","Slow"))
# 使用排序索引重新排列数据框
top_bottom_10_c2 <- top_bottom_10_c2[order(top_bottom_10_c2$Cluster), ]
#terms因子顺序
top_bottom_10_c2$Description <- factor(top_bottom_10_c2$Description, levels = top_bottom_10_c2$Description)

#展示的基因，我们选择每个terms展示5个基因，实际情况可以展示自己关注的基因
top_bottom_10_c2$geneID  <- sapply(strsplit(top_bottom_10_c2$geneID , "/"), function(x) paste(x[1:5], collapse = "/"))

library(ggplot2)

#ggplot作图(就是柱状图)
Fig4E_p2 = ggplot(top_bottom_10_c2, aes(x = abs_NES, y = rev(Description), fill = Cluster))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(x=0.1,y=rev(Description),label = Description),size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#8887b9", "#95c08b"))+
  geom_text(data = top_bottom_10_c2,
            aes(x = 0.1, y = rev(Description), label = geneID, color = Cluster),
            size = 4,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3)+
  scale_color_manual(values = c("#8887b9", "#95c08b"))+
  scale_y_discrete(expand = c(0.1,0))+
  ggtitle("curated gene enrichment")


write.csv(top_bottom_10_c2, file = "Fig4E.csv", quote = F, row.names = F)






