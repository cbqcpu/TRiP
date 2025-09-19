
library(ggplot2)
library(dplyr)
library(tidyr)
library(pracma)
library(gridExtra)

base_theme <- theme_bw() + theme(
  legend.position = "none",  # Remove legend
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
  panel.grid.major = element_blank(),  # Remove major grid lines
  panel.grid.minor = element_blank(),  # Remove minor grid lines
  panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
)

source("code/polysome_seq_functions.R")
load("data/HEK293T_preprocessed.RData")

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

auc_diff_ranked <- rank(auc_diff, ties.method = "first")
high_coupling_indices <- which(auc_diff_ranked > length(auc_diff_ranked) - 1000)
low_coupling_indices <- which(auc_diff_ranked <= 1000)

high_coupling_value <- min(auc_diff[high_coupling_indices])
low_coupling_value <- max(auc_diff[low_coupling_indices])
high_coupling_genes <- merge_df_counts_select[high_coupling_indices,1]
low_coupling_genes <- merge_df_counts_select[low_coupling_indices,1]


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
Fig2B <- ggplot() +
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
  base_theme


MRL_new_high_coupling <- MRL_new[high_coupling_indices, ]
MRL_old_high_coupling <- MRL_old[high_coupling_indices, ]

MRL_new_low_coupling <- MRL_new[low_coupling_indices, ]
MRL_old_low_coupling <- MRL_old[low_coupling_indices, ]


# Function to calculate standard error of the mean
sem <- function(x) sd(x, na.rm = TRUE)

# High Coupling
mean_new_high <- colMeans(MRL_new_high_coupling)
sem_new_high <- apply(MRL_new_high_coupling, 2, sem)
mean_old_high <- colMeans(MRL_old_high_coupling)
sem_old_high <- apply(MRL_old_high_coupling, 2, sem)

data_high_coupling <- data.frame(
  Time = rep(time_points, 2),
  Mean = c(mean_new_high, mean_old_high),
  SEM = c(sem_new_high, sem_old_high),
  Group = rep(c("New", "Old"), each = length(time_points))
)
data_high_coupling$Time <- factor(data_high_coupling$Time, levels = time_points)

# Low Coupling
mean_new_low <- colMeans(MRL_new_low_coupling)
sem_new_low <- apply(MRL_new_low_coupling, 2, sem)
mean_old_low <- colMeans(MRL_old_low_coupling)
sem_old_low <- apply(MRL_old_low_coupling, 2, sem)

data_low_coupling <- data.frame(
  Time = rep(time_points, 2),
  Mean = c(mean_new_low, mean_old_low),
  SEM = c(sem_new_low, sem_old_low),
  Group = rep(c("New", "Old"), each = length(time_points))
)
data_low_coupling$Time <- factor(data_low_coupling$Time, levels = time_points)

# Plotting Functions
plot_with_errorbar <- function(data, title) {
  ggplot(data, aes(x = Time, y = Mean, group = Group, color = Group)) +
    geom_point(aes(color = Group), size = 4, shape = 1, stroke = 1.5) +  # Increase point size here
    geom_line(aes(color = Group), linewidth = 1.5) + # Increase line size here
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
    base_theme +
    scale_color_manual(values = c("New" = "#316892", "Old" = "#92cff1")) + # Adjust the colors here
    labs(title = title, x = "Labeling Time", y = "MRL", color = "mRNAs")
}

# Create Plots
Fig2C_plot_high_coupling <- plot_with_errorbar(data_high_coupling, "High Coupling mRNAs")
Fig2C_plot_low_coupling <- plot_with_errorbar(data_low_coupling, "Low Coupling mRNAs")


MRL_diff <- MRL_new-MRL_old

MRL_diff_output <- MRL_diff
MRL_diff_output$gene_id <- merge_df_counts_select$gene_id
#write.table(MRL_diff_output, file = "gene_level_MRL_diff.csv", quote = FALSE, sep = ",")

MRL_diff$auc_diff <- auc_diff
MRL_diff <- MRL_diff[(!is.infinite(MRL_diff$NC_30min))&(!is.infinite(MRL_diff$NC_1h))&(!is.infinite(MRL_diff$NC_2h)), ]
# Calculate correlation and p-value for each time point
cor_test_30min <- cor.test(MRL_diff$auc_diff, MRL_diff$NC_30min)
cor_30min <- cor_test_30min$estimate
pval_30min <- cor_test_30min$p.value

cor_test_1h <- cor.test(MRL_diff$auc_diff, MRL_diff$NC_1h)
cor_1h <- cor_test_1h$estimate
pval_1h <- cor_test_1h$p.value

cor_test_2h <- cor.test(MRL_diff$auc_diff, MRL_diff$NC_2h)
cor_2h <- cor_test_2h$estimate
pval_2h <- cor_test_2h$p.value

# Plot for 30 min
FigS2C_plot1 <- ggplot(MRL_diff, aes(x = auc_diff, y = NC_30min)) +
  stat_density_2d(geom = "path", color = "#99c3e5") +
  labs(x = "AUC Difference", y = "MRL Difference", title = "NC 30min") +
  base_theme +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("r = ", round(cor_30min, 2), 
                          "\nP = ", format.pval(pval_30min, digits = 2, eps = 0.001)), 
           hjust = 1.1, vjust = 1.5, size = 3)

# Plot for 1h
FigS2C_plot2 <- ggplot(MRL_diff, aes(x = auc_diff, y = NC_1h)) +
  stat_density_2d(geom = "path", color = "#99c3e5") +
  labs(x = "AUC Difference", y = "MRL Difference", title = "NC 1h") +
  base_theme +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("r = ", round(cor_1h, 2), 
                          "\nP = ", format.pval(pval_1h, digits = 2, eps = 0.001)), 
           hjust = 1.1, vjust = 1.5, size = 3)

# Plot for 2h
FigS2C_plot3 <- ggplot(MRL_diff, aes(x = auc_diff, y = NC_2h)) +
  stat_density_2d(geom = "path", color = "#99c3e5") +
  labs(x = "AUC Difference", y = "MRL Difference", title = "NC 2h") +
  base_theme +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("r = ", round(cor_2h, 2), 
                          "\nP = ", format.pval(pval_2h, digits = 2, eps = 0.001)), 
           hjust = 1.1, vjust = 1.5, size = 3)



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


library(readr)
set.seed(1234)  # You can use any integer as the seed

MRL_df <- cbind(MRL_new, MRL_old)
row_max <- apply(MRL_df, 1, max)
MRL_new_normalized <- MRL_new / row_max
MRL_old_normalized <- MRL_old / row_max
colnames(MRL_new_normalized) <- c("new_30_min", "new_1h", "new_2h")
colnames(MRL_old_normalized) <- c("old_30_min", "old_1h", "old_2h")
data_matrix <- as.matrix(cbind(MRL_new_normalized, MRL_old_normalized))
rownames(data_matrix) <- merge_df_counts_select$gene_id


library(Mfuzz)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)

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


# Create a new dataframe including sorted data and cluster labels
data_for_analysis <- data.frame(sorted_data)
data_for_analysis$cluster <- sorted_clusters
data_for_analysis$gene_id <- rownames(data_for_analysis)
rownames(data_for_analysis) <- NULL

library(pracma)

# Provided time points
time_points <- c(0.5, 1, 2)  # Corresponding to 30_min, 1h, and 2h

# Function to calculate AUC for each row
calculate_auc <- function(new_vals, old_vals, time_points) {
  auc_new <- trapz(time_points, new_vals)  # AUC for 'new' columns
  auc_old <- trapz(time_points, old_vals)  # AUC for 'old' columns
  return(c(auc_new, auc_old))
}

# Apply the AUC calculation for each row
data_for_analysis <- data_for_analysis %>%
  rowwise() %>%
  mutate(
    AUC_new = calculate_auc(c_across(starts_with("new_")), c_across(starts_with("old_")), time_points)[1],
    AUC_old = calculate_auc(c_across(starts_with("new_")), c_across(starts_with("old_")), time_points)[2]
  ) %>%
  ungroup()

# View the updated dataframe
head(data_for_analysis)

# Assuming `data_for_analysis` already has `auc_diff` calculated
data_for_analysis$auc_diff <- data_for_analysis$AUC_new - data_for_analysis$AUC_old

library(dplyr)


# Step 5: Combine top and bottom rankings
#ranked_metric <- c(top_1000_ranks, bottom_1000_ranks)
ranked_metric <- data_for_analysis$auc_diff
# Combine the gene IDs
#ranked_genes <- bind_rows(top_1000_genes_cluster1, bottom_1000_genes_cluster2)
ranked_genes <- data_for_analysis
# Assign the custom rankings to the gene IDs
names(ranked_metric) <- ranked_genes$gene_id
ranked_metric <- sort(ranked_metric, decreasing = TRUE)

# Step 6: Map Ensembl IDs to Gene Symbols
library(clusterProfiler)
library(org.Hs.eg.db)

# Extract Ensembl IDs and remove version numbers (if present)
ensembl_ids <- names(ranked_metric)
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = ensembl_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

# Remove entries without gene symbols
valid_indices <- !is.na(gene_symbols)
ranked_metric <- ranked_metric[valid_indices]
names(ranked_metric) <- gene_symbols[valid_indices]

# Step 7: Perform GSEA analysis using the ranked_metric
library(GSEABase)   # For handling GMT files
library(enrichplot) # For visualization

# Load the GMT files for hallmark and transcription factor targets
gmt_file_hallmark <- "data/GSEA_sets/h.all.v2024.1.Hs.symbols.gmt"
gmt_c2 <- "data/GSEA_sets/c2.all.v2024.1.Hs.symbols.gmt"

gene_sets_hallmark <- read.gmt(gmt_file_hallmark)
gene_sets_c2 <- read.gmt(gmt_c2)

# Perform GSEA on Hallmark gene sets
gsea_result_hallmark <- GSEA(ranked_metric, 
                             TERM2GENE = gene_sets_hallmark, 
                             pvalueCutoff = 0.05, 
                             verbose = FALSE)


# Extract enrichment results as data frames
result_hallmark <- as.data.frame(gsea_result_hallmark@result)

result_hallmark <- result_hallmark %>%
  filter(qvalues < 1)


# Assuming `combined_results` contains the GSEA results from multiple gene sets
# with the columns 'ID', 'NES', 'qvalue', and 'Gene_Set'

# Step 1: Filter for positive NES (NES > 0) and sort by NES descending (top 10)
top_10_positive_hallmark <- result_hallmark %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  dplyr::slice(1:10)

# Step 2: Filter for negative NES (NES < 0) and sort by NES ascending (bottom 10)
bottom_10_negative_hallmark <- result_hallmark %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  dplyr::slice(1:10)

# Step 3: Combine the top 10 positive and bottom 10 negative results
top_bottom_10_hallmark <- bind_rows(top_10_positive_hallmark, bottom_10_negative_hallmark)
library(tools)
top_bottom_10_hallmark$Description <- top_bottom_10_hallmark$Description %>%
  # Step 1: Remove 'HALLMARK_' prefix
  gsub("HALLMARK_", "", .) %>%
  # Step 2: Replace underscores with spaces
  gsub("_", " ", .) %>%
  # Step 3: Capitalize each word
  toTitleCase()
#write.csv(top_bottom_10_hallmark, file = "data/top_bottom_10_hallmark_enrichment.csv", quote = FALSE, row.names = FALSE)
top_bottom_10_hallmark <- read_csv("data/figure1_data_sheet/figure1j_top_bottom_10_hallmark_enrichment_mod.csv")
top_bottom_10_hallmark <- top_bottom_10_hallmark %>%
  mutate(Cluster = ifelse(NES > 0, "Fast", "Slow")) %>%
  # Add 'abs_NES' column as the absolute value of 'NES'
  mutate(abs_NES = abs(NES))

# View the updated dataframe
print(top_bottom_10_hallmark)
#设置celltype展示的顺序
top_bottom_10_hallmark$Cluster <- factor(top_bottom_10_hallmark$Cluster, levels = c("Fast","Slow"))
# 使用排序索引重新排列数据框
top_bottom_10_hallmark <- top_bottom_10_hallmark[order(top_bottom_10_hallmark$Cluster), ]
#terms因子顺序
top_bottom_10_hallmark$Description <- factor(top_bottom_10_hallmark$Description, levels = top_bottom_10_hallmark$Description)

#展示的基因，我们选择每个terms展示5个基因，实际情况可以展示自己关注的基因
top_bottom_10_hallmark$geneID  <- sapply(strsplit(top_bottom_10_hallmark$geneID , "/"), function(x) paste(x[1:5], collapse = "/"))

library(ggplot2)

#ggplot作图(就是柱状图)
Fig2D_p1 = ggplot(top_bottom_10_hallmark, aes(x = abs_NES, y = rev(Description), fill = Cluster))+
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
  geom_text(data = top_bottom_10_hallmark,
            aes(x = 0.1, y = rev(Description), label = geneID, color = Cluster),
            size = 4,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3)+
  scale_color_manual(values = c("#8887b9", "#95c08b"))+
  scale_y_discrete(expand = c(0.1,0))+
  ggtitle("Hallmark gene enrichment")



# c2 gene set enrichment
# Perform GSEA on Hallmark gene sets
gsea_result_c2 <- GSEA(ranked_metric, 
                       TERM2GENE = gene_sets_c2, 
                       pvalueCutoff = 0.05, 
                       verbose = FALSE)


# Extract enrichment results as data frames
result_c2 <- as.data.frame(gsea_result_c2@result)

# Add a column to identify the source of each result
result_c2$Gene_Set <- "ontology gene sets"


result_c2 <- result_c2 %>%
  filter(qvalues < 1)


# Assuming `combined_results` contains the GSEA results from multiple gene sets
# with the columns 'ID', 'NES', 'qvalue', and 'Gene_Set'

# Step 1: Filter for positive NES (NES > 0) and sort by NES descending (top 10)
top_10_positive_c2 <- result_c2 %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  dplyr::slice(1:10)

# Step 2: Filter for negative NES (NES < 0) and sort by NES ascending (bottom 10)
bottom_10_negative_c2 <- result_c2 %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  dplyr::slice(1:10)

# Step 3: Combine the top 10 positive and bottom 10 negative results
top_bottom_10_c2 <- bind_rows(top_10_positive_c2, bottom_10_negative_c2)

top_bottom_10_c2$Description <- top_bottom_10_c2$Description %>%
  gsub("_", " ", .) %>%
  # Step 3: Capitalize each word
  toTitleCase()
#write.csv(top_bottom_10_c2, file = "data/top_bottom_10_c2_enrichment.csv", quote = FALSE, row.names = FALSE)
top_bottom_10_c2 <- read_csv("data/figure1_data_sheet/figureS1i_top_bottom_10_c2_enrichment_mod.csv")
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
FigS2D_p2 = ggplot(top_bottom_10_c2, aes(x = abs_NES, y = rev(Description), fill = Cluster))+
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




library(readr)
library(readxl)
load("data/protein_combined_df_gene_id_20240418_abndance.RData")
#load("data/protein_combined_df_gene_id.RData")
load("data/protein_length_20240418.RData")
load("data/HEK293T_counts_MRL.RData")
HEK293T_halflife <- read_excel("data/1-s2.0-S1097276521007498-mmc2.xlsx", sheet = "HEK293T")
colnames(protein_length)[colnames(protein_length)=='uniprotswissprot'] <- "accession"
half_life <- data.frame(accession = HEK293T_halflife$`Protein Id`, `half_life` = HEK293T_halflife$CHX_8...16/HEK293T_halflife$CHX_0...4)
half_life[half_life$half_life>1,2] <- 1

protein_combined_df_gene_id <- merge(protein_combined_df_gene_id, protein_length)
protein_combined_df_gene_id <- merge(protein_combined_df_gene_id, half_life)
mRNA_fractions <- merge_df_counts_select[,c("NC_2h_frac1","NC_2h_frac2","NC_2h_frac3","NC_2h_frac4")]
mRNA_expression <- rowSums(mRNA_fractions)/4
mRNA_expression <- data.frame(gene_id = merge_df_counts_select$gene_id, NC_2h_mRNA = mRNA_expression)
mRNA_expression$gene_id <- sub("\\..*", "", mRNA_expression$gene_id)
MRL_total <- (MRL_new+MRL_old)/2 

mRNA_expression$MRL <- MRL_total$NC_2h
MRL_diff <- MRL_new - MRL_old

library(pracma)
time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new <- apply(MRL_new, 1, calculate_auc)
auc_MRL_old <- apply(MRL_old, 1, calculate_auc)

auc_diff <- auc_MRL_new - auc_MRL_old

mRNA_expression$AUC_diff <- auc_diff
mRNA_expression$AUC_total <- auc_MRL_new + auc_MRL_old

mRNA_expression_new <- read_delim("data/HEK293T_NC.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
mRNA_expression_new <- mRNA_expression_new[,c("Geneid","NC-1","NC-2","NC-3")]
colnames(mRNA_expression_new) <- c("gene_id","NC-1","NC-2","NC-3")
mRNA_expression_new$counts_new <- rowSums(mRNA_expression_new[,2:4])/3
mRNA_expression_new$gene_id <- sub("\\..*", "", mRNA_expression_new$gene_id)

mRNA_expression_2 <- merge(mRNA_expression, mRNA_expression_new, by = "gene_id")

protein_mRNA_compare <- merge(mRNA_expression_2, protein_combined_df_gene_id, by = "gene_id")
protein_mRNA_compare <- protein_mRNA_compare %>% distinct(gene_id, .keep_all = TRUE)

library(ggplot2)

# Clean data to remove NA, NaN, and Inf for correlation and plotting
clean_data <- protein_mRNA_compare[
  complete.cases(protein_mRNA_compare) &
    is.finite(protein_mRNA_compare$AUC_diff) &
    is.finite(protein_mRNA_compare$MRL) &
    is.finite(log2(protein_mRNA_compare$abundance / protein_mRNA_compare$counts_new)), ]

# First Plot: AUC Difference vs Log(Abundance / Counts New)
auc_corr <- cor.test(
  clean_data$AUC_diff, 
  log2(clean_data$abundance / clean_data$counts_new)
)

AUC_diff_plot <- ggplot(clean_data, aes(x = AUC_diff, y = log2(abundance / counts_new))) +
  geom_point(color = "#4ba2dd", alpha = 0.1) +  # Points with user-friendly color
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +  # Add linear fitting line
  annotate("text", x = max(clean_data$AUC_diff, na.rm = TRUE) * 0.7, 
           y = max(log2(clean_data$abundance / clean_data$counts_new), na.rm = TRUE) * 0.9,
           label = paste0("r = ", round(auc_corr$estimate, 3)), 
           size = 5) +  # Annotate correlation coefficient
  labs(
    x = "AUC Difference",
    y = "Log2(Protein / mRNA)",
    title = "Scatter Plot of AUC Difference vs Log2(Protein / mRNA)"
  ) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )


library(ggplot2)
library(FNN)  # For k-nearest neighbors

# Create plotting dataframe
plot_data <- data.frame(
  AUC_diff = clean_data$AUC_diff,
  ratio_log2 = log2(clean_data$abundance / clean_data$counts_new)
)

# Calculate correlation and p-value
cor_results <- cor.test(plot_data$AUC_diff, plot_data$ratio_log2, method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10  # Number of neighbors to consider
neighbors <- get.knn(data = plot_data, k = k)
plot_data$density <- 1 / rowMeans(neighbors$nn.dist)  # Density = inverse of average distance to neighbors
plot_data$density <- plot_data$density / max(plot_data$density)  # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Plot with density-colored points and regression line
Fig2F_AUC_diff_plot <- ggplot(plot_data, aes(x = AUC_diff, y = ratio_log2, color = density)) +
  geom_point(size = 2, alpha = 0.7) +  # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +  # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.4, 1), name = "Density") +  # Custom color bar
  labs(
    x = "AUC Difference",
    y = "Log2(Protein / mRNA)",
    title = "Scatter Plot of AUC Difference vs Log2(Protein / mRNA)"
  ) +
  annotate("text", 
           x = max(plot_data$AUC_diff, na.rm = TRUE) * 0.7, 
           y = max(plot_data$ratio_log2, na.rm = TRUE) * 0.9,
           label = paste0("r = ", cor_value, "\nP = ", p_value), 
           size = 5, color = "black", fontface = "bold") +  # Correlation annotation
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )



# get the gene type by metascape
high_coupling_genes_metascape <- read_excel("data/metascape_result.ts30_g8y2_high_coupling.xlsx")
low_coupling_genes_metascape <- read_excel("data/metascape_result.t802z8gu5_low_coupling.xlsx")

source("code/RNA_function.R")

result_high <- summary_gene_function(high_coupling_genes_metascape)
result_low <- summary_gene_function(low_coupling_genes_metascape)

function_df <- merge(result_high, result_low, by = "Element", all = TRUE, suffixes = c(".high", ".low"))
function_df[is.na(function_df)] <- 0  # Replace NA with 0
colnames(function_df) <- c("Element", "High", "Low")
function_df$High_ratio <- function_df$High/sum(function_df$High)
function_df$Low_ratio <- function_df$Low/sum(function_df$Low)
function_df_enzymes <- function_df[function_df$Element %in% c("ENZYME proteins", "Enzymes", "Kinases"),-1]
function_df_enzymes <- data.frame(t(colSums(function_df[function_df$Element %in% c("ENZYME proteins", "Enzymes", "Kinases"),-1])))
function_df_enzymes$Element <- "Enzymes"

function_df_enzymes <- function_df_enzymes %>%
  dplyr::select(Element, everything())

function_df_select <- function_df[function_df$Element %in% c("Transcription factors", "Predicted secreted proteins", "CD markers", "Transporters"),]
function_df_select <- rbind(function_df_select,function_df_enzymes)


library(reshape2)
function_df_select <- function_df_select[, c("Element", "High_ratio", "Low_ratio")]
function_df_long <- melt(function_df_select, id.vars = "Element", variable.name = "Type", value.name = "Count")

# Plot
Fig2E <- ggplot(function_df_long, aes(x = Element, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )+
  labs(y = "Ratio", x = "Element", title = "Ratio of functional genes in high and low coupling") +
  scale_fill_manual(values = c("Low_ratio" = "#95c08b", "High_ratio" = "#8887b9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

result_high_location <- summary_subcellular_location(high_coupling_genes_metascape)
result_low_location <- summary_subcellular_location(low_coupling_genes_metascape)

result_high_location_cleaned <- result_high_location %>%
  # Remove everything from " (" onwards in the Element column
  mutate(Element = sub(" \\(.*$", "", Element)) %>%
  # Group by the cleaned Element column
  group_by(Element) %>%
  # Summarize the Count for each group
  summarize(Count = sum(Count)) %>%
  # Ungroup the data after summarizing
  ungroup()

result_low_location_cleaned <- result_low_location %>%
  # Remove everything from " (" onwards in the Element column
  mutate(Element = sub(" \\(.*$", "", Element)) %>%
  # Group by the cleaned Element column
  group_by(Element) %>%
  # Summarize the Count for each group
  summarize(Count = sum(Count)) %>%
  # Ungroup the data after summarizing
  ungroup()

location_df <- merge(result_high_location_cleaned, result_low_location_cleaned, by = "Element", all = TRUE, suffixes = c(".high", ".low"))
location_df[is.na(location_df)] <- 0  # Replace NA with 0
colnames(location_df) <- c("Element", "High", "Low")

location_df$High_ratio <- location_df$High/sum(location_df$High)
location_df$Low_ratio <- location_df$Low/sum(location_df$Low)
location_df <- location_df[, c("Element", "High_ratio", "Low_ratio")]
location_df_long <- melt(location_df, id.vars = "Element", variable.name = "Type", value.name = "Count")

# Plot
FigS2E <- ggplot(location_df_long, aes(x = Element, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )+
  labs(y = "Ratio", x = "Element", title = "Ratio of protein Subcellular Location in high and low coupling") +
  scale_fill_manual(values = c("Low_ratio" = "#95c08b", "High_ratio" = "#8887b9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# FigS2F-G
library(readr)
library(dplyr)
library(stringr)

# Example usage:
renamelist_NC30min <- 'data/for_isoform/renamelist_NC30min.csv'
renamelist_NC1h <- 'data/for_isoform/renamelist_NC1h.csv'
renamelist_NC2h <- 'data/for_isoform/renamelist_NC2h_lcd.csv'

NC1h_yeast_all_counts <- "data/for_isoform/NC1h_yeast_all_counts.csv"
NC2h_yeast_all_counts <- "data/for_isoform/NC2h_yeast_all_counts.csv" # This sample use: disk1/cbq/NGS_data/20230606_2h_newdata/NC2h_lcd_*
NC30min_yeast_all_counts <- "data/for_isoform/NC30min_yeast_all_counts.csv"

source("code/functions_HEK293_new_old_isoform_diff_2024.3.7.R")

all_samples_df_30min <- process_normalize_counts(renamelist_NC30min, NC30min_yeast_all_counts, 'data/for_isoform/cuffnorm_counts')
all_samples_df_NC1h <- process_normalize_counts(renamelist_NC1h, NC1h_yeast_all_counts, 'data/for_isoform/cuffnorm_counts')
all_samples_df_NC2h <- process_normalize_counts(renamelist_NC2h, NC2h_yeast_all_counts, 'data/for_isoform/cuffnorm_counts')

load("data/for_isoform/HEK293_transcript_identify_clear.RData")
identified_transcripts <- data.frame(transcript_identify_clear[,3:6])
colnames(identified_transcripts) <- c("gene_id_filter", "tracking_id", "gene_name_filter", "transcript_name")
all_samples_df_30min_filter <- all_samples_df_30min[all_samples_df_30min$tracking_id %in% identified_transcripts$tracking_id,]
all_samples_df_1h_filter <- all_samples_df_NC1h[all_samples_df_NC1h$tracking_id %in% identified_transcripts$tracking_id,]
all_samples_df_2h_filter <- all_samples_df_NC2h[all_samples_df_NC2h$tracking_id %in% identified_transcripts$tracking_id,]

all_samples_df_30min <- all_samples_df_30min_filter
all_samples_df_NC1h <- all_samples_df_1h_filter
all_samples_df_NC2h <- all_samples_df_2h_filter

NC30min_merge_reps <- merge_reps(all_samples_df_30min_filter, "NC30min", c(1,5,6))
NC1h_merge_reps <- merge_reps(all_samples_df_1h_filter, "NC1h", c(1,2,3))
NC2h_merge_reps <- merge_reps(all_samples_df_2h_filter, "NC2h", c(1,2,3))

source("code/polysome_seq_functions.R")
NC30min_new <- NC30min_merge_reps[,grep("new", colnames(NC30min_merge_reps), value = TRUE)]
NC30min_old <- NC30min_merge_reps[,grep("old", colnames(NC30min_merge_reps), value = TRUE)]
colnames(NC30min_new) <- c("NC30min_1","NC30min_2","NC30min_3","NC30min_4")
colnames(NC30min_old) <- c("NC30min_1","NC30min_2","NC30min_3","NC30min_4")
MRL_NC30min <- calculate_polysome_load(NC30min_new, NC30min_old, "NC30min")

NC1h_new <- NC1h_merge_reps[,grep("new", colnames(NC1h_merge_reps), value = TRUE)]
NC1h_old <- NC1h_merge_reps[,grep("old", colnames(NC1h_merge_reps), value = TRUE)]
colnames(NC1h_new) <- c("NC1h_1","NC1h_2","NC1h_3","NC1h_4")
colnames(NC1h_old) <- c("NC1h_1","NC1h_2","NC1h_3","NC1h_4")
MRL_NC1h <- calculate_polysome_load(NC1h_new, NC1h_old, "NC1h")

NC2h_new <- NC2h_merge_reps[,grep("new", colnames(NC2h_merge_reps), value = TRUE)]
NC2h_old <- NC2h_merge_reps[,grep("old", colnames(NC2h_merge_reps), value = TRUE)]
colnames(NC2h_new) <- c("NC2h_1","NC2h_2","NC2h_3","NC2h_4")
colnames(NC2h_old) <- c("NC2h_1","NC2h_2","NC2h_3","NC2h_4")
MRL_NC2h <- calculate_polysome_load(NC2h_new, NC2h_old, "NC2h")


MRL_diff_NC30min <- MRL_NC30min$new_polysome_load-MRL_NC30min$old_polysome_load
MRL_diff_NC1h <- MRL_NC1h$new_polysome_load-MRL_NC1h$old_polysome_load
MRL_diff_NC2h <- MRL_NC2h$new_polysome_load-MRL_NC2h$old_polysome_load


isoform_expression <- data.frame(expression = rowSums(cbind(all_samples_df_30min[,-c(1:4)], all_samples_df_NC1h[,-c(1:4)], all_samples_df_NC2h[,-c(1:4)]))/3)
isoform_expression$transcript_id <- all_samples_df_30min$tracking_id
isoform_expression_select <- isoform_expression[isoform_expression$expression>10,]

isoform_counts_new_select <- cbind(NC30min_new, NC1h_new, NC2h_new)
isoform_counts_old_select <- cbind(NC30min_old, NC1h_old, NC2h_old)

isoform_level_MRLs <- calculate_polysome_load(isoform_counts_new_select, isoform_counts_old_select, c("NC30min", "NC1h", "NC2h"))


library(pracma)
time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new_iso <- apply(isoform_level_MRLs$new_polysome_load, 1, calculate_auc)
auc_MRL_old_iso <- apply(isoform_level_MRLs$old_polysome_load, 1, calculate_auc)

auc_diff_iso <- auc_MRL_new_iso - auc_MRL_old_iso
auc_diff_iso_df <- data.frame(all_samples_df_30min_filter[,1:2], auc_diff_iso)
load("data/for_isoform/HEK293T_counts_MRL.RData")

auc_MRL_new <- apply(MRL_new, 1, calculate_auc)
auc_MRL_old <- apply(MRL_old, 1, calculate_auc)

auc_diff <- auc_MRL_new - auc_MRL_old

auc_diff_df <- data.frame(gene_id = merge_df_counts_select$gene_id, auc_diff)

auc_merge <- merge(auc_diff_iso_df, auc_diff_df, by = "gene_id")

auc_merge_clean <- auc_merge %>% filter(complete.cases(.))

MRL_diff_NC <- data.frame(MRL_diff_iso_NC30min = MRL_diff_NC30min$NC30min, MRL_diff_iso_NC1h = MRL_diff_NC1h$NC1h, MRL_diff_iso_NC2h = MRL_diff_NC2h$NC2h)
MRL_diff_NC$transcript_id <- all_samples_df_NC2h$tracking_id
MRL_diff_NC$gene_id <- all_samples_df_NC2h$gene_id
MRL_diff_NC$gene_name <- all_samples_df_NC2h$gene_name

gene_level_MRL_diff <- read_csv("data/for_isoform/gene_level_MRL_diff.csv")
#gene_level_MRL_diff <- gene_level_MRL_diff[,-1]
MRL_diff_NC_combine <- merge(MRL_diff_NC, gene_level_MRL_diff, by = "gene_id")
MRL_diff_NC_combine <- MRL_diff_NC_combine[MRL_diff_NC_combine$transcript_id %in% isoform_expression_select$transcript_id,]
MRL_diff_NC_combine <- merge(MRL_diff_NC_combine, isoform_expression_select, by = "transcript_id")
MRL_diff_NC_combine <- na.omit(MRL_diff_NC_combine)

auc_merge_clean$transcript_id <- auc_merge_clean$tracking_id
MRL_diff_NC_combine_with_AUC_diff <- merge(MRL_diff_NC_combine, auc_merge_clean, by = "transcript_id")

cor_test_S2F <- cor.test(MRL_diff_NC_combine$NC_2h, MRL_diff_NC_combine$MRL_diff_iso_NC2h)
FigS2F <- ggplot(data = MRL_diff_NC_combine, aes(x = NC_2h, y = MRL_diff_iso_NC2h)) +
  geom_point(color = "#56B4E9", alpha = 0.5) +  # Plot points with a specified color
  geom_density_2d(color = "#D55E00") +  # Add density contours with a distinct color
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Add linear fit with dashed line
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  ) +
  labs(x = "Gene level", y = "isoform level") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +  # Increase border line width
  coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(-7.5, 7.5))  # Set x and y limits


MRL_diff_NC_combine_with_AUC_diff <- MRL_diff_NC_combine_with_AUC_diff %>%
  group_by(gene_id.x) %>%
  mutate(highest_expression = if_else(expression == max(expression), 1, 0))

MRL_diff_NC_combine_with_AUC_diff <- MRL_diff_NC_combine_with_AUC_diff %>%
  group_by(gene_id.x) %>%
  mutate(highest_auc_diff_iso = auc_diff_iso[highest_expression == 1],
         diff_from_highest = auc_diff_iso - highest_auc_diff_iso)


non_highest_data <- MRL_diff_NC_combine_with_AUC_diff %>%
  filter(highest_expression == 0)

# Calculate counts and ratios
total_count <- nrow(non_highest_data)
count_greater_than_zero <- nrow(non_highest_data %>% filter(diff_from_highest > 0))
count_less_than_zero <- nrow(non_highest_data %>% filter(diff_from_highest < 0))
ratio_greater_than_zero <- count_greater_than_zero / total_count
ratio_less_than_zero <- count_less_than_zero / total_count

# Create the histogram plot
FigS2H_plot <- ggplot(non_highest_data, aes(x = diff_from_highest)) +
  geom_histogram(bins = 50, fill = "#E69F00", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  theme_bw() +
  labs(x = "Difference from Highest AUC_diff_iso_NC2h", y = "Count", title = "Histogram of Differences from Highest Expression Transcripts") +
  annotate("text", x = 0, y = max(table(cut(non_highest_data$diff_from_highest, breaks = 50))) * 0.75, 
           label = paste0("Ratio > 0: ", round(ratio_greater_than_zero, 3)), color = "blue", size = 5, hjust = 0.5) +
  annotate("text", x = min(non_highest_data$diff_from_highest) * 0.75, y = max(table(cut(non_highest_data$diff_from_highest, breaks = 50))), 
           label = paste0("Ratio < 0: ", round(ratio_less_than_zero, 3)), color = "blue", size = 5, hjust = 0)



library(scales) # For rescale function

# Define a function to calculate the weighted mean
weighted_mean <- function(x, w) {
  sum(x * w) / sum(w)
}

filtered_data <- MRL_diff_NC_combine_with_AUC_diff %>%
  filter(gene_name %in% c("ACTB", "EEF1B2", "GAPDH","PGK1", "PPIA", "RPL13A", "B2M",
                          "CCND1","CCND2","CCND3", "CDK4", "CDK6", "CDK2", 
                          "RPLP0", "RPLP1","RPLP2","RPL3","RPL4","RPL5","RPL6","RPL7","RPL8","RPL9","RPL10","RPL11","RPL12","RPL13"))

filtered_data <- filtered_data %>%
  group_by(gene_name) %>%
  mutate(scaled_expression = rescale(expression, to = c(0, 1)))

# Calculate weighted mean for each gene
weighted_means <- filtered_data %>%
  group_by(gene_name) %>%
  summarise(weighted_mean = weighted_mean(auc_diff_iso, scaled_expression))


FigS2G <- ggplot(filtered_data, aes(x = gene_name, y = auc_diff_iso, size = scaled_expression)) +
  geom_point(position = position_jitter(width = 0.2, height = 0), color = "#56B4E9") +
  geom_point(data = weighted_means, aes(x = gene_name, y = weighted_mean), color = "black", size = 5, shape = "-") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene Name", y = "Coupling", title = "Compare isoforms")

MRL_diff_NC_combine_with_AUC_diff <- MRL_diff_NC_combine_with_AUC_diff %>%
  group_by(gene_id.x) %>%
  mutate(highest_expression = if_else(expression == max(expression), 1, 0))

MRL_diff_NC_combine_with_AUC_diff <- MRL_diff_NC_combine_with_AUC_diff %>%
  group_by(gene_id.x) %>%
  mutate(highest_auc_diff_iso = auc_diff_iso[highest_expression == 1],
         diff_from_highest = auc_diff_iso - highest_auc_diff_iso)


non_highest_data <- MRL_diff_NC_combine_with_AUC_diff %>%
  filter(highest_expression == 0)



