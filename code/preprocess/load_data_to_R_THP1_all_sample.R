source("code/polysome_seq_functions_20230703.R")
dir_path_NC_30min <- "data/THP1/NC_30min"
NC_30min <- read_counts(dir_path_NC_30min)
dir_path_NC_1h <- "data/THP1/NC_1h"
NC_1h <- read_counts(dir_path_NC_1h)
dir_path_NC_2h <- "data/THP1/NC_2h"
NC_2h <- read_counts(dir_path_NC_2h)

save(NC_30min, NC_1h, NC_2h, file = "data/THP1/THP1_NC_30min_1h_2h.RData")

dir_path_NC_2h <- "data/THP1/NC_2h/"
NC_2h <- read_counts(dir_path_NC_2h)
dir_path_LPS_2h <- "data/THP1/LPS_2h/"
LPS_2h <- read_counts(dir_path_LPS_2h)
dir_path_LPS_12h <- "data/THP1/LPS_12h/"
LPS_12h <- read_counts(dir_path_LPS_12h)
dir_path_LPS_24h <- "data/THP1/LPS_24h/"
LPS_24h <- read_counts(dir_path_LPS_24h)

save(NC_2h, LPS_2h, LPS_12h, LPS_24h, file = "data/THP1/THP1_NC_LPS.RData")

# Load your data
# counts_matrix is a matrix of raw count data with genes in rows and samples in columns
# col_data is a data frame with sample information, such as condition
load("data/THP1/THP1_NC_30min_1h_2h.RData")
load("data/THP1/THP1_NC_LPS.RData")
# Load edgeR package
library(edgeR)
library(dplyr)

# Your counts_matrix preparation (from the code you provided)
counts_matrix <- cbind(NC_30min$`THP1-NC-30min-1-1`$ReadCount,
                       NC_30min$`THP1-NC-30min-2-1`$ReadCount,
                       NC_30min$`THP1-NC-30min-3-1`$ReadCount,
                       NC_30min$`THP1-NC-30min-1-2`$ReadCount,
                       NC_30min$`THP1-NC-30min-2-2`$ReadCount,
                       NC_30min$`THP1-NC-30min-3-2`$ReadCount,
                       NC_30min$`THP1-NC-30min-1-3`$ReadCount,
                       NC_30min$`THP1-NC-30min-2-3`$ReadCount,
                       NC_30min$`THP1-NC-30min-3-3`$ReadCount,
                       NC_30min$`THP1-NC-30min-1-4`$ReadCount,
                       NC_30min$`THP1-NC-30min-2-4`$ReadCount,
                       NC_30min$`THP1-NC-30min-3-4`$ReadCount,
                       NC_1h$`THP1-NC-1h-1-1`$ReadCount,
                       NC_1h$`THP1-NC-1h-2-1`$ReadCount,
                       NC_1h$`THP1-NC-1h-3-1`$ReadCount,
                       NC_1h$`THP1-NC-1h-1-2`$ReadCount,
                       NC_1h$`THP1-NC-1h-2-2`$ReadCount,
                       NC_1h$`THP1-NC-1h-3-2`$ReadCount,
                       NC_1h$`THP1-NC-1h-1-3`$ReadCount,
                       NC_1h$`THP1-NC-1h-2-3`$ReadCount,
                       NC_1h$`THP1-NC-1h-3-3`$ReadCount,
                       NC_1h$`THP1-NC-1h-1-4`$ReadCount,
                       NC_1h$`THP1-NC-1h-2-4`$ReadCount,
                       NC_1h$`THP1-NC-1h-3-4`$ReadCount,
                       NC_2h$`THP1-NC-2-1`$ReadCount,
                       NC_2h$`THP1-NC-3-1`$ReadCount,
                       NC_2h$`THP1-NC-4-1`$ReadCount,
                       NC_2h$`THP1-NC-2-2`$ReadCount,
                       NC_2h$`THP1-NC-3-2`$ReadCount,
                       NC_2h$`THP1-NC-4-2`$ReadCount,
                       NC_2h$`THP1-NC-2-3`$ReadCount,
                       NC_2h$`THP1-NC-3-3`$ReadCount,
                       NC_2h$`THP1-NC-4-3`$ReadCount,
                       NC_2h$`THP1-NC-2-4`$ReadCount,
                       NC_2h$`THP1-NC-3-4`$ReadCount,
                       NC_2h$`THP1-NC-4-4`$ReadCount,
                       LPS_2h$`THP1-2h-1-1`$ReadCount,
                       LPS_2h$`THP1-2h-3-1`$ReadCount,
                       LPS_2h$`THP1-2h-4-1`$ReadCount,
                       LPS_2h$`THP1-2h-1-2`$ReadCount,
                       LPS_2h$`THP1-2h-3-2`$ReadCount,
                       LPS_2h$`THP1-2h-4-2`$ReadCount,
                       LPS_2h$`THP1-2h-1-3`$ReadCount,
                       LPS_2h$`THP1-2h-3-3`$ReadCount,
                       LPS_2h$`THP1-2h-4-3`$ReadCount,
                       LPS_2h$`THP1-2h-1-4`$ReadCount,
                       LPS_2h$`THP1-2h-3-4`$ReadCount,
                       LPS_2h$`THP1-2h-4-4`$ReadCount,
                       LPS_12h$`THP1-12h-1-1`$ReadCount,
                       LPS_12h$`THP1-12h-3-1`$ReadCount,
                       LPS_12h$`THP1-12h-4-1`$ReadCount,
                       LPS_12h$`THP1-12h-1-2`$ReadCount,
                       LPS_12h$`THP1-12h-3-2`$ReadCount,
                       LPS_12h$`THP1-12h-4-2`$ReadCount,
                       LPS_12h$`THP1-12h-1-3`$ReadCount,
                       LPS_12h$`THP1-12h-3-3`$ReadCount,
                       LPS_12h$`THP1-12h-4-3`$ReadCount,
                       LPS_12h$`THP1-12h-1-4`$ReadCount,
                       LPS_12h$`THP1-12h-3-4`$ReadCount,
                       LPS_12h$`THP1-12h-4-4`$ReadCount,
                       LPS_24h$`THP1-24h-3-1`$ReadCount,
                       LPS_24h$`THP1-24h-4-1`$ReadCount,
                       LPS_24h$`THP1-24h-5-1`$ReadCount,
                       LPS_24h$`THP1-24h-3-2`$ReadCount,
                       LPS_24h$`THP1-24h-4-2`$ReadCount,
                       LPS_24h$`THP1-24h-5-2`$ReadCount,
                       LPS_24h$`THP1-24h-3-3`$ReadCount,
                       LPS_24h$`THP1-24h-4-3`$ReadCount,
                       LPS_24h$`THP1-24h-5-3`$ReadCount,
                       LPS_24h$`THP1-24h-3-4`$ReadCount,
                       LPS_24h$`THP1-24h-4-4`$ReadCount,
                       LPS_24h$`THP1-24h-5-4`$ReadCount)

# Assuming NC_30min$`THP1-NC-30min-1-1`$Name contains gene_id
gene_ids <- NC_30min$`THP1-NC-30min-1-1`$Name

counts_matrix <- data.frame(counts_matrix)
counts_matrix$gene_id <- gene_ids

counts_matrix_merge <- counts_matrix %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

gene_ids_unique <- counts_matrix_merge$gene_id
counts_matrix_merge <- counts_matrix_merge[,-1]
counts_matrix_merge <- as.matrix(counts_matrix_merge)
rownames(counts_matrix_merge) <- gene_ids_unique

# Assign column names to counts_matrix
colnames(counts_matrix_merge) <- c("30min_frac1_rep1", "30min_frac1_rep2", "30min_frac1_rep3", "30min_frac2_rep1", "30min_frac2_rep2", "30min_frac2_rep3", "30min_frac3_rep1", "30min_frac3_rep2", "30min_frac3_rep3", "30min_frac4_rep1", "30min_frac4_rep2", "30min_frac4_rep3",
                             "1h_frac1_rep1", "1h_frac1_rep2", "1h_frac1_rep3", "1h_frac2_rep1", "1h_frac2_rep2", "1h_frac2_rep3", "1h_frac3_rep1", "1h_frac3_rep2", "1h_frac3_rep3", "1h_frac4_rep1", "1h_frac4_rep2", "1h_frac4_rep3",
                             "2h_frac1_rep1", "2h_frac1_rep2", "2h_frac1_rep3", "2h_frac2_rep1", "2h_frac2_rep2", "2h_frac2_rep3", "2h_frac3_rep1", "2h_frac3_rep2", "2h_frac3_rep3", "2h_frac4_rep1", "2h_frac4_rep2", "2h_frac4_rep3",
                             "LPS_2h_frac1_rep1", "LPS_2h_frac1_rep2", "LPS_2h_frac1_rep3", "LPS_2h_frac2_rep1", "LPS_2h_frac2_rep2", "LPS_2h_frac2_rep3", "LPS_2h_frac3_rep1", "LPS_2h_frac3_rep2", "LPS_2h_frac3_rep3", "LPS_2h_frac4_rep1", "LPS_2h_frac4_rep2", "LPS_2h_frac4_rep3",
                             "LPS_12h_frac1_rep1", "LPS_12h_frac1_rep2", "LPS_12h_frac1_rep3", "LPS_12h_frac2_rep1", "LPS_12h_frac2_rep2", "LPS_12h_frac2_rep3", "LPS_12h_frac3_rep1", "LPS_12h_frac3_rep2", "LPS_12h_frac3_rep3", "LPS_12h_frac4_rep1", "LPS_12h_frac4_rep2", "LPS_12h_frac4_rep3",
                             "LPS_24h_frac1_rep1", "LPS_24h_frac1_rep2", "LPS_24h_frac1_rep3", "LPS_24h_frac2_rep1", "LPS_24h_frac2_rep2", "LPS_24h_frac2_rep3", "LPS_24h_frac3_rep1", "LPS_24h_frac3_rep2", "LPS_24h_frac3_rep3", "LPS_24h_frac4_rep1", "LPS_24h_frac4_rep2", "LPS_24h_frac4_rep3")

group <- factor(c(rep("30min_frac1", 3), rep("30min_frac2", 3), rep("30min_frac3", 3), rep("30min_frac4", 3),
                  rep("1h_frac1", 3), rep("1h_frac2", 3), rep("1h_frac3", 3), rep("1h_frac4", 3),
                  rep("2h_frac1", 3), rep("2h_frac2", 3), rep("2h_frac3", 3), rep("2h_frac4", 3),
                  rep("LPS_2h_frac1", 3), rep("LPS_2h_frac2", 3), rep("LPS_2h_frac3", 3), rep("LPS_2h_frac4", 3),
                  rep("LPS_12h_frac1", 3), rep("LPS_12h_frac2", 3), rep("LPS_12h_frac3", 3), rep("LPS_12h_frac4", 3),
                  rep("LPS_24h_frac1", 3), rep("LPS_24h_frac2", 3), rep("LPS_24h_frac3", 3), rep("LPS_24h_frac4", 3)))


# Create a DGEList object
dge <- DGEList(counts = counts_matrix_merge, group = group)

# Calculate normalization factors
dge <- calcNormFactors(dge)

# Estimate common dispersion
dge <- estimateCommonDisp(dge)

# Estimate tagwise dispersion
dge <- estimateTagwiseDisp(dge)

# Fit gene-wise models
dge <- glmFit(dge)

# Conduct likelihood ratio tests
dge <- glmLRT(dge)

# Now extract topTags
tt <- topTags(dge, n=Inf)


# Now you can plot
plotSmear(dge, low=FALSE) # Plot the initial MA plot

# Extract logFC and average logCPM from the topTags results
df <- as.data.frame(tt)

# Identify genes with low and high logCPM
low_logCPM_genes <- df$logCPM < -4
high_logCPM_genes <- df$logCPM >= -4

# Add points for genes with low logCPM in blue
with(df[low_logCPM_genes, ], points(logCPM, logFC, pch=20, col="blue"))

# Add points for genes with high logCPM in red
with(df[high_logCPM_genes, ], points(logCPM, logFC, pch=20, col="red"))

# Identify genes to delete and to keep
genes_to_delete <- rownames(dge)[df$logCPM < -4]
genes_to_keep <- rownames(dge)[df$logCPM >= -4]

genes_logCPM <- df
genes_logCPM$gene_id <- rownames(genes_logCPM)

# Save gene lists to files
#write.csv(genes_to_delete, "genes_to_delete_THP1.csv", row.names = FALSE)
#write.csv(genes_to_keep, "genes_to_keep_THP1.csv", row.names = FALSE)
#write.table(genes_logCPM, file = "gene_logCPM_THP1.csv", sep = ",", quote = FALSE, col.names = TRUE)

# Calculate deletion ratio
deletion_ratio <- length(genes_to_delete) / (length(genes_to_delete) + length(genes_to_keep))
print(paste("Deletion ratio:", deletion_ratio))


source("code/polysome_seq_functions_20230703.R")

NC1h_yeast_all_counts <- read_csv("data/THP1/NC_1h/THP1-NC-1h_yeast_all_counts.csv")
NC2h_yeast_all_counts <- read_csv("data/THP1/NC_2h/THP1-NC_yeast_all_counts.csv")
NC30min_yeast_all_counts <- read_csv("data/THP1/NC_30min/THP1-NC-30min_yeast_all_counts.csv")
LPS_2h_yeast_all_counts <- read_csv("data/THP1/LPS_2h/THP1-2h_yeast_all_counts.csv")
LPS_12h_yeast_all_counts <- read_csv("data/THP1/LPS_12h/THP1-12h_yeast_all_counts.csv")
LPS_24h_yeast_all_counts <- read_csv("data/THP1/LPS_24h/THP1-24h_yeast_all_counts.csv")

NC_30min_merge <- merge_reps_with_cv(NC_30min, "THP1-NC-30min", '30min', c("1","2","3"), NC30min_yeast_all_counts)
NC_1h_merge <- merge_reps_with_cv(NC_1h, "THP1-NC-1h", '1h', c("1","2","3"), NC1h_yeast_all_counts)
NC_2h_merge <- merge_reps_with_cv(NC_2h, "THP1-NC", '2h', c("2","3","4"), NC2h_yeast_all_counts)
LPS_2h_merge <- merge_reps_with_cv(LPS_2h, "THP1-2h", 'THP1-2h', c("1","3","4"), LPS_2h_yeast_all_counts)
LPS_12h_merge <- merge_reps_with_cv(LPS_12h, "THP1-12h", 'THP1-12h', c("1","3","4"), LPS_12h_yeast_all_counts)
LPS_24h_merge <- merge_reps_with_cv(LPS_24h, "THP1-24h", 'THP1-24h', c("3","4","5"), LPS_24h_yeast_all_counts)


TC_rate <- data.frame(NC_30min_merge$gene_id$Name, NC_30min_merge$TC_rate_merge, NC_1h_merge$TC_rate_merge, NC_2h_merge$TC_rate_merge, LPS_2h_merge$TC_rate_merge, LPS_12h_merge$TC_rate_merge, LPS_24h_merge$TC_rate_merge)
colnames(TC_rate) <- c("gene_id", paste0("NC_30min_", colnames(NC_30min_merge$TC_rate_merge)),
                                paste0("NC_1h_", colnames(NC_1h_merge$TC_rate_merge)),
                                paste0("NC_2h_", colnames(NC_2h_merge$TC_rate_merge)),
                                paste0("LPS_2h_", colnames(LPS_2h_merge$TC_rate_merge)),
                                paste0("LPS_12h_", colnames(LPS_12h_merge$TC_rate_merge)),
                                paste0("LPS_24h_", colnames(LPS_24h_merge$TC_rate_merge)))

counts <- data.frame(NC_30min_merge$gene_id$Name, NC_30min_merge$data_norm_counts, NC_1h_merge$data_norm_counts, NC_2h_merge$data_norm_counts, LPS_2h_merge$data_norm_counts, LPS_12h_merge$data_norm_counts, LPS_24h_merge$data_norm_counts)
colnames(counts) <- c("gene_id", paste0("NC_30min_", colnames(NC_30min_merge$data_norm_counts)),
                               paste0("NC_1h_", colnames(NC_1h_merge$data_norm_counts)),
                               paste0("NC_2h_", colnames(NC_2h_merge$data_norm_counts)),
                               paste0("LPS_2h_", colnames(LPS_2h_merge$data_norm_counts)),
                               paste0("LPS_12h_", colnames(LPS_12h_merge$data_norm_counts)),
                               paste0("LPS_24h_", colnames(LPS_24h_merge$data_norm_counts)))

new_counts <- data.frame(NC_30min_merge$gene_id$Name, NC_30min_merge$data_raw_counts_new_norm, NC_1h_merge$data_raw_counts_new_norm, NC_2h_merge$data_raw_counts_new_norm, LPS_2h_merge$data_raw_counts_new_norm, LPS_12h_merge$data_raw_counts_new_norm, LPS_24h_merge$data_raw_counts_new_norm)
colnames(new_counts) <- c("gene_id", paste0("NC_30min_", colnames(NC_30min_merge$data_raw_counts_new_norm)),
                                   paste0("NC_1h_", colnames(NC_1h_merge$data_raw_counts_new_norm)),
                                   paste0("NC_2h_", colnames(NC_2h_merge$data_raw_counts_new_norm)),
                                   paste0("LPS_2h_", colnames(LPS_2h_merge$data_raw_counts_new_norm)),
                                   paste0("LPS_12h_", colnames(LPS_12h_merge$data_raw_counts_new_norm)),
                                   paste0("LPS_24h_", colnames(LPS_24h_merge$data_raw_counts_new_norm)))

merge_df_TC_rate <- aggregate(. ~gene_id, data = TC_rate, sum)
merge_df_counts <- aggregate(. ~gene_id, data = counts, sum)
merge_df_counts_new <- aggregate(. ~gene_id, data = new_counts, sum)


save(merge_df_TC_rate, merge_df_counts, merge_df_counts_new, genes_logCPM, file = "data/THP1_preprocessed.RData")

#merge_df_new_percentage <- data.frame(gene_id = merge_df_counts_new$gene_id, merge_df_counts_new[,-1]/merge_df_counts[,-1])
#merge_df_new_percentage <- merge_df_new_percentage %>%
#  filter(!if_any(everything(), is.nan))

# First, find the intersection of gene IDs to ensure they exist in your dataframe
#common_genes <- intersect(merge_df_new_percentage$gene_id, genes_to_keep)

#merge_df_TC_rate_select <- merge_df_TC_rate[merge_df_TC_rate$gene_id %in% common_genes, ]
#merge_df_counts_select <- merge_df_counts[merge_df_counts$gene_id %in% common_genes, ]
#merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts_new$gene_id %in% common_genes, ]
#merge_df_new_percentage_select <- merge_df_new_percentage[merge_df_new_percentage$gene_id %in% common_genes, ]




