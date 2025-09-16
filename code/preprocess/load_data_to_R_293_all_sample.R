source("code/polysome_seq_functions_20230703.R")
# Directory path
dir_path_NC_30min <- "data/time_point/NC30min"
NC_30min <- read_counts(dir_path_NC_30min)
dir_path_NC_1h <- "data/time_point/NC1h"
NC_1h <- read_counts(dir_path_NC_1h)
dir_path_NC_2h <- "data/time_point/NC2h_new2"
NC_2h <- read_counts(dir_path_NC_2h)

save(NC_30min, NC_1h, NC_2h, file = "data/time_point/NC_30min_1h_2h.RData")

dir_path_EGI <- "data/Inhibitors/4EGI2h/"
EGI <- read_counts(dir_path_EGI)
dir_path_actD <- "data/Inhibitors/actD2h/"
actD <- read_counts(dir_path_actD)
dir_path_puro <- "data/Inhibitors/Puro2h/"
puro <- read_counts(dir_path_puro)
dir_path_STM2457x2h <- "data/Inhibitors/STM2457x2h"
STM2457 <- read_counts(dir_path_STM2457x2h)

save(EGI, actD, puro, STM2457, file = "data/Inhibitors/Inhibitors.RData")


# Load your data
# counts_matrix is a matrix of raw count data with genes in rows and samples in columns
# col_data is a data frame with sample information, such as condition
load("data/time_point/NC_30min_1h_2h.RData")
load("data/Inhibitors/Inhibitors.RData")
load("data/YTHDFs_KO/YTHDF_KO.RData")
# Load edgeR package
library(edgeR)
library(dplyr)

# Your counts_matrix preparation (from the code you provided)
counts_matrix <- cbind(NC_1h$`NC1h-1-1`$ReadCount,
                       NC_1h$`NC1h-2-1`$ReadCount,
                       NC_1h$`NC1h-3-1`$ReadCount,
                       NC_1h$`NC1h-1-2`$ReadCount,
                       NC_1h$`NC1h-2-2`$ReadCount,
                       NC_1h$`NC1h-3-2`$ReadCount,
                       NC_1h$`NC1h-1-3`$ReadCount,
                       NC_1h$`NC1h-2-3`$ReadCount,
                       NC_1h$`NC1h-3-3`$ReadCount,
                       NC_1h$`NC1h-1-4`$ReadCount,
                       NC_1h$`NC1h-2-4`$ReadCount,
                       NC_1h$`NC1h-3-4`$ReadCount,
                       NC_2h$`NC2h-1-1`$ReadCount,
                       NC_2h$`NC2h-2-1`$ReadCount,
                       NC_2h$`NC2h-3-1`$ReadCount,
                       NC_2h$`NC2h-1-2`$ReadCount,
                       NC_2h$`NC2h-2-2`$ReadCount,
                       NC_2h$`NC2h-3-2`$ReadCount,
                       NC_2h$`NC2h-1-3`$ReadCount,
                       NC_2h$`NC2h-2-3`$ReadCount,
                       NC_2h$`NC2h-3-3`$ReadCount,
                       NC_2h$`NC2h-1-4`$ReadCount,
                       NC_2h$`NC2h-2-4`$ReadCount,
                       NC_2h$`NC2h-3-4`$ReadCount,
                       NC_30min$`NC30min-1-1`$ReadCount,
                       NC_30min$`NC30min-5-1`$ReadCount,
                       NC_30min$`NC30min-6-1`$ReadCount,
                       NC_30min$`NC30min-1-2`$ReadCount,
                       NC_30min$`NC30min-5-2`$ReadCount,
                       NC_30min$`NC30min-6-2`$ReadCount,
                       NC_30min$`NC30min-1-3`$ReadCount,
                       NC_30min$`NC30min-5-3`$ReadCount,
                       NC_30min$`NC30min-6-3`$ReadCount,
                       NC_30min$`NC30min-1-4`$ReadCount,
                       NC_30min$`NC30min-5-4`$ReadCount,
                       NC_30min$`NC30min-6-4`$ReadCount,
                       actD$`actDx2h-1-1`$ReadCount,
                       actD$`actDx2h-2-1`$ReadCount,
                       actD$`actDx2h-3-1`$ReadCount,
                       actD$`actDx2h-1-2`$ReadCount,
                       actD$`actDx2h-2-2`$ReadCount,
                       actD$`actDx2h-3-2`$ReadCount,
                       actD$`actDx2h-1-3`$ReadCount,
                       actD$`actDx2h-2-3`$ReadCount,
                       actD$`actDx2h-3-3`$ReadCount,
                       actD$`actDx2h-1-4`$ReadCount,
                       actD$`actDx2h-2-4`$ReadCount,
                       actD$`actDx2h-3-4`$ReadCount,
                       EGI$`4EGIx2h-1-1`$ReadCount,
                       EGI$`4EGIx2h-2-1`$ReadCount,
                       EGI$`4EGIx2h-3-1`$ReadCount,
                       EGI$`4EGIx2h-1-2`$ReadCount,
                       EGI$`4EGIx2h-2-2`$ReadCount,
                       EGI$`4EGIx2h-3-2`$ReadCount,
                       EGI$`4EGIx2h-1-3`$ReadCount,
                       EGI$`4EGIx2h-2-3`$ReadCount,
                       EGI$`4EGIx2h-3-3`$ReadCount,
                       EGI$`4EGIx2h-1-4`$ReadCount,
                       EGI$`4EGIx2h-2-4`$ReadCount,
                       EGI$`4EGIx2h-3-4`$ReadCount,
                       puro$`Purox2h-1-1`$ReadCount,
                       puro$`Purox2h-4-1`$ReadCount,
                       puro$`Purox2h-5-1`$ReadCount,
                       puro$`Purox2h-1-2`$ReadCount,
                       puro$`Purox2h-4-2`$ReadCount,
                       puro$`Purox2h-5-2`$ReadCount,
                       puro$`Purox2h-1-3`$ReadCount,
                       puro$`Purox2h-4-3`$ReadCount,
                       puro$`Purox2h-5-3`$ReadCount,
                       puro$`Purox2h-1-4`$ReadCount,
                       puro$`Purox2h-4-4`$ReadCount,
                       puro$`Purox2h-5-4`$ReadCount,
                       STM2457$`STM2457x20uMx2h-1-1`$ReadCount,
                       STM2457$`STM2457x20uMx2h-3-1`$ReadCount,
                       STM2457$`STM2457x20uMx2h-4-1`$ReadCount,
                       STM2457$`STM2457x20uMx2h-1-2`$ReadCount,
                       STM2457$`STM2457x20uMx2h-3-2`$ReadCount,
                       STM2457$`STM2457x20uMx2h-4-2`$ReadCount,
                       STM2457$`STM2457x20uMx2h-1-3`$ReadCount,
                       STM2457$`STM2457x20uMx2h-3-3`$ReadCount,
                       STM2457$`STM2457x20uMx2h-4-3`$ReadCount,
                       STM2457$`STM2457x20uMx2h-1-4`$ReadCount,
                       STM2457$`STM2457x20uMx2h-3-4`$ReadCount,
                       STM2457$`STM2457x20uMx2h-4-4`$ReadCount,
                       YTHDF1$`YTHDF1-1-1`$ReadCount,
                       YTHDF1$`YTHDF1-2-1`$ReadCount,
                       YTHDF1$`YTHDF1-3-1`$ReadCount,
                       YTHDF1$`YTHDF1-1-2`$ReadCount,
                       YTHDF1$`YTHDF1-2-2`$ReadCount,
                       YTHDF1$`YTHDF1-3-2`$ReadCount,
                       YTHDF1$`YTHDF1-1-3`$ReadCount,
                       YTHDF1$`YTHDF1-2-3`$ReadCount,
                       YTHDF1$`YTHDF1-3-3`$ReadCount,
                       YTHDF1$`YTHDF1-1-4`$ReadCount,
                       YTHDF1$`YTHDF1-2-4`$ReadCount,
                       YTHDF1$`YTHDF1-3-4`$ReadCount,
                       YTHDF2$`YTHDF2-1-1`$ReadCount,
                       YTHDF2$`YTHDF2-2-1`$ReadCount,
                       YTHDF2$`YTHDF2-3-1`$ReadCount,
                       YTHDF2$`YTHDF2-1-2`$ReadCount,
                       YTHDF2$`YTHDF2-2-2`$ReadCount,
                       YTHDF2$`YTHDF2-3-2`$ReadCount,
                       YTHDF2$`YTHDF2-1-3`$ReadCount,
                       YTHDF2$`YTHDF2-2-3`$ReadCount,
                       YTHDF2$`YTHDF2-3-3`$ReadCount,
                       YTHDF2$`YTHDF2-1-4`$ReadCount,
                       YTHDF2$`YTHDF2-2-4`$ReadCount,
                       YTHDF2$`YTHDF2-3-4`$ReadCount,
                       YTHDF3$`YTHDF3-1-1`$ReadCount,
                       YTHDF3$`YTHDF3-2-1`$ReadCount,
                       YTHDF3$`YTHDF3-3-1`$ReadCount,
                       YTHDF3$`YTHDF3-1-2`$ReadCount,
                       YTHDF3$`YTHDF3-2-2`$ReadCount,
                       YTHDF3$`YTHDF3-3-2`$ReadCount,
                       YTHDF3$`YTHDF3-1-3`$ReadCount,
                       YTHDF3$`YTHDF3-2-3`$ReadCount,
                       YTHDF3$`YTHDF3-3-3`$ReadCount,
                       YTHDF3$`YTHDF3-1-4`$ReadCount,
                       YTHDF3$`YTHDF3-2-4`$ReadCount,
                       YTHDF3$`YTHDF3-3-4`$ReadCount)

# Assuming NC_30min$`THP1-NC-30min-1-1`$Name contains gene_id
gene_ids <- NC_1h$`NC1h-1-1`$Name

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
colnames(counts_matrix_merge) <- c("1h_frac1_rep1","1h_frac1_rep2","1h_frac1_rep3","1h_frac2_rep1","1h_frac2_rep2","1h_frac2_rep3","1h_frac3_rep1","1h_frac3_rep2","1h_frac3_rep3","1h_frac4_rep1","1h_frac4_rep2","1h_frac4_rep3",
                                   "2h_frac1_rep1","2h_frac1_rep2","2h_frac1_rep3","2h_frac2_rep1","2h_frac2_rep2","2h_frac2_rep3","2h_frac3_rep1","2h_frac3_rep2","2h_frac3_rep3","2h_frac4_rep1","2h_frac4_rep2","2h_frac4_rep3",
                                   "30min_frac1_rep1","30min_frac1_rep2","30min_frac1_rep3","30min_frac2_rep1","30min_frac2_rep2","30min_frac2_rep3","30min_frac3_rep1","30min_frac3_rep2","30min_frac3_rep3","30min_frac4_rep1","30min_frac4_rep2","30min_frac4_rep3",
                                   "actD_frac1_rep1","actD_frac1_rep2","actD_frac1_rep3","actD_frac2_rep1","actD_frac2_rep2","actD_frac2_rep3","actD_frac3_rep1","actD_frac3_rep2","actD_frac3_rep3","actD_frac4_rep1","actD_frac4_rep2","actD_frac4_rep3",
                                   "EGI_frac1_rep1","EGI_frac1_rep2","EGI_frac1_rep3","EGI_frac2_rep1","EGI_frac2_rep2","EGI_frac2_rep3","EGI_frac3_rep1","EGI_frac3_rep2","EGI_frac3_rep3","EGI_frac4_rep1","EGI_frac4_rep2","EGI_frac4_rep3",
                                   "puro_frac1_rep1","puro_frac1_rep2","puro_frac1_rep3","puro_frac2_rep1","puro_frac2_rep2","puro_frac2_rep3","puro_frac3_rep1","puro_frac3_rep2","puro_frac3_rep3","puro_frac4_rep1","puro_frac4_rep2","puro_frac4_rep3",
                                   "STM2457_frac1_rep1","STM2457_frac1_rep2","STM2457_frac1_rep3","STM2457_frac2_rep1","STM2457_frac2_rep2","STM2457_frac2_rep3","STM2457_frac3_rep1","STM2457_frac3_rep2","STM2457_frac3_rep3","STM2457_frac4_rep1","STM2457_frac4_rep2","STM2457_frac4_rep3",
                                   "YTHDF1_frac1_rep1","YTHDF1_frac1_rep2","YTHDF1_frac1_rep3","YTHDF1_frac2_rep1","YTHDF1_frac2_rep2","YTHDF1_frac2_rep3","YTHDF1_frac3_rep1","YTHDF1_frac3_rep2","YTHDF1_frac3_rep3","YTHDF1_frac4_rep1","YTHDF1_frac4_rep2","YTHDF1_frac4_rep3",
                                   "YTHDF2_frac1_rep1","YTHDF2_frac1_rep2","YTHDF2_frac1_rep3","YTHDF2_frac2_rep1","YTHDF2_frac2_rep2","YTHDF2_frac2_rep3","YTHDF2_frac3_rep1","YTHDF2_frac3_rep2","YTHDF2_frac3_rep3","YTHDF2_frac4_rep1","YTHDF2_frac4_rep2","YTHDF2_frac4_rep3",
                                   "YTHDF3_frac1_rep1","YTHDF3_frac1_rep2","YTHDF3_frac1_rep3","YTHDF3_frac2_rep1","YTHDF3_frac2_rep2","YTHDF3_frac2_rep3","YTHDF3_frac3_rep1","YTHDF3_frac3_rep2","YTHDF3_frac3_rep3","YTHDF3_frac4_rep1","YTHDF3_frac4_rep2","YTHDF3_frac4_rep3")

group <- factor(c(rep("1h_frac1", 3), rep("1h_frac2", 3), rep("1h_frac3", 3), rep("1h_frac4", 3),
                  rep("2h_frac1", 3), rep("2h_frac2", 3), rep("2h_frac3", 3), rep("2h_frac4", 3),
                  rep("30min_frac1", 3), rep("30min_frac2", 3), rep("30min_frac3", 3), rep("30min_frac4", 3),
                  rep("actD_frac1", 3), rep("actD_frac2", 3), rep("actD_frac3", 3), rep("actD_frac4", 3),
                  rep("EGI_frac1", 3), rep("EGI_frac2", 3), rep("EGI_frac3", 3), rep("EGI_frac4", 3),
                  rep("puro_frac1", 3), rep("puro_frac2", 3), rep("puro_frac3", 3), rep("puro_frac4", 3),
                  rep("STM2457_frac1", 3), rep("STM2457_frac2", 3), rep("STM2457_frac3", 3), rep("STM2457_frac4", 3),
                  rep("YTHDF1_frac1", 3), rep("YTHDF1_frac2", 3), rep("YTHDF1_frac3", 3), rep("YTHDF1_frac4", 3),
                  rep("YTHDF2_frac1", 3), rep("YTHDF2_frac2", 3), rep("YTHDF2_frac3", 3), rep("YTHDF2_frac4", 3),
                  rep("YTHDF3_frac1", 3), rep("YTHDF3_frac2", 3), rep("YTHDF3_frac3", 3), rep("YTHDF3_frac4", 3)))

########### re-run start here ############
#load("data/THP1/filter_genes_workspace.RData")

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
#write.csv(genes_to_delete, "genes_to_delete_293T.csv", row.names = FALSE)
#write.csv(genes_to_keep, "genes_to_keep_293T.csv", row.names = FALSE)
#write.table(genes_logCPM, file = "gene_logCPM_293T.csv", sep = ",", quote = FALSE, col.names = TRUE)

# Calculate deletion ratio
deletion_ratio <- length(genes_to_delete) / (length(genes_to_delete) + length(genes_to_keep))
print(paste("Deletion ratio:", deletion_ratio))

#rm(LPS_2h,LPS_12h,LPS_24h,NC_1h,NC_2h,NC_30min)
#save.image(file = "data/THP1/filter_genes_workspace.RData")

source("code/polysome_seq_functions_20230703.R")

NC1h_yeast_all_counts <- read_csv("data/time_point/NC1h/NC1h_yeast_all_counts.csv")
NC2h_yeast_all_counts <- read_csv("data/time_point/NC2h_new2/NC2h_yeast_all_counts.csv")
NC30min_yeast_all_counts <- read_csv("data/time_point/NC30min/NC30min_yeast_all_counts.csv")
EGI_yeast_all_counts <- read_csv("data/Inhibitors/4EGI2h/4EGI2h_yeast_all_counts.csv")
actD_yeast_all_counts <- read_csv("data/Inhibitors/actD2h/actD2h_yeast_all_counts.csv")
puro_yeast_all_counts <- read_csv("data/Inhibitors/Puro2h/Puro2h_yeast_all_counts.csv")
STM2457_yeast_all_counts <- read_csv("data/Inhibitors/STM2457x2h/STM2457x2h_yeast_all_counts.csv")
YTHDF1_yeast_all_counts <- read_csv("data/YTHDFs_KO/YTHDF1/YTHDF1_yeast_all_counts.csv")
YTHDF2_yeast_all_counts <- read_csv("data/YTHDFs_KO/YTHDF2/YTHDF2_yeast_all_counts.csv")
YTHDF3_yeast_all_counts <- read_csv("data/YTHDFs_KO/YTHDF3/YTHDF3_yeast_all_counts.csv")

NC_30min_merge <- merge_reps_with_cv(NC_30min, "NC30min", '30min', c("1","5","6"), NC30min_yeast_all_counts)
NC_1h_merge <- merge_reps_with_cv(NC_1h, "NC1h", '1h', c("1","2","3"), NC1h_yeast_all_counts)
NC_2h_merge <- merge_reps_with_cv(NC_2h, "NC2h", '2h', c("1","2","3"), NC2h_yeast_all_counts)
EGI_merge <- merge_reps_with_cv(EGI, "4EGIx2h", '4EGIx2h', c("1","2","3"), EGI_yeast_all_counts)
actD_merge <- merge_reps_with_cv(actD, "actDx2h", 'actDx2h', c("1","2","3"), actD_yeast_all_counts)
puro_merge <- merge_reps_with_cv(puro, "Purox2h", 'Purox2h', c("1","4","5"), puro_yeast_all_counts)
STM2457_merge <- merge_reps_with_cv(STM2457, "STM2457x20uMx2h", 'STM2457x20uMx2h', c("1","3","4"), STM2457_yeast_all_counts)
YTHDF1_merge <- merge_reps_with_cv(YTHDF1, "YTHDF1", 'YTHDF1', c("1","2","3"), YTHDF1_yeast_all_counts)
YTHDF2_merge <- merge_reps_with_cv(YTHDF2, "YTHDF2", 'YTHDF2', c("1","2","3"), YTHDF2_yeast_all_counts)
YTHDF3_merge <- merge_reps_with_cv(YTHDF3, "YTHDF3", 'YTHDF3', c("1","2","3"), YTHDF3_yeast_all_counts)



TC_rate <- data.frame(NC_30min_merge$gene_id$Name, NC_30min_merge$TC_rate_merge, NC_1h_merge$TC_rate_merge, NC_2h_merge$TC_rate_merge, EGI_merge$TC_rate_merge, actD_merge$TC_rate_merge, puro_merge$TC_rate_merge, STM2457_merge$TC_rate_merge, YTHDF1_merge$TC_rate_merge, YTHDF2_merge$TC_rate_merge, YTHDF3_merge$TC_rate_merge)
colnames(TC_rate) <- c("gene_id", paste0("NC_30min_", colnames(NC_30min_merge$TC_rate_merge)),
                                paste0("NC_1h_", colnames(NC_1h_merge$TC_rate_merge)),
                                paste0("NC_2h_", colnames(NC_2h_merge$TC_rate_merge)),
                                paste0("EGI_", colnames(EGI_merge$TC_rate_merge)),
                                paste0("actD_", colnames(actD_merge$TC_rate_merge)),
                                paste0("puro_", colnames(puro_merge$TC_rate_merge)),
                                paste0("STM2457_", colnames(STM2457_merge$TC_rate_merge)),
                                paste0("YTHDF1_", colnames(YTHDF1_merge$TC_rate_merge)),
                                paste0("YTHDF2_", colnames(YTHDF2_merge$TC_rate_merge)),
                                paste0("YTHDF3_", colnames(YTHDF3_merge$TC_rate_merge)))

counts <- data.frame(NC_30min_merge$gene_id$Name, NC_30min_merge$data_norm_counts, NC_1h_merge$data_norm_counts, NC_2h_merge$data_norm_counts, EGI_merge$data_norm_counts, actD_merge$data_norm_counts, puro_merge$data_norm_counts, STM2457_merge$data_norm_counts, YTHDF1_merge$data_norm_counts, YTHDF2_merge$data_norm_counts, YTHDF3_merge$data_norm_counts)
colnames(counts) <- c("gene_id", paste0("NC_30min_", colnames(NC_30min_merge$data_norm_counts)),
                                 paste0("NC_1h_", colnames(NC_1h_merge$data_norm_counts)),
                                 paste0("NC_2h_", colnames(NC_2h_merge$data_norm_counts)),
                                 paste0("EGI_", colnames(EGI_merge$data_norm_counts)),
                                 paste0("actD_", colnames(actD_merge$data_norm_counts)),
                                 paste0("puro_", colnames(puro_merge$data_norm_counts)),
                                 paste0("STM2457_", colnames(STM2457_merge$data_norm_counts)),
                                 paste0("YTHDF1_", colnames(YTHDF1_merge$data_norm_counts)),
                                 paste0("YTHDF2_", colnames(YTHDF2_merge$data_norm_counts)),
                                 paste0("YTHDF3_", colnames(YTHDF3_merge$data_norm_counts)))


new_counts <- data.frame(NC_30min_merge$gene_id$Name, NC_30min_merge$data_raw_counts_new_norm, NC_1h_merge$data_raw_counts_new_norm, NC_2h_merge$data_raw_counts_new_norm, EGI_merge$data_raw_counts_new_norm, actD_merge$data_raw_counts_new_norm, puro_merge$data_raw_counts_new_norm, STM2457_merge$data_raw_counts_new_norm, YTHDF1_merge$data_raw_counts_new_norm, YTHDF2_merge$data_raw_counts_new_norm, YTHDF3_merge$data_raw_counts_new_norm)
colnames(new_counts) <- c("gene_id", paste0("NC_30min_", colnames(NC_30min_merge$data_raw_counts_new_norm)),
                                paste0("NC_1h_", colnames(NC_1h_merge$data_raw_counts_new_norm)),
                                paste0("NC_2h_", colnames(NC_2h_merge$data_raw_counts_new_norm)),
                                paste0("EGI_", colnames(EGI_merge$data_raw_counts_new_norm)),
                                paste0("actD_", colnames(actD_merge$data_raw_counts_new_norm)),
                                paste0("puro_", colnames(puro_merge$data_raw_counts_new_norm)),
                                paste0("STM2457_", colnames(STM2457_merge$data_raw_counts_new_norm)),
                                paste0("YTHDF1_", colnames(YTHDF1_merge$data_raw_counts_new_norm)),
                                paste0("YTHDF2_", colnames(YTHDF2_merge$data_raw_counts_new_norm)),
                                paste0("YTHDF3_", colnames(YTHDF3_merge$data_raw_counts_new_norm)))

merge_df_TC_rate <- aggregate(. ~gene_id, data = TC_rate, mean)
merge_df_counts <- aggregate(. ~gene_id, data = counts, sum)
merge_df_counts_new <- aggregate(. ~gene_id, data = new_counts, sum)

save(merge_df_TC_rate, merge_df_counts, merge_df_counts_new, genes_logCPM, file = "data/HEK293T_preprocessed.RData")

merge_df_new_percentage <- data.frame(gene_id = merge_df_counts_new$gene_id, merge_df_counts_new[,-1]/merge_df_counts[,-1])
merge_df_new_percentage <- merge_df_new_percentage %>%
  filter(!if_any(everything(), is.nan))

# First, find the intersection of gene IDs to ensure they exist in your dataframe
common_genes <- intersect(merge_df_new_percentage$gene_id, genes_to_keep)

merge_df_TC_rate_select <- merge_df_TC_rate[merge_df_TC_rate$gene_id %in% common_genes, ]
merge_df_counts_select <- merge_df_counts[merge_df_counts$gene_id %in% common_genes, ]
merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts_new$gene_id %in% common_genes, ]
merge_df_new_percentage_select <- merge_df_new_percentage[merge_df_new_percentage$gene_id %in% common_genes, ]



p1 <- plot_violin(merge_df_TC_rate_select, "NC_30min")
p2 <- plot_violin(merge_df_TC_rate_select, "NC_1h")
p3 <- plot_violin(merge_df_TC_rate_select, "NC_2h")
p4 <- plot_violin(merge_df_TC_rate_select, "EGI")
p5 <- plot_violin(merge_df_TC_rate_select, "actD")
p6 <- plot_violin(merge_df_TC_rate_select, "puro")
p7 <- plot_violin(merge_df_TC_rate_select, "STM2457")
p8 <- plot_violin(merge_df_TC_rate_select, "YTHDF1")
p9 <- plot_violin(merge_df_TC_rate_select, "YTHDF2")
p10 <- plot_violin(merge_df_TC_rate_select, "YTHDF3")

# Arrange the three plots horizontally
#pdf(file = "results/Figure1/timepoints.pdf", width = 10, height = 3)
print(plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 3))
#dev.off()

p1 <- plot_violin_new_ratio(merge_df_new_percentage_select, "NC_30min")
p2 <- plot_violin_new_ratio(merge_df_new_percentage_select, "NC_1h")
p3 <- plot_violin_new_ratio(merge_df_new_percentage_select, "NC_2h")
p4 <- plot_violin_new_ratio(merge_df_new_percentage_select, "EGI")
p5 <- plot_violin_new_ratio(merge_df_new_percentage_select, "actD")
p6 <- plot_violin_new_ratio(merge_df_new_percentage_select, "puro")
p7 <- plot_violin_new_ratio(merge_df_new_percentage_select, "STM2457")
p8 <- plot_violin_new_ratio(merge_df_new_percentage_select, "YTHDF1")
p9 <- plot_violin_new_ratio(merge_df_new_percentage_select, "YTHDF2")
p10 <- plot_violin_new_ratio(merge_df_new_percentage_select, "YTHDF3")

# Arrange the three plots horizontally
#pdf(file = "results/Figure1/timepoints_new_mRNA_ratio.pdf", width = 10, height = 3)
print(plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, ncol = 3))
#dev.off()

#save(merge_df_TC_rate_select, merge_df_counts_select, merge_df_counts_new_select, merge_df_new_percentage_select, file = "data/THP1_filtered_gene_20231109.RData")




