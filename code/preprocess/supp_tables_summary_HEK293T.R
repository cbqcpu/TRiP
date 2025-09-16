source("code/polysome_seq_functions_20230703.R")
load("data/HEK293T_preprocessed.RData")
load("data/RNA_features_gene_level_20240418.RData")

merge_df_counts_select <- merge_df_counts[rowMeans(merge_df_counts[,-1])>1000,]

merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts$gene_id %in% merge_df_counts_select$gene_id,]

merge_df_counts_old_select <- data.frame(gene_id = merge_df_counts_select$gene_id, merge_df_counts_select[,-1] - merge_df_counts_new_select[,-1])

#MRL_total <- calculate_polysome_load(merge_df_counts_select, c("NC_30min", "NC_1h", "NC_2h"))
MRLs <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h", "STM2457", "YTHDF1", "YTHDF2", "YTHDF3"))
MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]

MRL_diff <- MRL_new-MRL_old
data_combined <- data.frame(
  NC_2h = c(MRL_old$NC_2h, MRL_new$NC_2h),
  Group = rep(c("MRL_old", "MRL_new"), 
              times = c(length(MRL_old$NC_2h), length(MRL_new$NC_2h)))
)

library(pracma)
time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new <- apply(MRL_new[,1:3], 1, calculate_auc)
auc_MRL_old <- apply(MRL_old[,1:3], 1, calculate_auc)

auc_diff <- auc_MRL_new - auc_MRL_old

colnames(MRL_new) <- paste0("MRL_new_", colnames(MRL_new))
colnames(MRL_old) <- paste0("MRL_old_", colnames(MRL_old))
colnames(MRL_diff) <- paste0("MRL_diff_", colnames(MRL_diff))

gene_id_and_name <- data.frame(gene_id = RNA_features_gene_level$gene_id, gene_name = RNA_features_gene_level$gene_name)

MRL_HEK293T <- cbind(MRL_new, MRL_old, MRL_diff)
MRL_HEK293T$auc_diff <- auc_diff
MRL_HEK293T$gene_id <- merge_df_counts_select$gene_id
MRL_HEK293T <- merge(MRL_HEK293T, gene_id_and_name, by = "gene_id")

write.csv(MRL_HEK293T, file = "data/HEK293T_supp.csv", quote = F, row.names = F)
write.csv(RNA_features_gene_level, file = "data/HEK293T_RNA_features_gene_level_supp.csv", quote = F, row.names = F)

