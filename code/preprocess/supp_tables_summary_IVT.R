#from /storage/disk3/cbq/NGS_data/20240608_DP_lib_coupling/

library(pracma)
library(gridExtra)
library(ggplot2)
library(readxl)
library(scales)
source("code/polysome_seq_functions_20230703.R")
egfp_unmod_1 <- read_csv("data/GSM3130435_egfp_unmod_1.csv")
aligned_sequence <- read_csv("data/aligned_sequence.csv")
egfp_unmod_1_RL <- egfp_unmod_1[,c("utr", "rl")]
colnames(egfp_unmod_1_RL) <- c("sequence", "RL")
aligned_sequence_value <- aligned_sequence[,5:ncol(aligned_sequence)]
aligned_sequence_value_normalized <- aligned_sequence_value %>%
  mutate(across(everything(), ~ . / last(.)))

aligned_sequence_select <- aligned_sequence_value_normalized[aligned_sequence$file_tag %in% c("test", "mRNA_1273", "BNT-162b2"),]
lib_1h_rep1 <- aligned_sequence_select[,c("1h-1-1", "1h-1-2", "1h-1-3", "1h-1-4")]
lib_1h_rep2 <- aligned_sequence_select[,c("1h-2-1", "1h-2-2", "1h-2-3", "1h-2-4")]
lib_1h_rep3 <- aligned_sequence_select[,c("1h-3-1", "1h-3-2", "1h-3-3", "1h-3-4")]
colnames(lib_1h_rep1) <- c("1h_frac1", "1h_frac2", "1h_frac3", "1h_frac4")
colnames(lib_1h_rep2) <- c("1h_frac1", "1h_frac2", "1h_frac3", "1h_frac4")
colnames(lib_1h_rep3) <- c("1h_frac1", "1h_frac2", "1h_frac3", "1h_frac4")
lib_1h_rep1_MRL <- calculate_polysome_load_old(lib_1h_rep1, "1h")
lib_1h_rep2_MRL <- calculate_polysome_load_old(lib_1h_rep2, "1h")
lib_1h_rep3_MRL <- calculate_polysome_load_old(lib_1h_rep3, "1h")
lib_1h_rep1_mRNA <- rowSums(lib_1h_rep1)
lib_1h_rep2_mRNA <- rowSums(lib_1h_rep2)
lib_1h_rep3_mRNA <- rowSums(lib_1h_rep3)

lib_2h_rep1 <- aligned_sequence_select[,c("2h-1-1", "2h-1-2", "2h-1-3", "2h-1-4")]
lib_2h_rep2 <- aligned_sequence_select[,c("2h-2-1", "2h-2-2", "2h-2-3", "2h-2-4")]
lib_2h_rep3 <- aligned_sequence_select[,c("2h-3-1", "2h-3-2", "2h-3-3", "2h-3-4")]
colnames(lib_2h_rep1) <- c("2h_frac1", "2h_frac2", "2h_frac3", "2h_frac4")
colnames(lib_2h_rep2) <- c("2h_frac1", "2h_frac2", "2h_frac3", "2h_frac4")
colnames(lib_2h_rep3) <- c("2h_frac1", "2h_frac2", "2h_frac3", "2h_frac4")
lib_2h_rep1_MRL <- calculate_polysome_load_old(lib_2h_rep1, "2h")
lib_2h_rep2_MRL <- calculate_polysome_load_old(lib_2h_rep2, "2h")
lib_2h_rep3_MRL <- calculate_polysome_load_old(lib_2h_rep3, "2h")
lib_2h_rep1_mRNA <- rowSums(lib_2h_rep1)
lib_2h_rep2_mRNA <- rowSums(lib_2h_rep2)
lib_2h_rep3_mRNA <- rowSums(lib_2h_rep3)

lib_12h_rep1 <- aligned_sequence_select[,c("12h-2-1", "12h-2-2", "12h-2-3", "12h-2-4")]
lib_12h_rep2 <- aligned_sequence_select[,c("12h-3-1", "12h-3-2", "12h-3-3", "12h-3-4")]
lib_12h_rep3 <- aligned_sequence_select[,c("12h-4-1", "12h-4-2", "12h-4-3", "12h-4-4")]
colnames(lib_12h_rep1) <- c("12h_frac1", "12h_frac2", "12h_frac3", "12h_frac4")
colnames(lib_12h_rep2) <- c("12h_frac1", "12h_frac2", "12h_frac3", "12h_frac4")
colnames(lib_12h_rep3) <- c("12h_frac1", "12h_frac2", "12h_frac3", "12h_frac4")
lib_12h_rep1_MRL <- calculate_polysome_load_old(lib_12h_rep1, "12h")
lib_12h_rep2_MRL <- calculate_polysome_load_old(lib_12h_rep2, "12h")
lib_12h_rep3_MRL <- calculate_polysome_load_old(lib_12h_rep3, "12h")
lib_12h_rep1_mRNA <- rowSums(lib_12h_rep1[,3:4])
lib_12h_rep2_mRNA <- rowSums(lib_12h_rep2[,3:4])
lib_12h_rep3_mRNA <- rowSums(lib_12h_rep3[,3:4])

lib_1h_MRL <- (lib_1h_rep1_MRL$`1h` + lib_1h_rep2_MRL$`1h` + lib_1h_rep3_MRL$`1h`)/3
lib_2h_MRL <- (lib_2h_rep1_MRL$`2h` + lib_2h_rep2_MRL$`2h` + lib_2h_rep3_MRL$`2h`)/3
lib_12h_MRL <- (lib_12h_rep1_MRL$`12h` + lib_12h_rep3_MRL$`12h`)/2
lib_1h_mRNA <- (lib_1h_rep1_mRNA + lib_1h_rep2_mRNA + lib_1h_rep3_mRNA)/3
lib_2h_mRNA <- (lib_2h_rep1_mRNA + lib_2h_rep2_mRNA + lib_2h_rep3_mRNA)/3
lib_12h_mRNA <- (lib_12h_rep1_mRNA + lib_12h_rep3_mRNA)/2
mRNA_deg <- lib_12h_mRNA/lib_2h_mRNA

lib_MRLs <- data.frame(lib_1h_MRL, lib_2h_MRL, lib_12h_MRL)
time_points <- c(1, 2, 12)  # 1h, 2h, 12h
calculate_auc <- function(row) {
  trapz(time_points, row)
}


lib_MRLs_AUC <- apply(lib_MRLs, 1, calculate_auc)

MRL_diff_2h <- lib_2h_MRL - lib_1h_MRL
MRL_diff_12h <- lib_12h_MRL - lib_1h_MRL

coupling_data <- data.frame(sequence = aligned_sequence[aligned_sequence$file_tag %in% c("test", "mRNA_1273", "BNT-162b2"),"sequence"], 
                            lib_1h_MRL = lib_1h_MRL, lib_2h_MRL = lib_2h_MRL, lib_12h_MRL = lib_12h_MRL, AUC = lib_MRLs_AUC, MRL_diff = MRL_diff_2h,
                            lib_1h_mRNA = lib_1h_mRNA, lib_2h_mRNA = lib_2h_mRNA, lib_12h_mRNA = lib_12h_mRNA, mRNA_deg = mRNA_deg)


#write.table(coupling_data, file = "coupling_data_20240613.csv", quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
coupling_data <- merge(coupling_data, egfp_unmod_1_RL, by = "sequence")

sequence_selected <- read_excel("data/sequence_selected.xlsx")
mean_intensity <- read_csv("data/mean_intensity2.csv")

mean_intensity_summary <- mean_intensity %>%
  group_by(Condition = sub("-\\d$", "", Condition)) %>%
  summarise(mean_value = mean(mean_intensity))

sequence_selected_clean <- sequence_selected[,c(1,2,25)]
colnames(sequence_selected_clean) <- c("Condition", "sequence", "group")

sequence_selected_clean <- merge(sequence_selected_clean, mean_intensity_summary, by = "Condition")
sequence_selected_clean <- merge(sequence_selected_clean, coupling_data, by = "sequence")


coupling_data_with_GFP <- merge(coupling_data, sequence_selected_clean[, 1:4], by = "sequence", all.x = TRUE)
write.table(coupling_data_with_GFP, file = "data/IVT_supp.csv", quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
