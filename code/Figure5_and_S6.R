library(pracma)
library(gridExtra)
library(ggplot2)
library(readxl)
library(scales)
library(dplyr)
library(tidyr)
library(ggpubr)
library(FNN)
library(ComplexHeatmap)
library(circlize)

source("code/polysome_seq_functions.R")
egfp_unmod_1 <- read_csv("data/GSM3130435_egfp_unmod_1.csv")
aligned_sequence <- read_csv("data/aligned_sequence.csv")
test_set_2000_MFE <- read_csv("data/coupling_data_20240614_with_gfp_MFE_NUPACK.csv")

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
time_points <- c(1, 2, 12)
calculate_auc <- function(row) {
  trapz(time_points, row)
}
lib_MRLs_AUC <- apply(lib_MRLs, 1, calculate_auc)
MRL_diff_2h <- lib_2h_MRL - lib_1h_MRL
coupling_data <- data.frame(sequence = aligned_sequence[aligned_sequence$file_tag %in% c("test", "mRNA_1273", "BNT-162b2"),"sequence"],
                            lib_1h_MRL = lib_1h_MRL, lib_2h_MRL = lib_2h_MRL, lib_12h_MRL = lib_12h_MRL, AUC = lib_MRLs_AUC, MRL_diff = MRL_diff_2h,
                            lib_1h_mRNA = lib_1h_mRNA, lib_2h_mRNA = lib_2h_mRNA, lib_12h_mRNA = lib_12h_mRNA, mRNA_deg = mRNA_deg)
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
sequence_selected_clean_MRL_high <- sequence_selected_clean[sequence_selected_clean$group=="High",]
sequence_selected_clean_MRL_high <- sequence_selected_clean_MRL_high[sequence_selected_clean_MRL_high$Condition!="UTR-23",]
coupling_data_with_GFP <- merge(coupling_data, sequence_selected_clean[, 1:4], by = "sequence", all.x = TRUE)
coupling_data_with_GFP_select <- coupling_data_with_GFP[coupling_data_with_GFP$mRNA_deg<0.5,]

# Fig5I
MRL_2h_no_difference_high <- coupling_data_with_GFP[(coupling_data_with_GFP$lib_2h_MRL>1.5)&(coupling_data_with_GFP$lib_2h_MRL<2.6)&(coupling_data_with_GFP$group=="High")&!is.na(coupling_data_with_GFP$mean_value),]
MRL_2h_no_difference_low <- coupling_data_with_GFP[(coupling_data_with_GFP$lib_2h_MRL>2)&(coupling_data_with_GFP$lib_2h_MRL<3)&(coupling_data_with_GFP$group=="Low")&!is.na(coupling_data_with_GFP$mean_value),]
MRL_2h_no_difference <- rbind(MRL_2h_no_difference_high, MRL_2h_no_difference_low)
MRL_2h_no_difference <- MRL_2h_no_difference[!(MRL_2h_no_difference$Condition %in% c("UTR-23", "UTR-76")), ]
t_test_MRL_2h_no_difference_p1 <- t.test(
  MRL_2h_no_difference[MRL_2h_no_difference$group=="High",]$lib_2h_MRL,
  MRL_2h_no_difference[MRL_2h_no_difference$group=="Low",]$lib_2h_MRL
)
p_value_MRL_2h_no_difference_p1 <- t_test_MRL_2h_no_difference_p1$p.value
Fig5I_1_MRL_2h_no_difference_p1 <- ggplot(MRL_2h_no_difference, aes(x = group, y = lib_2h_MRL, alpha = 0.5)) +
  geom_boxplot(aes(fill = group)) +
  scale_fill_manual(values = c("High" = "#bd0026", "Low" = "#253494")) +
  annotate("text", x = 1.5, y = max(MRL_2h_no_difference$lib_2h_MRL),
           label = paste("p-value:", p_value_MRL_2h_no_difference_p1), size = 5) +
  labs(title = "Comparison of Mean Values",
       x = "Group",
       y = "2h MRL") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(Fig5I_1_MRL_2h_no_difference_p1)
t_test_MRL_2h_no_difference_p2 <- t.test(
  MRL_2h_no_difference[MRL_2h_no_difference$group=="High",]$mean_value,
  MRL_2h_no_difference[MRL_2h_no_difference$group=="Low",]$mean_value
)
p_value_MRL_2h_no_difference_p2 <- t_test_MRL_2h_no_difference_p2$p.value
Fig5I_3_MRL_2h_no_difference_p2 <- ggplot(MRL_2h_no_difference, aes(x = group, y = mean_value, alpha = 0.5)) +
  geom_boxplot(aes(fill = group)) +
  scale_fill_manual(values = c("High" = "#bd0026", "Low" = "#253494")) +
  annotate("text", x = 1.5, y = max(MRL_2h_no_difference$mean_value),
           label = paste("p-value:", formatC(p_value_MRL_2h_no_difference_p2, format = "e", digits = 2)), size = 5) +
  labs(title = "Comparison of Mean Values",
       x = "Group",
       y = "GFP Intensity A.U.") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(Fig5I_3_MRL_2h_no_difference_p2)
t_test_MRL_2h_no_difference_p3 <- t.test(
  MRL_2h_no_difference[MRL_2h_no_difference$group=="High",]$MRL_diff,
  MRL_2h_no_difference[MRL_2h_no_difference$group=="Low",]$MRL_diff
)
p_value_MRL_2h_no_difference_p3 <- t_test_MRL_2h_no_difference_p3$p.value
Fig5I_2_MRL_2h_no_difference_p3 <- ggplot(MRL_2h_no_difference, aes(x = group, y = MRL_diff, alpha = 0.5)) +
  geom_boxplot(aes(fill = group)) +
  scale_fill_manual(values = c("High" = "#bd0026", "Low" = "#253494")) +
  annotate("text", x = 1.5, y = max(MRL_2h_no_difference$MRL_diff),
           label = paste("p-value:", formatC(p_value_MRL_2h_no_difference_p3, format = "e", digits = 2)), size = 5) +
  labs(title = "Comparison of Mean Values",
       x = "Group",
       y = "MRL diff.") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(Fig5I_2_MRL_2h_no_difference_p3)

# Fig5F
cor_results <- cor.test(test_set_2000_MFE$MRL_diff, test_set_2000_MFE$MFE, method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)
k <- 10
neighbors <- get.knn(data = test_set_2000_MFE[, c("MRL_diff", "MFE")], k = k)
test_set_2000_MFE$density <- 1 / rowMeans(neighbors$nn.dist)
test_set_2000_MFE$density <- test_set_2000_MFE$density / max(test_set_2000_MFE$density)
custom_colors <- c("#54a5de", "white", "#f15389")
Fig5F_compare_MFE_regression <- ggplot(test_set_2000_MFE, aes(x = MRL_diff, y = MFE, color = density)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.4, 1), name = "Density") +
  labs(
    title = "Scatter Plot of MFE vs MRL_diff",
    x = "MRL_diff",
    y = "Minimum Free Energy (MFE)"
  ) +
  annotate("text", x = 1, y = 1.4, label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(Fig5F_compare_MFE_regression)

# Fig5G
cor_results <- cor.test(coupling_data_with_GFP_select$MRL_diff, coupling_data_with_GFP_select$mRNA_deg, method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)
k <- 10
neighbors <- get.knn(data = coupling_data_with_GFP_select[, c("MRL_diff", "mRNA_deg")], k = k)
coupling_data_with_GFP_select$density <- 1 / rowMeans(neighbors$nn.dist)
coupling_data_with_GFP_select$density <- coupling_data_with_GFP_select$density / max(coupling_data_with_GFP_select$density)
custom_colors <- c("#54a5de", "white", "#f15389")
Fig5G_plot1 <- ggplot(coupling_data_with_GFP_select, aes(x = MRL_diff, y = mRNA_deg, color = density)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.4, 1), name = "Density") +
  labs(
    x = "MRL Difference",
    y = "mRNA (12 h / 2 h)",
    title = "Scatter Plot of MRL-diff vs. mRNA Degradation"
  ) +
  annotate("text", x = 0.4, y = 0.4, label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(Fig5G_plot1)

# Fig5D
pattern_data <- coupling_data %>%
  mutate(
    pattern = case_when(
      lib_2h_MRL > lib_1h_MRL & lib_2h_MRL > lib_12h_MRL ~ "Rise and Fall",
      lib_12h_MRL > lib_2h_MRL & lib_2h_MRL > lib_1h_MRL ~ "Monotonic Increase",
      TRUE ~ "Other"
    )
  )
rise_and_fall_count <- sum(pattern_data$pattern == "Rise and Fall")
monotonic_increase_count <- sum(pattern_data$pattern == "Monotonic Increase")
other_count <- sum(pattern_data$pattern == "Other")
long_pattern_data <- pattern_data %>%
  pivot_longer(
    cols = starts_with("lib_"),
    names_to = "time_point",
    values_to = "MRL"
  ) %>%
  mutate(
    time_point = recode(time_point,
                        "lib_1h_MRL" = "1h",
                        "lib_2h_MRL" = "2h",
                        "lib_12h_MRL" = "12h"),
    time_point = factor(time_point, levels = c("1h", "2h", "12h"))
  ) %>%
  filter(!is.na(time_point) & !is.na(MRL))
summary_data <- long_pattern_data %>%
  group_by(pattern, time_point) %>%
  summarise(
    mean_MRL = mean(MRL),
    sd_MRL = sd(MRL),
    .groups = 'drop'
  ) %>%
  mutate(
    pattern_label = case_when(
      pattern == "Rise and Fall" ~ paste0("Rise and Fall (n=", rise_and_fall_count, ")"),
      pattern == "Monotonic Increase" ~ paste0("Monotonic Increase (n=", monotonic_increase_count, ")"),
      pattern == "Other" ~ paste0("Other (n=", other_count, ")")
    )
  )
Fig5D_p_two_patterns <- ggplot(summary_data, aes(x = time_point, y = mean_MRL, group = pattern_label, color = pattern_label)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_MRL - sd_MRL, ymax = mean_MRL + sd_MRL), width = 0.2, size = 0.8) +
  scale_color_manual(values = c("steelblue", "darkorange", "gray50")) +
  labs(
    x = "Time Point",
    y = "Mean MRL ± SD",
    color = "Pattern",
    title = "Mean MRL ± SD Over Time for Rise and Fall, Monotonic Increase, and Other"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(Fig5D_p_two_patterns)

# Fig5E
pattern_data <- coupling_data %>%
  mutate(
    pattern = case_when(
      lib_2h_MRL > lib_1h_MRL & lib_2h_MRL > lib_12h_MRL ~ "Rise and Fall",
      lib_12h_MRL > lib_2h_MRL & lib_2h_MRL > lib_1h_MRL ~ "Monotonic Increase",
      TRUE ~ "Other"
    ),
    pattern = factor(pattern, levels = c("Rise and Fall", "Monotonic Increase", "Other"))
  )
pattern_data <- pattern_data %>%
  mutate(MRL_diff_2h_1h = lib_2h_MRL - lib_1h_MRL)
comparison_data <- pattern_data %>%
  filter(pattern %in% c("Rise and Fall", "Monotonic Increase"))
Fig5E <- ggplot(comparison_data, aes(x = pattern, y = MRL_diff_2h_1h, fill = pattern)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) +
  scale_fill_manual(values = c("Rise and Fall" = "#95c08b", "Monotonic Increase" = "#8887b9")) +
  labs(
    title = "Comparison of MRL Differences (2h - 1h)",
    x = "Pattern",
    y = expression(Delta ~ "MRL (2h - 1h)")
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(Fig5E)

# Fig5C
pattern_data <- coupling_data %>%
  mutate(
    pattern = case_when(
      lib_2h_MRL > lib_1h_MRL & lib_2h_MRL > lib_12h_MRL ~ "Rise and Fall",
      lib_12h_MRL > lib_2h_MRL & lib_2h_MRL > lib_1h_MRL ~ "Monotonic Increase",
      TRUE ~ "Other"
    ),
    pattern = factor(pattern, levels = c("Rise and Fall", "Monotonic Increase", "Other"))
  )
pattern_data <- pattern_data %>%
  arrange(pattern)
heatmap_data <- pattern_data %>%
  dplyr::select(lib_1h_MRL, lib_2h_MRL, lib_12h_MRL) %>%
  as.matrix()
heatmap_data_scaled <- t(apply(heatmap_data, 1, function(x) (x - mean(x)) / sd(x)))
pattern_colors <- c("Rise and Fall" = "#95c08b", "Monotonic Increase" = "#8887b9", "Other" = "#b0b0b0")
pattern_annotations <- rowAnnotation(
  Pattern = pattern_data$pattern,
  col = list(Pattern = pattern_colors),
  annotation_legend_param = list(Pattern = list(title = "Pattern"))
)
Fig5C <- Heatmap(heatmap_data_scaled,
                 name = "MRL (Z-score)",
                 left_annotation = pattern_annotations,
                 col = colorRamp2(c(min(heatmap_data_scaled), mean(heatmap_data_scaled), max(heatmap_data_scaled)),
                                  c("white", "orange", "darkred")),
                 show_row_names = FALSE,
                 cluster_columns = FALSE,
                 cluster_rows = FALSE,
                 column_title = "Time Points (1h, 2h, 12h)",
                 row_title = "Sequences"
)
print(Fig5C)

# FigS6D
top_10_percent <- read_csv("data/test_set_top_10_percent_MFE_NUPACK.csv")
bottom_10_percent <- read_csv("data/test_set_bottom_10_percent_MFE_NUPACK.csv")
top_10_percent <- top_10_percent %>%
  mutate(Group = "Top 10%")
bottom_10_percent <- bottom_10_percent %>%
  mutate(Group = "Bottom 10%")
combined_data <- bind_rows(top_10_percent, bottom_10_percent)
FigS6D_compare_MFE <- ggboxplot(combined_data, x = "Group", y = "MFE", fill = "Group",
                                palette = c("Top 10%" = "#bd0026", "Bottom 10%" = "#253494"),
                                title = "Comparison of MFE Between Top 10% and Bottom 10%",
                                xlab = "MRL-diff", ylab = "Minimum Free Energy (MFE)",
                                alpha = 0.5) +
  stat_compare_means(method = "t.test", label = "p.format") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(FigS6D_compare_MFE)

# FigS6A
model <- lm(lib_12h_MRL ~ RL, data = coupling_data_with_GFP)
r_squared <- summary(model)$r.squared
cor_test_FigS6A <- cor.test(coupling_data_with_GFP$RL, coupling_data_with_GFP$lib_12h_MRL)
FigS6A_p3 <- ggplot() +
  geom_point(data = subset(coupling_data_with_GFP, is.na(group)),
             aes(x = RL, y = lib_12h_MRL, color = "NA"), alpha = 0.1) +
  geom_point(data = subset(coupling_data_with_GFP, !is.na(group)),
             aes(x = RL, y = lib_12h_MRL, color = group), alpha = 0.5) +
  scale_color_manual(
    values = c("High" = "#bd0026", "Low" = "#253494", "NA" = "#4ba2dd"),
    na.value = "#4ba2dd"
  ) +
  geom_smooth(data = coupling_data_with_GFP, aes(x = RL, y = lib_12h_MRL), method = "lm", se = FALSE, color = "#08306b") +
  annotate("text", x = Inf, y = Inf, label = paste("R^2 =", round(r_squared, 4)), hjust = 1.1, vjust = 1.1, size = 5) +
  labs(title = "RL vs lib_12h_MRL", x = "RL", y = "lib_12h_MRL") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )
print(FigS6A_p3)





library(dplyr)
library(stringr)
library(ggplot2)

# Load data
RandomDesign_50_seq_1m <- read_csv("data/RandomDesign_50_seq_1m.csv")

# Rescale MRL_2h and MRL_diff
RandomDesign_50_seq_1m <- RandomDesign_50_seq_1m %>%
  mutate(
    MRL_2h = MRL_2h / 10,
    MRL_diff = MRL_diff / 10
  )

# Define function for reverse complement
reverse_complement <- function(seq) {
  complement <- chartr("ACGT", "TGCA", seq)
  paste(rev(strsplit(complement, NULL)[[1]]), collapse = "")
}

# Define restriction sites and their reverse complements
motifs <- c("CGTCTC", "TCTAGA")
reverse_motifs <- sapply(motifs, reverse_complement)
all_motifs <- c(motifs, reverse_motifs)

# Filter sequences without motifs and within MRL_2h range
filtered_data <- RandomDesign_50_seq_1m %>%
  mutate(contains_motif = str_detect(sequence, paste(all_motifs, collapse = "|"))) %>%
  filter(!contains_motif, MRL_2h >= 0.9, MRL_2h <= 1.1)

# Select top and bottom 10 sequences by MRL_diff
top_10 <- filtered_data %>%
  arrange(desc(MRL_diff)) %>%
  head(10) %>%
  mutate(group = "top_10")
bottom_10 <- filtered_data %>%
  arrange(MRL_diff) %>%
  head(10) %>%
  mutate(group = "bottom_10")

# Load and summarize mean intensity data
mean_intensity <- read_csv("data/mean_intensity_IVT_selected_expression.csv")
mean_intensity_summary <- mean_intensity %>%
  group_by(Condition = sub("-\\d$", "", Condition)) %>%
  summarise(mean_value = mean(mean_intensity)) %>%
  filter(Condition != "NC", !str_detect(Condition, "low_MRL_diff_10"))

# Combine and merge data
combined_data <- bind_rows(top_10, bottom_10) %>%
  group_by(group) %>%
  mutate(Condition = paste0(group, "-", row_number())) %>%
  ungroup()
combined_data_gfp <- merge(combined_data, mean_intensity_summary, by = "Condition")

# Filter specific conditions for plotting
top_10_data <- combined_data_gfp %>%
  filter(group == "top_10", !Condition %in% c("top_10-7", "top_10-8")) %>%
  arrange(MRL_2h)
bottom_10_data <- combined_data_gfp %>%
  filter(group == "bottom_10", !Condition %in% c("bottom_10-8", "bottom_10-10")) %>%
  arrange(MRL_2h)

# Create Figure 5K: Density plot with highlighted points
Fig5K_predicted <- ggplot(RandomDesign_50_seq_1m, aes(x = MRL_diff, y = MRL_2h)) +
  geom_density_2d() +
  geom_point(data = top_10_data, aes(x = MRL_diff, y = MRL_2h), color = "#ae017e", size = 3, alpha = 0.5) +
  geom_point(data = bottom_10_data, aes(x = MRL_diff, y = MRL_2h), color = "#225ea8", size = 3, alpha = 0.5) +
  labs(
    x = "MRL Difference",
    y = "MRL 2 Hours",
    title = "Density Plot of MRL_diff vs MRL_2h"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

# Output the plot
print(Fig5K_predicted)


library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Reading and summarizing the data
mean_intensity <- read_csv("data/mean_intensity_flow_cytometry.csv")

mean_intensity_summary <- mean_intensity %>%
  group_by(Condition = sub("-\\d$", "", Condition)) %>%
  summarise(mean_value = mean(mean_intensity))


# Assuming 'mean_intensity_summary' is your provided data frame
# Create a new grouping variable based on the Condition
mean_intensity_summary <- mean_intensity_summary %>%
  mutate(Group = case_when(
    Condition %in% paste0("", 31:50) ~ "Group 1 High MRL High coupling",
    Condition %in% paste0("", c(3,4,10)) ~ "Group 1 High MRL High coupling",
    Condition %in% paste0("", 51:70) ~ "Group 2 High MRL Low coupling",
    Condition %in% paste0("", c(12,14,16)) ~ "Group 2 High MRL Low coupling",
    TRUE ~ "Other"
  )) %>%
  filter(Group != "Other")

mean_intensity_summary <- mean_intensity_summary %>%
  filter(Condition!="38")


p_value <- t.test(mean_value ~ Group, data = mean_intensity_summary %>%
                    filter(Group %in% c("Group 1 High MRL High coupling", "Group 2 High MRL Low coupling")))$p.value

median_values <- mean_intensity_summary %>%
  group_by(Group) %>%
  summarize(median_value = median(mean_value, na.rm = TRUE))


# Create the plot with direct p-value annotation
Fig5L_p <- ggplot(mean_intensity_summary, aes(x = Group, y = mean_value, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Comparison of GFP Expression by MRL-diff",
    x = "Group",
    y = "GFP Intensity (A.U.)"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  # Add p-value annotation
  annotate("text", x = 1.5, y = max(mean_intensity_summary$mean_value) * 1.1,
           label = paste0("p = ", format.pval(p_value, digits = 3)),
           size = 5, color = "black") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )




LNP_mouse_liver_GFP_quant <- read_csv("data/LNP_mouse_liver_GFP_quant.csv")

# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(viridis)
library(ggpubr) # For adding p-value annotations

# Step 1: Separate Group and Replicate
LNP_mouse_liver_GFP_quant <- LNP_mouse_liver_GFP_quant %>%
  mutate(
    Group = str_extract(Sample, "^[^-]+"), # Extract main group (before the '-')
    Replicate = str_extract(Sample, "[^-]+$") # Extract replicate (after the '-')
  )

# Ensure Group order is BNT, SEQ59, SEQ67, PBS
LNP_mouse_liver_GFP_quant <- LNP_mouse_liver_GFP_quant %>%
  mutate(Group = factor(Group, levels = c("BNT", "SEQ59", "SEQ67", "PBS")))

# Step 2: Perform t-tests
# Filter and prepare data for SEQ59 vs BNT
data_SEQ59_BNT <- LNP_mouse_liver_GFP_quant %>%
  filter(Group %in% c("SEQ59", "BNT"))

t_test_SEQ59_BNT <- t.test(`GFP norm GAPDH` ~ Group, data = data_SEQ59_BNT)

# Filter and prepare data for SEQ67 vs BNT
data_SEQ67_BNT <- LNP_mouse_liver_GFP_quant %>%
  filter(Group %in% c("SEQ67", "BNT"))

t_test_SEQ67_BNT <- t.test(`GFP norm GAPDH` ~ Group, data = data_SEQ67_BNT)

# Extract p-values
p_val_SEQ59_BNT <- t_test_SEQ59_BNT$p.value
p_val_SEQ67_BNT <- t_test_SEQ67_BNT$p.value

# Step 3: Create Boxplot with p-value Annotations
Fig5M <- GFP_expression_mouse_liver <- ggplot(LNP_mouse_liver_GFP_quant, aes(x = Group, y = `GFP norm GAPDH`, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  scale_fill_viridis_d() + # Automatically handles sufficient colors
  theme_bw() +  # Use a minimal theme
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  ) +
  labs(
    title = "Comparison of GFP norm GAPDH across Groups",
    x = "Group",
    y = "GFP norm GAPDH"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(
    comparisons = list(c("SEQ59", "BNT"), c("SEQ67", "BNT")),
    method = "t.test",
    label = "p.format" # Format the p-values nicely
  )



train_MRL <- read_csv("data/MRL_diff_output_model_20241028/train_MRL.csv")
train_MRL_diff <- read_csv("data/MRL_diff_output_model_20241028/train_MRL_diff.csv")

valid_MRL <- read_csv("data/MRL_diff_output_model_20241028/valid_MRL.csv")
valid_MRL_diff <- read_csv("data/MRL_diff_output_model_20241028/valid_MRL_diff.csv")

library(ggplot2)

# Function to calculate correlation and format it
add_corr <- function(data, x, y) {
  corr <- cor(data[[x]], data[[y]], use = "complete.obs")
  paste0("Corr: ", round(corr, 2))
}

# Plot 1: train_MRL$Norm_MRL vs. train_MRL$MRL_2h_pred
cor_test_FigS6E_train_MRL <- cor.test(train_MRL$Norm_MRL,train_MRL$MRL_2h_pred)
FigS6E_train_MRL <- ggplot(train_MRL, aes(x = Norm_MRL, y = MRL_2h_pred)) +
  geom_point(color = "#4ba2dd", alpha = 0.3) +
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +
  annotate(
    "text", x = Inf, y = Inf, label = add_corr(train_MRL, "Norm_MRL", "MRL_2h_pred"), 
    hjust = 1.1, vjust = 2, size = 5, color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "Train MRL",
    x = "Normalized MRL",
    y = "MRL 2h Prediction"
  )

# Plot 2: train_MRL_diff$Norm_MRL_diff vs. train_MRL_diff$MRL_diff_pred
cor_test_FigS6E_train_MRL_diff <- cor.test(train_MRL_diff$Norm_MRL_diff,train_MRL_diff$MRL_diff_pred)
FigS6E_train_MRL_diff <- ggplot(train_MRL_diff, aes(x = Norm_MRL_diff, y = MRL_diff_pred)) +
  geom_point(color = "#4ba2dd", alpha = 0.3) +
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +
  annotate(
    "text", x = Inf, y = Inf, label = add_corr(train_MRL_diff, "Norm_MRL_diff", "MRL_diff_pred"), 
    hjust = 1.1, vjust = 2, size = 5, color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "Train MRL-diff",
    x = "Normalized MRL Difference",
    y = "MRL Difference Prediction"
  )

# Plot 3: valid_MRL$Norm_MRL vs. valid_MRL$MRL_2h_pred
cor_test_FigS6E_valid_MRL <- cor.test(valid_MRL$Norm_MRL,valid_MRL$MRL_2h_pred)
FigS6E_valid_MRL <- ggplot(valid_MRL, aes(x = Norm_MRL, y = MRL_2h_pred)) +
  geom_point(color = "#4ba2dd", alpha = 0.3) +
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +
  annotate(
    "text", x = Inf, y = Inf, label = add_corr(valid_MRL, "Norm_MRL", "MRL_2h_pred"), 
    hjust = 1.1, vjust = 2, size = 5, color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "Validation MRL 2h",
    x = "Normalized MRL",
    y = "MRL 2h Prediction"
  )

# Plot 4: valid_MRL_diff$Norm_MRL_diff vs. valid_MRL_diff$MRL_diff_pred
cor_test_FigS6E_valid_MRL_diff <- cor.test(valid_MRL_diff$Norm_MRL_diff,valid_MRL_diff$MRL_diff_pred)
FigS6E_valid_MRL_diff <- ggplot(valid_MRL_diff, aes(x = Norm_MRL_diff, y = MRL_diff_pred)) +
  geom_point(color = "#4ba2dd", alpha = 0.3) +
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +
  annotate(
    "text", x = Inf, y = Inf, label = add_corr(valid_MRL_diff, "Norm_MRL_diff", "MRL_diff_pred"), 
    hjust = 1.1, vjust = 2, size = 5, color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "Validation MRL-diff",
    x = "Normalized MRL Difference",
    y = "MRL Difference Prediction"
  )



independent_test_withGFP_MRL <- read_csv("data/MRL_diff_output_model_20241028/independent_test_withGFP_MRL.csv")
independent_test_withGFP_MRL_diff <- read_csv("data/MRL_diff_output_model_20241028/independent_test_withGFP_MRL_diff.csv")

# Define a helper function to compute correlation and p-value
add_corr_and_p <- function(data, x_col, y_col) {
  # Calculate correlation coefficient
  cor_test <- cor.test(data[[x_col]], data[[y_col]])
  r_value <- cor_test$estimate  # Correlation coefficient
  p_value <- cor_test$p.value   # P-value
  
  # Format the label with both r and p
  label <- paste0("r = ", format(r_value, digits = 3), "\n",
                  "p = ", format(p_value, digits = 3, scientific = TRUE))
  return(label)
}

# First plot with labels at top-left
cor_test_FigS6E_independent_test_withGFP_MRL_plot <- cor.test(independent_test_withGFP_MRL$Norm_MRL,independent_test_withGFP_MRL$MRL_2h_pred)
FigS6E_independent_test_withGFP_MRL_plot <- ggplot(independent_test_withGFP_MRL, aes(x = Norm_MRL, y = MRL_2h_pred)) +
  geom_point(color = "#4ba2dd", alpha = 0.3) +
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +
  annotate(
    "text", x = -Inf, y = Inf, 
    label = add_corr_and_p(independent_test_withGFP_MRL, "Norm_MRL", "MRL_2h_pred"), 
    hjust = -0.1, vjust = 1.1, size = 5, color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "independent_test MRL",
    x = "Normalized MRL",
    y = "MRL Prediction"
  )

# Second plot with labels at top-left
FigS6E_independent_test_withGFP_MRL_diff_plot <- ggplot(independent_test_withGFP_MRL_diff, aes(x = Norm_MRL_diff, y = MRL_diff_pred)) +
  geom_point(color = "#4ba2dd", alpha = 0.3) +
  geom_smooth(method = "lm", color = "#08306b", se = FALSE) +
  annotate(
    "text", x = -Inf, y = Inf, 
    label = add_corr_and_p(independent_test_withGFP_MRL_diff, "Norm_MRL_diff", "MRL_diff_pred"), 
    hjust = -0.1, vjust = 1.1, size = 5, color = "black"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "independent_test MRL-diff",
    x = "Normalized MRL Difference",
    y = "MRL Difference Prediction"
  )





