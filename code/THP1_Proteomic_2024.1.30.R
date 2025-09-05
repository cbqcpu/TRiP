THP1_proteomic <- read_excel("data/THP1_Proteomic/41467_2021_26000_MOESM4_ESM.xlsx", skip = 2)
idmapping_2023_11_08 <- read_delim("data/THP1_Proteomic/idmapping_2023_11_08.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
colnames(idmapping_2023_11_08) <- c("Accession", "gene_id")
THP1_proteomic_id <- merge(THP1_proteomic, idmapping_2023_11_08, "Accession")
THP1_proteomic_id_protein_name <- data.frame(gene_id = THP1_proteomic_id$gene_id, Accession = THP1_proteomic_id$Accession, gene_name = THP1_proteomic_id$GN...54)
THP1_proteomic_0h <- data.frame(gene_id = THP1_proteomic_id$gene_id, rep1 = THP1_proteomic_id$`0 hr rep 1`, rep2 = THP1_proteomic_id$`0 hr rep 2`, rep3 = THP1_proteomic_id$`0 hr rep 3`)
THP1_proteomic_2h <- data.frame(gene_id = THP1_proteomic_id$gene_id, rep1 = THP1_proteomic_id$`2 hr rep 1`, rep2 = THP1_proteomic_id$`2 hr rep 2`, rep3 = THP1_proteomic_id$`2 hr rep 3`)
THP1_proteomic_4h <- data.frame(gene_id = THP1_proteomic_id$gene_id, rep1 = THP1_proteomic_id$`4 hr rep 1`, rep2 = THP1_proteomic_id$`4 hr rep 2`, rep3 = THP1_proteomic_id$`4 hr rep 3`)
THP1_proteomic_6h <- data.frame(gene_id = THP1_proteomic_id$gene_id, rep1 = THP1_proteomic_id$`6 hr rep 1`, rep2 = THP1_proteomic_id$`6 hr rep 2`, rep3 = THP1_proteomic_id$`6 hr rep 3`)
THP1_proteomic_12h <- data.frame(gene_id = THP1_proteomic_id$gene_id, rep1 = THP1_proteomic_id$`12 hr rep 1`, rep2 = THP1_proteomic_id$`12 hr rep 2`, rep3 = THP1_proteomic_id$`12 hr rep 3`)
THP1_proteomic_24h <- data.frame(gene_id = THP1_proteomic_id$gene_id, rep1 = THP1_proteomic_id$`24 hr rep 1`, rep2 = THP1_proteomic_id$`24 hr rep 2`, rep3 = THP1_proteomic_id$`24 hr rep 3`)

# Calculate mean for each time point and add it as a new column
THP1_proteomic_0h$mean_0h <- rowMeans(THP1_proteomic_0h[, c("rep1", "rep2", "rep3")])
THP1_proteomic_2h$mean_2h <- rowMeans(THP1_proteomic_2h[, c("rep1", "rep2", "rep3")])
THP1_proteomic_4h$mean_4h <- rowMeans(THP1_proteomic_4h[, c("rep1", "rep2", "rep3")])
THP1_proteomic_6h$mean_6h <- rowMeans(THP1_proteomic_6h[, c("rep1", "rep2", "rep3")])
THP1_proteomic_12h$mean_12h <- rowMeans(THP1_proteomic_12h[, c("rep1", "rep2", "rep3")])
THP1_proteomic_24h$mean_24h <- rowMeans(THP1_proteomic_24h[, c("rep1", "rep2", "rep3")])

THP1_proteomic_time <- data.frame(gene_id = THP1_proteomic_0h$gene_id, gene_name = THP1_proteomic_id_protein_name$gene_name, LPS_0h = THP1_proteomic_0h$mean_0h,
                                  LPS_2h = THP1_proteomic_2h$mean_2h,
                                  LPS_4h = THP1_proteomic_4h$mean_4h,
                                  LPS_6h = THP1_proteomic_6h$mean_6h,
                                  LPS_12h = THP1_proteomic_12h$mean_12h,
                                  LPS_24h = THP1_proteomic_24h$mean_24h)

THP1_proteomic_time_high_coupling <- THP1_proteomic_time[THP1_proteomic_time$gene_id %in% high_coupling_genes,]
THP1_proteomic_time_low_coupling <- THP1_proteomic_time[THP1_proteomic_time$gene_id %in% low_coupling_genes,]

library(ComplexHeatmap)

# Assuming 'df' is your dataframe
# Remove the 'gene_id' column
THP1_proteomic_time_high_coupling_for_heatmap <- THP1_proteomic_time_high_coupling[, c(-1,-2)]  # removes the first column

# Scale the data by row
THP1_proteomic_time_high_coupling_for_heatmap_scaled <- t(scale(t(THP1_proteomic_time_high_coupling_for_heatmap)))

# Draw the heatmap
Heatmap(THP1_proteomic_time_high_coupling_for_heatmap_scaled,
        cluster_rows = FALSE,
        cluster_columns = FALSE)


THP1_proteomic_time_low_coupling_for_heatmap <- THP1_proteomic_time_low_coupling[, c(-1,-2)]  # removes the first column

# Scale the data by row
THP1_proteomic_time_low_coupling_for_heatmap_scaled <- t(scale(t(THP1_proteomic_time_low_coupling_for_heatmap)))

# Draw the heatmap
Heatmap(THP1_proteomic_time_low_coupling_for_heatmap_scaled,
        cluster_rows = FALSE,
        cluster_columns = FALSE)

library(ComplexHeatmap)

# Assuming 'THP1_proteomic_time_high_coupling' and 'THP1_proteomic_time_low_coupling' are your dataframes

# Remove the 'gene_id' column
high_data <- THP1_proteomic_time_high_coupling[, c(-1,-2)]
low_data <- THP1_proteomic_time_low_coupling[, c(-1,-2)]

high_data_long <- high_data %>%
  pivot_longer(cols = everything(), names_to = "Time_Point", values_to = "Mean") %>%
  mutate(Time_Point = factor(Time_Point, levels = names(high_data)))

# Calculate the mean for each time point
high_data_summary <- high_data_long %>%
  group_by(Time_Point) %>%
  summarize(Mean = mean(Mean))

# Draw the line plot with ggplot
ggplot(high_data_summary, aes(x = Time_Point, y = Mean, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Time Point", y = "Mean Expression", title = "Mean Expression Over Time")


low_data_long <- low_data %>%
  pivot_longer(cols = everything(), names_to = "Time_Point", values_to = "Mean") %>%
  mutate(Time_Point = factor(Time_Point, levels = names(low_data)))

# Calculate the mean for each time point
low_data_summary <- low_data_long %>%
  group_by(Time_Point) %>%
  summarize(Mean = mean(Mean))

# Draw the line plot with ggplot
ggplot(low_data_summary, aes(x = Time_Point, y = Mean, group = 1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Time Point", y = "Mean Expression", title = "Mean Expression Over Time")




