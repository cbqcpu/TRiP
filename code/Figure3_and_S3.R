source("R/polysome_seq_functions_20230703.R")
load("data/HEK293T_preprocessed.RData")
load("data/RNA_features_gene_level_20240418.RData")
library(readr)
library(pracma)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(effsize)
source("R/RNA_features_new.R")

merge_df_counts_select <- merge_df_counts[rowMeans(merge_df_counts[,-1])>1000,]
merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts$gene_id %in% merge_df_counts_select$gene_id,]
merge_df_counts_old_select <- data.frame(gene_id = merge_df_counts_select$gene_id, merge_df_counts_select[,-1] - merge_df_counts_new_select[,-1])

#MRL_total <- calculate_polysome_load(merge_df_counts_select, c("NC_30min", "NC_1h", "NC_2h"))
MRLs <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_30min", "NC_1h", "NC_2h"))
MRL_new <- MRLs[[1]]
MRL_old <- MRLs[[2]]
MRL_total <- (MRL_new + MRL_old)/2
MRL_diff <- MRL_new-MRL_old

MRL_diff_output <- MRL_diff
MRL_diff_output$gene_id <- merge_df_counts_select$gene_id
#write.table(MRL_diff_output, file = "data/gene_level_MRL_diff.csv", quote = FALSE, sep = ",",row.names = FALSE)


time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new <- apply(MRL_new, 1, calculate_auc)
auc_MRL_old <- apply(MRL_old, 1, calculate_auc)

auc_diff <- auc_MRL_new - auc_MRL_old
#auc_diff <- MRL_diff
auc_diff_ranked <- rank(auc_diff, ties.method = "first")
high_coupling_indices <- which(auc_diff_ranked > length(auc_diff_ranked) - 2000)
low_coupling_indices <- which(auc_diff_ranked <= 2000)


# Calculate the values for the top and bottom 2000 indices
high_coupling_value <- min(auc_diff[high_coupling_indices])
low_coupling_value <- max(auc_diff[low_coupling_indices])
high_coupling_genes <- merge_df_counts_select[high_coupling_indices,1]
low_coupling_genes <- merge_df_counts_select[low_coupling_indices,1]

#save(high_coupling_genes, low_coupling_genes, file = "HEK293_high_low_coupling_genes.RData")
#write.table(high_coupling_genes, file = "./high_coupling_genes.csv", quote = FALSE, sep = ',')
#write.table(low_coupling_genes, file = "./low_coupling_genes.csv", quote = FALSE, sep = ',')
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

auc_diff_gene_id <- data.frame(gene_id = merge_df_counts_select$gene_id, coupling = auc_diff)
MRL_diff_gene_id <- data.frame(gene_id = merge_df_counts_select$gene_id, coupling = MRL_diff)
RNA_feature_auc_diff <- merge(RNA_features_gene_level, auc_diff_gene_id, by = "gene_id")
RNA_feature_MRL_diff <- merge(RNA_features_gene_level, MRL_diff_gene_id, by = "gene_id")

# Create the density plot with adjusted vertical lines
#pdf(file = "results/Figure1/AUC_diff_293T.pdf", width = 3.5, height = 2.5)
ggplot() +
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
#dev.off()


############################ RNA features analsis ####################################

RNA_feature_high_coupling <- RNA_features_gene_level[RNA_features_gene_level$gene_id %in% high_coupling_genes,]
RNA_feature_low_coupling <- RNA_features_gene_level[RNA_features_gene_level$gene_id %in% low_coupling_genes,]

# Assuming the column names are stored in a variable called `colnames_vector`
colnames_vector <- colnames(RNA_features_gene_level)

# Filter out 'gene_id' and 'gene_name'
# Create a mapping from original column names to publication-friendly names
column_name_mapping <- c(
  "transcript_length" = "Transcript Length",
  "5utr_gc_content" = "5' UTR GC Content",
  "5utr_length" = "5' UTR Length",
  "5utrcap_MFE" = "5' UTR Cap MFE",
  "5utr_utr5_cap_gc_content" = "5' UTR Cap GC Content",
  "utr5_MFE_nupack" = "5' UTR MFE Nupack",
  "start_kozak_score" = "Start Kozak Score",
  "start_uORF_count" = "Start uORF Count",
  "cds_gc_content" = "CDS GC Content",
  "cds_length" = "CDS Length",
  "cds_MFE_min" = "CDS MFE Min",
  "cds_avg_codon_freq" = "CDS Avg Codon Frequency",
  "cds_min_codon_freq" = "CDS Min Codon Frequency",
  "3utr_gc_content" = "3' UTR GC Content",
  "3utr_length" = "3' UTR Length",
  "3utr_MFE_min" = "3' UTR MFE Min",
  "3utr_au_element_count" = "3' UTR AU Element Count",
  "3utr_au_element_frac" = "3' UTR AU Element Fraction",
  "3utr_max_au_length" = "3' UTR Max AU Length",
  "exons_gc_content" = "Exons GC Content",
  "exons_length" = "Exons Length",
  "exons_exonct" = "Exons Exon Count",
  "introns_gc_content" = "Introns GC Content",
  "introns_length" = "Introns Length",
  "m6A_5utr_average_level" = "m6A 5' UTR Avg Level",
  "m6A_cds_average_level" = "m6A CDS Avg Level",
  "m6A_3utr_average_level" = "m6A 3' UTR Avg Level",
  "m6A_exons_average_level" = "m6A Exons Avg Level",
  "m6A_introns_average_level" = "m6A Introns Avg Level",
  "RNA_half_life" = "RNA Half-Life"
)

# Filter out unwanted columns
filtered_colnames <- colnames_vector[!colnames_vector %in% c("gene_id", "gene_name",
                                                             "5utr_MFE_min", "start_uORF_overlap", "m6A_5utr_site_num", "m6A_cds_site_num", "m6A_3utr_site_num",
                                                             "m6A_exons_site_num", "m6A_introns_site_num")]

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each column name (excluding 'gene_id' and 'gene_name') to plot
for(i in seq_along(filtered_colnames)) {
  feature_column <- filtered_colnames[i]
  
  # Print the feature_column name to keep track of which plot is being generated
  print(paste("Generating plot for:", feature_column))
  
  # Use the friendly name if it exists, otherwise use the original name
  friendly_name <- column_name_mapping[feature_column]
  if (is.null(friendly_name)) {
    friendly_name <- feature_column
  }
  
  # Call your compare_boxplot function for each feature
  plot_list[[feature_column]] <- compare_boxplot(
    RNA_feature_high_coupling, 
    RNA_feature_low_coupling, 
    feature_column, 
    name1 = "High", 
    name2 = "Low"
  ) + labs(title = friendly_name, y = friendly_name)
  
  # If you're running this interactively, you might want to uncomment the next line to pause between plots.
  # readline(prompt="Press [enter] to continue")
}


# Split the list into two parts
FigS3A_plots_part1 <- plot_list

calculate_cliffs_delta_for_all_features <- function(high_coupling_df, low_coupling_df) {
  feature_names <- colnames(high_coupling_df)
  cliffs_deltas <- setNames(numeric(length(feature_names)), feature_names)
  
  for (feature in feature_names) {
    cliffs_deltas[feature] <- cliff.delta(high_coupling_df[[feature]], low_coupling_df[[feature]])$estimate
  }
  
  return(cliffs_deltas)
}

# Filter columns for high and low coupling data frames
RNA_feature_high_coupling_filtered <- RNA_feature_high_coupling %>%
  dplyr::select(all_of(filtered_colnames))

RNA_feature_low_coupling_filtered <- RNA_feature_low_coupling %>%
  dplyr::select(all_of(filtered_colnames))

# Calculate Cliff's delta for all features
cliffs_deltas <- calculate_cliffs_delta_for_all_features(RNA_feature_high_coupling_filtered, RNA_feature_low_coupling_filtered)

# Convert the cliffs_deltas vector into a data frame for ggplot
cliffs_deltas_df <- data.frame(
  Feature = names(cliffs_deltas),
  Delta = cliffs_deltas
)

# Map original feature names to publication-friendly names
cliffs_deltas_df$Feature <- sapply(cliffs_deltas_df$Feature, function(x) {
  if (x %in% names(column_name_mapping)) {
    return(column_name_mapping[x])
  } else {
    return(x)
  }
})

# Sort the data frame by the Delta values in descending order for better visualization
cliffs_deltas_df <- cliffs_deltas_df[order(cliffs_deltas_df$Delta, decreasing = TRUE),]
#cliffs_deltas_df_293 <- cliffs_deltas_df
#save(cliffs_deltas_df_293, file = "../../figure4/compare_293_THP1/data/cliffs_delta_293.RData")

# Create the bar plot
Fig3B <- ggplot(cliffs_deltas_df, aes(x = reorder(Feature, Delta), y = Delta, fill = Delta > 0)) +
  geom_bar(stat = "identity", show.legend = FALSE, width = 0.8) + # Adjust width for border line here
  coord_flip() + # Flip coordinates to make it horizontal; helps in reading feature names
  scale_fill_manual(values = c("TRUE" = "#f768a1", "FALSE" = "#2b8cbe")) +
  labs(title = "Cliff's Delta Values for RNA Features",
       x = "Features",
       y = "Cliff's Delta") +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )



library(ranger)
library(pROC)
library(gbm)
library(rpart)
library(xgboost)
library(ggplot2)
library(dplyr)

# Assuming RNA_feature_coupling_clean is your dataset
# Make sure your data is loaded here
load("data/high_low_coupling_2000_gene.RData")
RNA_feature_low_coupling <- RNA_feature_low_coupling %>%
  dplyr::select(all_of(filtered_colnames))

RNA_feature_high_coupling <- RNA_feature_high_coupling %>%
  dplyr::select(all_of(filtered_colnames))

RNA_feature_low_coupling$coupling <- "low"
RNA_feature_high_coupling$coupling <- "high"

RNA_feature_coupling_clean <- rbind(RNA_feature_low_coupling, RNA_feature_high_coupling)
RNA_feature_coupling_clean$`3utr_pentamer_count` <- NULL
# Prepare the data by converting 'coupling' from 'high'/'low' to 0/1
RNA_feature_coupling_clean$coupling <- ifelse(RNA_feature_coupling_clean$coupling == "high", 0, 1)

# Check for NA values across the dataset
na_columns <- sapply(RNA_feature_coupling_clean, function(x) any(is.na(x)))
na_columns <- na_columns[na_columns]

# Check for NaN values across the dataset
nan_columns <- sapply(RNA_feature_coupling_clean, function(x) any(is.nan(x)))
nan_columns <- nan_columns[nan_columns]

# Check for Inf or -Inf values across the dataset
inf_columns <- sapply(RNA_feature_coupling_clean, function(x) any(is.infinite(x)))
inf_columns <- inf_columns[inf_columns]

# Print the results
if (length(na_columns) > 0) {
  cat("Columns with NA values:\n")
  print(names(na_columns))
} else {
  cat("No columns with NA values.\n")
}

if (length(nan_columns) > 0) {
  cat("Columns with NaN values:\n")
  print(names(nan_columns))
} else {
  cat("No columns with NaN values.\n")
}

if (length(inf_columns) > 0) {
  cat("Columns with Inf or -Inf values:\n")
  print(names(inf_columns))
} else {
  cat("No columns with Inf or -Inf values.\n")
}

# Select features and the target variable
RNA_features <- RNA_feature_coupling_clean[, !(names(RNA_feature_coupling_clean) %in% c("gene_id", "gene_name"))]

# Normalize the features (excluding the target variable 'coupling')
feature_columns <- setdiff(names(RNA_features), "coupling")
RNA_features[, feature_columns] <- scale(RNA_features[, feature_columns])

# Split the data into training and testing sets
set.seed(123)  # for reproducibility
training_indices <- sample(1:nrow(RNA_features), 0.7 * nrow(RNA_features))
train_data <- RNA_features[training_indices, ]
test_data <- RNA_features[-training_indices, ]

# Convert 'coupling' to numeric and handle NAs
train_data <- na.omit(train_data)
train_data$coupling <- as.numeric(train_data$coupling)
test_data <- na.omit(test_data)
test_data$coupling <- as.numeric(test_data$coupling)

names(train_data) <- make.names(names(train_data), unique = TRUE)
names(test_data) <- make.names(names(test_data), unique = TRUE)

# Prepare data for xgboost
xgb_data <- xgb.DMatrix(data = as.matrix(train_data[, -which(names(train_data) == "coupling")]), label = train_data$coupling)

# Train XGBoost model
xgb_model <- xgboost(data = xgb_data, objective = "binary:logistic", nrounds = 500, verbose = 0)

# Train the other models
glm_model <- glm(coupling ~ ., data = train_data, family = binomial(link = "logit"))
rf_model <- ranger(coupling ~ ., data = train_data, probability = TRUE)
gbm_model <- gbm(coupling ~ ., data = train_data, distribution = "bernoulli", n.trees = 500)
rpart_model <- rpart(coupling ~ ., data = train_data, method = "class")

# Initialize a data frame to store the ROC data
roc_data <- data.frame()

# List of models for looping
models <- list("RF" = rf_model, "GLM" = glm_model, "GBM" = gbm_model, "XGB" = xgb_model, "DT" = rpart_model)

for(model_name in names(models)) {
  model <- models[[model_name]]
  
  # Predict using the model
  if(model_name == "XGB") {
    pred <- predict(model, as.matrix(test_data[, -which(names(test_data) == "coupling")]))
  } else if(model_name == "RF") {
    pred <- predict(model, data = test_data, type = "response")$predictions[, 2]
  } else if(model_name == "GBM") {
    pred <- predict(model, newdata = test_data, n.trees = 500, type = "response")
  } else if(model_name == "DT") {
    pred <- predict(model, newdata = test_data, type = "prob")[, 2]
  } else {
    pred <- predict(model, newdata = test_data, type = "response")
  }
  
  # Calculate ROC and AUC
  roc_res <- roc(response = test_data$coupling, predictor = pred)
  auc_res <- auc(roc_res)
  
  # Prepare the data for ggplot
  roc_data <- rbind(roc_data, data.frame(
    model = model_name,
    sensitivity = roc_res$sensitivities,
    specificity = 1 - roc_res$specificities,
    auc = rep(auc_res, length(roc_res$sensitivities))
  ))
}


auc_labels <- roc_data %>%
  group_by(model) %>%
  summarize(
    avg_specificity = mean(specificity),
    avg_sensitivity = mean(sensitivity),
    auc = dplyr::first(auc)  # Use dplyr::first() to avoid conflicts
  )

auc_label_positions <- data.frame(
  model = c("RF", "GLM", "GBM", "XGB", "DT"),
  x = rep(0.8, 5),  # All labels will have the same x-coordinate
  y = seq(0.1, 0.5, length.out = 5)  # Distribute labels vertically
)

# Combine AUC values with the label positions
auc_labels <- merge(auc_label_positions, auc_labels, by = "model")


# Plot ROC curves using ggplot2
Fig3C <- ggplot(data = roc_data, aes(x = specificity, y = sensitivity, color = model)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  geom_text(data = auc_labels, aes(label = paste("AUC =", round(auc, 3)), x = x, y = y), hjust = 0, vjust = 0) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  ggtitle("ROC Curves") +
  scale_color_manual(values = c("#6086ff", "#4cb7a5", "#f29e8e", "#ff0051", "#062ede")) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )+
  theme(legend.title = element_blank())


# Train the Random Forest model with importance calculation enabled
rf_model <- ranger(coupling ~ ., data = train_data, probability = TRUE, importance = 'impurity')

# Extract feature importance
importance_rf <- importance(rf_model)

# Convert to a data frame for ggplot
# Create a data frame from the importance vector
importance_df <- data.frame(
  Feature = names(importance_rf),
  Importance = importance_rf
)
column_name_mapping <- c(
  "transcript_length" = "Transcript Length",
  "5utr_gc_content" = "5' UTR GC Content",
  "5utr_length" = "5' UTR Length",
  "5utrcap_MFE" = "5' UTR Cap MFE",
  "5utr_utr5_cap_gc_content" = "5' UTR Cap GC Content",
  "utr5_MFE_nupack" = "5' UTR MFE Nupack",
  "start_kozak_score" = "Start Kozak Score",
  "start_uORF_count" = "Start uORF Count",
  "cds_gc_content" = "CDS GC Content",
  "cds_length" = "CDS Length",
  "cds_MFE_min" = "CDS MFE Min",
  "cds_avg_codon_freq" = "CDS Avg Codon Frequency",
  "cds_min_codon_freq" = "CDS Min Codon Frequency",
  "3utr_gc_content" = "3' UTR GC Content",
  "3utr_length" = "3' UTR Length",
  "3utr_MFE_min" = "3' UTR MFE Min",
  "3utr_au_element_count" = "3' UTR AU Element Count",
  "3utr_au_element_frac" = "3' UTR AU Element Fraction",
  "3utr_max_au_length" = "3' UTR Max AU Length",
  "exons_gc_content" = "Exons GC Content",
  "exons_length" = "Exons Length",
  "exons_exonct" = "Exons Exon Count",
  "introns_gc_content" = "Introns GC Content",
  "introns_length" = "Introns Length",
  "m6A_5utr_average_level" = "m6A 5' UTR Avg Level",
  "m6A_cds_average_level" = "m6A CDS Avg Level",
  "m6A_3utr_average_level" = "m6A 3' UTR Avg Level",
  "m6A_exons_average_level" = "m6A Exons Avg Level",
  "m6A_introns_average_level" = "m6A Introns Avg Level",
  "RNA_half_life" = "RNA Half-Life"
)

# Function to map feature names considering prefix issues
map_feature_name <- function(feature) {
  clean_feature <- gsub("^X", "", feature)
  if (clean_feature %in% names(column_name_mapping)) {
    return(column_name_mapping[[clean_feature]])
  } else {
    return(feature)
  }
}


# Update feature names in the Random Forest importance plot
importance_df <- importance_df[order(-importance_df$Importance),]

# Map original feature names to publication-friendly names
importance_df$Feature <- sapply(importance_df$Feature, map_feature_name)

# Plotting feature importance using ggplot2
FigS3E <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#6baed6") + # Set bar fill color here
  coord_flip() +  # Flips the axes for better visualization of feature names
  xlab("Features") +
  ylab("Importance") +
  ggtitle("Feature Importance in Random Forest Model") +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

# Update feature names in the XGBoost importance plot
importance_xgb <- xgb.importance(model = xgb_model)
importance_df_xgb <- data.frame(
  Feature = importance_xgb$Feature,
  Importance = importance_xgb$Gain
)

# Ordering the features by importance
importance_df_xgb <- importance_df_xgb[order(-importance_df_xgb$Importance),]

# Map original feature names to publication-friendly names
importance_df_xgb$Feature <- sapply(importance_df_xgb$Feature, map_feature_name)

# Plotting feature importance using ggplot2
Fig3E <- ggplot(importance_df_xgb, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#6baed6") + # Set bar fill color here
  coord_flip() +  # Flips the axes for better visualization of feature names
  xlab("Features") +
  ylab("Importance") +
  ggtitle("Feature Importance in XGBoost Model") +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )



library(ranger)
library(gbm)
library(rpart)
library(xgboost)
library(ggplot2)
library(dplyr)
library(caret) # For scaling

# Assuming the dataset is loaded
load("data/RNA_feature_auc_diff.RData")

RNA_feature_filtered <- RNA_feature_auc_diff %>%
  dplyr::select(all_of(filtered_colnames))

RNA_feature_coupling_clean <- RNA_feature_filtered
RNA_feature_coupling_clean$coupling <- RNA_feature_auc_diff$coupling

# Scale 'coupling' to a range of 0-1
RNA_feature_coupling_clean$coupling <- (RNA_feature_coupling_clean$coupling - min(RNA_feature_coupling_clean$coupling)) / (max(RNA_feature_coupling_clean$coupling) - min(RNA_feature_coupling_clean$coupling))

# Check for NA values across the dataset
na_columns <- sapply(RNA_feature_coupling_clean, function(x) any(is.na(x)))
na_columns <- na_columns[na_columns]

# Check for NaN values across the dataset
nan_columns <- sapply(RNA_feature_coupling_clean, function(x) any(is.nan(x)))
nan_columns <- nan_columns[nan_columns]

# Check for Inf or -Inf values across the dataset
inf_columns <- sapply(RNA_feature_coupling_clean, function(x) any(is.infinite(x)))
inf_columns <- inf_columns[inf_columns]

# Print the results
if (length(na_columns) > 0) {
  cat("Columns with NA values:\n")
  print(names(na_columns))
} else {
  cat("No columns with NA values.\n")
}

if (length(nan_columns) > 0) {
  cat("Columns with NaN values:\n")
  print(names(nan_columns))
} else {
  cat("No columns with NaN values.\n")
}

if (length(inf_columns) > 0) {
  cat("Columns with Inf or -Inf values:\n")
  print(names(inf_columns))
} else {
  cat("No columns with Inf or -Inf values.\n")
}

# Select features and the target variable
RNA_features <- RNA_feature_coupling_clean[, !(names(RNA_feature_coupling_clean) %in% c("gene_id", "gene_name"))]

# Normalize the features (excluding the target variable 'coupling')
feature_columns <- setdiff(names(RNA_features), "coupling")
RNA_features[, feature_columns] <- scale(RNA_features[, feature_columns])

# Split the data into training and testing sets
set.seed(123)  # for reproducibility
training_indices <- sample(1:nrow(RNA_features), 0.7 * nrow(RNA_features))
train_data <- RNA_features[training_indices, ]
test_data <- RNA_features[-training_indices, ]

# Ensure column names are valid for models (important for ranger)
names(train_data) <- make.names(names(train_data), unique = TRUE)
names(test_data) <- make.names(names(test_data), unique = TRUE)

# Model Training for Regression
# XGBoost for regression
xgb_data <- xgb.DMatrix(data = as.matrix(train_data[, -which(names(train_data) == "coupling")]), label = train_data$coupling)

# Prepare data for predictions
xgb_test_matrix <- xgb.DMatrix(data = as.matrix(test_data[, -which(names(test_data) == "coupling")]))

params <- list(objective = "reg:squarederror", 
               eta = 0.1, 
               max_depth = 6)

# Perform cross-validation
cv_results <- xgb.cv(params = params, 
                     data = xgb_data, 
                     nrounds = 1000, 
                     nfold = 5, 
                     showsd = TRUE, 
                     stratified = FALSE, 
                     print_every_n = 10, 
                     early_stopping_rounds = 100, 
                     maximize = FALSE)

# Optimal number of rounds
optimal_nrounds <- cv_results$best_iteration

# Retrain model with optimal nrounds
xgb_model_optimized <- xgboost(data = xgb_data, 
                               params = params, 
                               nrounds = optimal_nrounds, 
                               verbose = 0)

# Predict on the test set with the optimized XGBoost model
optimized_xgb_predictions <- predict(xgb_model_optimized, xgb_test_matrix)

# Calculate the correlation between the optimized predictions and actual values
cor_optimized_xgb <- cor(optimized_xgb_predictions, test_data$coupling)

# Print the correlation coefficient
cat("Correlation coefficient for the optimized XGBoost model:", cor_optimized_xgb, "\n")

# Plotting the correlation with ggplot2
library(ggplot2)

# Create a dataframe for plotting
data_for_plot <- data.frame(Actual = test_data$coupling, Predicted = optimized_xgb_predictions)


library(ggplot2)
library(FNN)  # For k-nearest neighbors

# Calculate correlation and p-value
cor_results <- cor.test(data_for_plot$Actual, data_for_plot$Predicted, method = "pearson")
cor_value <- round(cor_results$estimate, 3)
p_value <- signif(cor_results$p.value, 3)

# Calculate density using k-nearest neighbors
k <- 10  # Number of neighbors to consider
neighbors <- get.knn(data = data_for_plot[, c("Actual", "Predicted")], k = k)
data_for_plot$density <- 1 / rowMeans(neighbors$nn.dist)  # Density = inverse of average distance to neighbors
data_for_plot$density <- data_for_plot$density / max(data_for_plot$density)  # Normalize density

# Custom color palette
custom_colors <- c("#54a5de", "white", "#f15389")

# Plot with density-colored points and regression line
Fig3D_actual_vs_predicted_plot <- ggplot(data_for_plot, aes(x = Actual, y = Predicted, color = density)) +
  geom_point(size = 2, alpha = 0.7) +  # Points colored by density
  geom_smooth(method = "lm", color = "#08306b", linewidth = 1, linetype = "solid", se = FALSE) +  # Regression line
  scale_color_gradientn(colors = custom_colors, values = c(0, 0.2, 1), name = "Density") +  # Custom color bar
  labs(
    x = "Actual Values",
    y = "Predicted Values",
    title = sprintf("Correlation between Actual and Predicted Values\nr = %.3f, P = %.3g", cor_value, p_value)
  ) +
  annotate("text", x = 1, y = 0.5, label = paste0("r = ", cor_value, "\nP = ", p_value),
           hjust = 1, vjust = 1, size = 5, color = "black", fontface = "bold") +  # Correlation annotation
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )



source("code/polysome_seq_functions.R")
load("data/HEK293T_preprocessed.RData")

merge_df_counts_select <- merge_df_counts[merge_df_counts$gene_id %in% genes_logCPM[genes_logCPM$logCPM>0,"gene_id"],]
merge_df_counts_new_select <- merge_df_counts_new[merge_df_counts$gene_id %in% genes_logCPM[genes_logCPM$logCPM>0,"gene_id"],]

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
  theme_classic()


MRL_new_high_coupling <- MRL_new[high_coupling_indices, ]
MRL_old_high_coupling <- MRL_old[high_coupling_indices, ]

MRL_new_low_coupling <- MRL_new[low_coupling_indices, ]
MRL_old_low_coupling <- MRL_old[low_coupling_indices, ]

library(ggplot2)
library(dplyr)
library(tidyr)

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


# Function to convert p-values to star labels
convert_p_to_stars <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

library(dplyr)
library(ggplot2)

# Step 1: Select high coupling genes
merge_df_TC_rate_high <- merge_df_TC_rate %>%
  filter(gene_id %in% high_coupling_genes)

# Step 2: Select specific columns
selected_columns <- c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4",
                      "EGI_frac1", "EGI_frac2", "EGI_frac3", "EGI_frac4",
                      "actD_frac1", "actD_frac2", "actD_frac3", "actD_frac4",
                      "puro_frac1", "puro_frac2", "puro_frac3", "puro_frac4",
                      "STM2457_frac1", "STM2457_frac2", "STM2457_frac3", "STM2457_frac4",
                      "STM2457_12h_frac1", "STM2457_12h_frac2", "STM2457_12h_frac3", "STM2457_12h_frac4",
                      "YTHDF1_frac1", "YTHDF1_frac2", "YTHDF1_frac3", "YTHDF1_frac4",
                      "YTHDF2_frac1", "YTHDF2_frac2", "YTHDF2_frac3", "YTHDF2_frac4",
                      "YTHDF3_frac1", "YTHDF3_frac2", "YTHDF3_frac3", "YTHDF3_frac4")
merge_df_TC_rate_high <- merge_df_TC_rate_high %>%
  dplyr::select(all_of(selected_columns))

# Step 3: Calculate mean values (assuming you need the mean across all genes)
mean_values <- colMeans(merge_df_TC_rate_high, na.rm = TRUE)

# Step 4: Create the plot
# Reshaping data for plotting
mean_values_df <- data.frame(
  Fraction = rep(paste0("frac", 1:4), length(mean_values) / 4),
  Treatment = rep(c("NC_2h", "EGI", "actD", "puro", "STM2457", "STM2457_12h", "YTHDF1", "YTHDF2", "YTHDF3"), each = 4),
  MeanValue = mean_values
)

# Reorder the Treatment factor to put NC_2h first
mean_values_df$Treatment <- factor(mean_values_df$Treatment, levels = c("NC_2h", "EGI", "actD", "puro", "STM2457", "STM2457_12h", "YTHDF1", "YTHDF2", "YTHDF3"))


# Step 1: Select high coupling genes
merge_df_TC_rate_low <- merge_df_TC_rate %>%
  filter(gene_id %in% low_coupling_genes)

# Step 2: Select specific columns
selected_columns <- c("NC_2h_frac1", "NC_2h_frac2", "NC_2h_frac3", "NC_2h_frac4",
                      "EGI_frac1", "EGI_frac2", "EGI_frac3", "EGI_frac4",
                      "actD_frac1", "actD_frac2", "actD_frac3", "actD_frac4",
                      "puro_frac1", "puro_frac2", "puro_frac3", "puro_frac4",
                      "STM2457_frac1", "STM2457_frac2", "STM2457_frac3", "STM2457_frac4",
                      "STM2457_12h_frac1", "STM2457_12h_frac2", "STM2457_12h_frac3", "STM2457_12h_frac4",
                      "YTHDF1_frac1", "YTHDF1_frac2", "YTHDF1_frac3", "YTHDF1_frac4",
                      "YTHDF2_frac1", "YTHDF2_frac2", "YTHDF2_frac3", "YTHDF2_frac4",
                      "YTHDF3_frac1", "YTHDF3_frac2", "YTHDF3_frac3", "YTHDF3_frac4")
merge_df_TC_rate_low <- merge_df_TC_rate_low %>%
  dplyr::select(all_of(selected_columns))

# Step 3: Calculate mean values (assuming you need the mean across all genes)
mean_values <- colMeans(merge_df_TC_rate_low, na.rm = TRUE)

# Step 4: Create the plot
# Reshaping data for plotting
mean_values_df <- data.frame(
  Fraction = rep(paste0("frac", 1:4), length(mean_values) / 4),
  Treatment = rep(c("NC_2h", "EGI", "actD", "puro", "STM2457", "STM2457_12h", "YTHDF1", "YTHDF2", "YTHDF3"), each = 4),
  MeanValue = mean_values
)

# Reorder the Treatment factor to put NC_2h first
mean_values_df$Treatment <- factor(mean_values_df$Treatment, levels = c("NC_2h", "EGI", "actD", "puro", "STM2457", "STM2457_12h", "YTHDF1", "YTHDF2", "YTHDF3"))


NC_cols <- grep("^NC_2h", names(merge_df_counts_new_select), value = TRUE)
STM2457_2h_cols <- grep("^STM2457_f", names(merge_df_counts_new_select), value = TRUE)
STM2457_12h_cols <- grep("^STM2457_12h", names(merge_df_counts_new_select), value = TRUE)
NC_counts <- merge_df_counts_new_select[,NC_cols]
STM2457_2h_counts <- merge_df_counts_new_select[,STM2457_2h_cols]
STM2457_12h_counts <- merge_df_counts_new_select[,STM2457_12h_cols]
NC_counts_norm <- sweep(NC_counts, 2, colSums(NC_counts), "/")
STM2457_2h_counts_norm <- sweep(STM2457_2h_counts, 2, colSums(STM2457_2h_counts), "/")
STM2457_12h_counts_norm <- sweep(STM2457_12h_counts, 2, colSums(STM2457_12h_counts), "/")

NC_counts_norm_mean <- rowMeans(NC_counts_norm)
STM2457_2h_counts_norm_mean <- rowMeans(STM2457_2h_counts_norm)
STM2457_12h_counts_norm_mean <- rowMeans(STM2457_12h_counts_norm)

NC_2h_RL <- calculate_polysome_load_old(merge_df_counts_new_select, "NC_2h")
STM2457_2h_RL <- calculate_polysome_load_old(merge_df_counts_new_select, "STM2457")
STM2457_12h_RL <- calculate_polysome_load_old(merge_df_counts_new_select, "STM2457_12h")

counts_fc_2h <- STM2457_2h_counts_norm_mean/NC_counts_norm_mean
RL_fc_2h <- STM2457_2h_RL/NC_2h_RL

counts_fc_12h <- STM2457_12h_counts_norm_mean/NC_counts_norm_mean
RL_fc_12h <- STM2457_12h_RL/NC_2h_RL

# Create a data frame
data <- data.frame(
  auc_diff = auc_diff,
  STM2457_2h = log10(RL_fc_2h$STM2457),
  STM2457_12h = log10(RL_fc_12h$STM2457)
)


# Source the HEK293T_RNA_features.R
load("data/RNA_features_gene_level_20240418.RData")
RL_fc_2h$gene_id <- merge_df_counts_new_select$gene_id

RL_fc_2h_tail <- RL_fc_2h[log10(RL_fc_2h$STM2457)<= -0.1,]
RL_fc_2h_high <- RL_fc_2h[log10(RL_fc_2h$STM2457)> 0.1,]

RNA_feature_of_gene <- RNA_features_gene_level
RNA_feature_of_gene_select <- RNA_feature_of_gene[RNA_feature_of_gene$gene_id %in% merge_df_counts_new_select$gene_id,]

RL_fc_2h_tail_feature <- merge(RNA_feature_of_gene_select, RL_fc_2h_tail, by.x = "gene_id")
RL_fc_2h_high_feature <- merge(RNA_feature_of_gene_select, RL_fc_2h_high, by.x = "gene_id")

t.test(RL_fc_2h_tail_feature$m6A_cds_average_level, RL_fc_2h_high_feature$m6A_cds_average_level)

RNA_feature_of_gene_select <- RNA_feature_of_gene[RNA_feature_of_gene$gene_id %in% merge_df_counts_new_select$gene_id,]
RL_fc_2h$gene_id <- merge_df_counts_new_select$gene_id
RL_fc_12h$gene_id <- merge_df_counts_new_select$gene_id
RNA_feature_of_gene_select <- merge(RNA_feature_of_gene_select, RL_fc_2h, by.x = "gene_id")
RNA_feature_of_gene_select <- merge(RNA_feature_of_gene_select, RL_fc_12h, by.x = "gene_id")


##########  用不同时间点的new/old mRNA的MRL比值与AUC_diff 比较，先证明两个是正比的
MRL_diff <- MRL_new-MRL_old

MRL_diff_output <- MRL_diff
MRL_diff_output$gene_id <- merge_df_counts_select$gene_id
#write.table(MRL_diff_output, file = "gene_level_MRL_diff.csv", quote = FALSE, sep = ",")

MRL_diff$auc_diff <- auc_diff
MRL_diff <- MRL_diff[(!is.infinite(MRL_diff$NC_30min))&(!is.infinite(MRL_diff$NC_1h))&(!is.infinite(MRL_diff$NC_2h)), ]

NC_2h_RL_new <- calculate_polysome_load_old(merge_df_counts_new_select, "NC_2h")
STM2457_2h_RL_new <- calculate_polysome_load_old(merge_df_counts_new_select, "STM2457")
STM2457_12h_RL_new <- calculate_polysome_load_old(merge_df_counts_new_select, "STM2457_12h")

NC_2h_RL_old <- calculate_polysome_load_old(merge_df_counts_old_select, "NC_2h")
STM2457_2h_RL_old <- calculate_polysome_load_old(merge_df_counts_old_select, "STM2457")
STM2457_12h_RL_old <- calculate_polysome_load_old(merge_df_counts_old_select, "STM2457_12h")

NC_2h_RL_diff <- NC_2h_RL_new - NC_2h_RL_old
STM2457_2h_RL_diff <- STM2457_2h_RL_new - STM2457_2h_RL_old
STM2457_12h_RL_diff <- STM2457_12h_RL_new - STM2457_12h_RL_old


MRL_STM2457 <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_2h", "STM2457", "STM2457_12h"))
MRL_STM2457_new <- MRL_STM2457[[1]]
MRL_STM2457_old <- MRL_STM2457[[2]]

MRL_STM2457_diff <- MRL_STM2457_new - MRL_STM2457_old

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Fit a linear model
model1 <- lm(STM2457 ~ NC_2h, data = MRL_STM2457_diff)

# Calculate residuals and standard deviation of residuals
MRL_STM2457_diff$residuals1 <- resid(model1)
std_dev1 <- sd(MRL_STM2457_diff$residuals1)

# Identify upper and lower outliers
MRL_STM2457_diff <- MRL_STM2457_diff %>%
  mutate(
    rank_residuals1 = rank(-abs(residuals1)),  # Rank by absolute values of residuals
    highlight_upper1 = rank(-residuals1) <= 725,  # Top 100 positive residuals
    highlight_lower1 = rank(residuals1) <= 725    # Top 100 negative residuals
  )

MRL_STM2457_diff$gene_id <- merge_df_counts_new_select$gene_id

# Create the plot and highlight upper and lower outliers with specified colors
Fig3F_plot1 <- ggplot(MRL_STM2457_diff, aes(x = NC_2h, y = STM2457)) +
  geom_point(aes(color = ifelse(highlight_upper1, "#e95966", ifelse(highlight_lower1, "#00bad5", "grey")), alpha = 0.01)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "NC 2h", y = "STM2457", title = "NC 2h vs STM2457 with Outliers Highlighted") +
  scale_color_identity() +
  coord_cartesian(xlim = c(-7.5, 5.5), ylim = c(-7.5, 5.5)) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )



MRL_STM2457_diff$gene_id <- merge_df_counts_new_select$gene_id
RNA_feature_of_gene_select_STM2457_2h_high <- RNA_feature_of_gene_select[RNA_feature_of_gene_select$gene_id %in% MRL_STM2457_diff[MRL_STM2457_diff$highlight_upper1,]$gene_id,]
RNA_feature_of_gene_select_STM2457_2h_low <- RNA_feature_of_gene_select[RNA_feature_of_gene_select$gene_id %in% MRL_STM2457_diff[MRL_STM2457_diff$highlight_lower1,]$gene_id,]

# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Perform t-tests
t_test_result_cds <- t.test(RNA_feature_of_gene_select_STM2457_2h_high$m6A_cds_average_level, RNA_feature_of_gene_select_STM2457_2h_low$m6A_cds_average_level)
t_test_result_utr5 <- t.test(RNA_feature_of_gene_select_STM2457_2h_high$m6A_5utr_average_level, RNA_feature_of_gene_select_STM2457_2h_low$m6A_5utr_average_level)
t_test_result_utr3 <- t.test(RNA_feature_of_gene_select_STM2457_2h_high$m6A_3utr_average_level, RNA_feature_of_gene_select_STM2457_2h_low$m6A_3utr_average_level)

# Create combined datasets for plotting
data_for_plot_cds <- rbind(
  data.frame(Group = "Up", m6A_level_cds = RNA_feature_of_gene_select_STM2457_2h_high$m6A_cds_average_level),
  data.frame(Group = "Down", m6A_level_cds = RNA_feature_of_gene_select_STM2457_2h_low$m6A_cds_average_level)
)

data_for_plot_utr5 <- rbind(
  data.frame(Group = "Up", m6A_level_5utr = RNA_feature_of_gene_select_STM2457_2h_high$m6A_5utr_average_level),
  data.frame(Group = "Down", m6A_level_5utr = RNA_feature_of_gene_select_STM2457_2h_low$m6A_5utr_average_level)
)

data_for_plot_utr3 <- rbind(
  data.frame(Group = "Up", m6A_level_3utr = RNA_feature_of_gene_select_STM2457_2h_high$m6A_3utr_average_level),
  data.frame(Group = "Down", m6A_level_3utr = RNA_feature_of_gene_select_STM2457_2h_low$m6A_3utr_average_level)
)

# Create the box plots
# CDS Plot
box_plot_cds <- ggplot(data_for_plot_cds, aes(x = Group, y = m6A_level_cds, fill = Group)) +
  geom_boxplot() +
  labs(title = "m6A level CDS", x = "Group", y = "m6A Level") +
  scale_fill_manual(values = c("Up" = "#e95966", "Down" = "#00bad5")) +
  annotate("text", x = 1.5, y = 0.8, 
           label = paste("p =", t_test_result_cds$p.value), 
           vjust = 1.5, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )


# UTR5 Plot
box_plot_utr5 <- ggplot(data_for_plot_utr5, aes(x = Group, y = m6A_level_5utr, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "m6A Level 5'UTR", x = "Group", y = "m6A Level UTR5") +
  scale_fill_manual(values = c("Up" = "#e95966", "Down" = "#00bad5")) +
  annotate("text", x = 1.5, y = 0.8, 
           label = paste("p =", t_test_result_utr5$p.value), 
           vjust = 1.5, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

# UTR3 Plot
box_plot_utr3 <- ggplot(data_for_plot_utr3, aes(x = Group, y = m6A_level_3utr, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "m6A Level UTR3", x = "Group", y = "m6A Level") +
  scale_fill_manual(values = c("Up" = "#e95966", "Down" = "#00bad5")) +
  annotate("text", x = 1.5, y = 0.8, 
           label = paste("p =", t_test_result_utr3$p.value), 
           vjust = 1.5, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )


Fig3G <- grid.arrange(box_plot_utr5, box_plot_cds, box_plot_utr3, ncol = 3)


RNA_feature_of_gene_select_STM2457_2h_high <- RNA_feature_of_gene_select[RNA_feature_of_gene_select$gene_id %in% MRL_STM2457_diff[MRL_STM2457_diff$highlight_upper1,]$gene_id,]
RNA_feature_of_gene_select_STM2457_2h_low <- RNA_feature_of_gene_select[RNA_feature_of_gene_select$gene_id %in% MRL_STM2457_diff[MRL_STM2457_diff$highlight_lower1,]$gene_id,]


MRL_YTHDF <- calculate_polysome_load(merge_df_counts_new_select, merge_df_counts_old_select, c("NC_2h", "YTHDF1", "YTHDF2", "YTHDF3"))
MRL_YTHDF_new <- MRL_YTHDF[[1]]
MRL_YTHDF_old <- MRL_YTHDF[[2]]

MRL_YTHDF_diff <- MRL_YTHDF_new - MRL_YTHDF_old
MRL_YTHDF_diff[is.na(MRL_YTHDF_diff)] <- 0
# Load necessary libraries
library(ggplot2)
library(dplyr)

#YTHDF1
# Fit a linear model
model1 <- lm(YTHDF1 ~ NC_2h, data = MRL_YTHDF_diff)
summary(model1)
# Calculate residuals and standard deviation of residuals
MRL_YTHDF_diff$residuals1 <- resid(model1)
std_dev1 <- sd(MRL_YTHDF_diff$residuals1)

# Identify upper and lower outliers
MRL_YTHDF_diff <- MRL_YTHDF_diff %>%
  mutate(
    rank_residuals1 = rank(-abs(residuals1)),  # Rank by absolute values of residuals
    highlight_upper1 = rank(-residuals1) <= 725,  # Top 100 positive residuals
    highlight_lower1 = rank(residuals1) <= 725    # Top 100 negative residuals
  )

# Create the plot and highlight upper and lower outliers with specified colors
Fig3H_plot1 <- ggplot(MRL_YTHDF_diff, aes(x = NC_2h, y = YTHDF1)) +
  geom_point(aes(color = ifelse(highlight_upper1, "#e95966", ifelse(highlight_lower1, "#00bad5", "grey")), alpha = 0.01)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "NC 2h", y = "YTHDF1", title = "NC 2h vs YTHDF1 with Outliers Highlighted") +
  scale_color_identity() +
  coord_cartesian(xlim = c(-7.5, 5.5), ylim = c(-7.5, 5.5)) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )



#YTHDF2
# Fit a second linear model
model2 <- lm(YTHDF2 ~ NC_2h, data = MRL_YTHDF_diff)
summary(model2)
# Calculate residuals and standard deviation of residuals for the second model
MRL_YTHDF_diff$residuals2 <- resid(model2)
std_dev2 <- sd(MRL_YTHDF_diff$residuals2)

MRL_YTHDF_diff <- MRL_YTHDF_diff %>%
  mutate(
    rank_residuals2 = rank(-abs(residuals2)),  # Rank by absolute values of residuals
    highlight_upper2 = rank(-residuals2) <= 725,  # Top 100 positive residuals
    highlight_lower2 = rank(residuals2) <= 725    # Top 100 negative residuals
  )

# Create the plot for the second model and highlight upper and lower outliers with specified colors
FigS3F_plot2 <- ggplot(MRL_YTHDF_diff, aes(x = NC_2h, y = YTHDF2)) +
  geom_point(aes(color = ifelse(highlight_upper2, "#e95966", ifelse(highlight_lower2, "#00bad5", "grey")), alpha = 0.01)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "NC 2h", y = "YTHDF2", title = "NC 2h vs YTHDF2 with Outliers Highlighted") +
  scale_color_identity() +
  coord_cartesian(xlim = c(-7.5, 5.5), ylim = c(-7.5, 5.5)) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )



#YTHDF3
# Fit a second linear model
model3 <- lm(YTHDF3 ~ NC_2h, data = MRL_YTHDF_diff)
summary(model3)
# Calculate residuals and standard deviation of residuals for the second model
MRL_YTHDF_diff$residuals3 <- resid(model3)
std_dev3 <- sd(MRL_YTHDF_diff$residuals3)

# Identify upper and lower outliers for the second model

MRL_YTHDF_diff <- MRL_YTHDF_diff %>%
  mutate(
    rank_residuals3 = rank(-abs(residuals3)),  # Rank by absolute values of residuals
    highlight_upper3 = rank(-residuals3) <= 725,  # Top 100 positive residuals
    highlight_lower3 = rank(residuals3) <= 725    # Top 100 negative residuals
  )
# Create the plot for the second model and highlight upper and lower outliers with specified colors
FigS3G_plot3 <- ggplot(MRL_YTHDF_diff, aes(x = NC_2h, y = YTHDF3)) +
  geom_point(aes(color = ifelse(highlight_upper3, "#e95966", ifelse(highlight_lower3, "#00bad5", "grey")), alpha = 0.01)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Add diagonal line
  labs(x = "NC 2h", y = "YTHDF3", title = "NC 2h vs YTHDF3 with Outliers Highlighted") +
  scale_color_identity() +
  coord_cartesian(xlim = c(-7.5, 5.5), ylim = c(-7.5, 5.5)) +
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )



# YTHDF1
MRL_YTHDF_diff$gene_id <- merge_df_counts_new_select$gene_id
RNA_feature_of_gene_select_YTHDF1_2h_high <- RNA_feature_of_gene_select[RNA_feature_of_gene_select$gene_id %in% MRL_YTHDF_diff[MRL_YTHDF_diff$highlight_upper1,]$gene_id,]
RNA_feature_of_gene_select_YTHDF1_2h_low <- RNA_feature_of_gene_select[RNA_feature_of_gene_select$gene_id %in% MRL_YTHDF_diff[MRL_YTHDF_diff$highlight_lower1,]$gene_id,]


# Perform t-tests
t_test_result_cds <- t.test(RNA_feature_of_gene_select_YTHDF1_2h_high$m6A_cds_average_level, RNA_feature_of_gene_select_YTHDF1_2h_low$m6A_cds_average_level)
t_test_result_utr5 <- t.test(RNA_feature_of_gene_select_YTHDF1_2h_high$m6A_5utr_average_level, RNA_feature_of_gene_select_YTHDF1_2h_low$m6A_5utr_average_level)
t_test_result_utr3 <- t.test(RNA_feature_of_gene_select_YTHDF1_2h_high$m6A_3utr_average_level, RNA_feature_of_gene_select_YTHDF1_2h_low$m6A_3utr_average_level)

# Create combined datasets for plotting
data_for_plot_cds <- rbind(
  data.frame(Group = "Up", m6A_level_cds = RNA_feature_of_gene_select_YTHDF1_2h_high$m6A_cds_average_level),
  data.frame(Group = "Down", m6A_level_cds = RNA_feature_of_gene_select_YTHDF1_2h_low$m6A_cds_average_level)
)

data_for_plot_utr5 <- rbind(
  data.frame(Group = "Up", m6A_level_5utr = RNA_feature_of_gene_select_YTHDF1_2h_high$m6A_5utr_average_level),
  data.frame(Group = "Down", m6A_level_5utr = RNA_feature_of_gene_select_YTHDF1_2h_low$m6A_5utr_average_level)
)

data_for_plot_utr3 <- rbind(
  data.frame(Group = "Up", m6A_level_3utr = RNA_feature_of_gene_select_YTHDF1_2h_high$m6A_3utr_average_level),
  data.frame(Group = "Down", m6A_level_3utr = RNA_feature_of_gene_select_YTHDF1_2h_low$m6A_3utr_average_level)
)

# Create the box plots
box_plot_cds <- ggplot(data_for_plot_cds, aes(x = Group, y = m6A_level_cds, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "m6A level CDS", x = "Group", y = "m6A Level") +
  scale_fill_manual(values = c("Up" = "#e95966", "Down" = "#00bad5")) +
  annotate("text", x = 1.5, y = 0.8, 
           label = paste("p =", t_test_result_cds$p.value), 
           vjust = 1.5, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

box_plot_utr5 <- ggplot(data_for_plot_utr5, aes(x = Group, y = m6A_level_5utr, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "m6A Level 5'UTR", x = "Group", y = "m6A Level UTR5") +
  scale_fill_manual(values = c("Up" = "#e95966", "Down" = "#00bad5")) +
  annotate("text", x = 1.5, y = 0.8, 
           label = paste("p =", t_test_result_utr5$p.value), 
           vjust = 1.5, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

box_plot_utr3 <- ggplot(data_for_plot_utr3, aes(x = Group, y = m6A_level_3utr, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "m6A Level UTR3", x = "Group", y = "m6A Level") +
  scale_fill_manual(values = c("Up" = "#e95966", "Down" = "#00bad5")) +
  annotate("text", x = 1.5, y = 0.8, 
           label = paste("p =", t_test_result_utr3$p.value), 
           vjust = 1.5, size = 3)+
  theme_bw() +  # Use a minimal theme
  theme(
    legend.position = "none",  # Remove legend
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),  # Change axis numbers to black and larger
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", linewidth = 1)  # Change frame color to black
  )

# Output to PDF
#pdf(file = "result_new_20241116/m6A level of YTHDF1 up and down.pdf", width = 6, height = 3.5)
Fig3I <- grid.arrange(box_plot_utr5, box_plot_cds, box_plot_utr3, ncol = 3)
#dev.off()



library(readr)
library(dplyr)
library(stringr)

# Example usage:
renamelist_NC30min <- '/storage/cbq/cbq_backup2/disk1_old/NGS_data/20230328_30min_1h_2h/renamelist_NC30min.csv'
renamelist_NC1h <- '/storage/cbq/cbq_backup2/disk1_old/NGS_data/20230328_30min_1h_2h/renamelist_NC1h.csv'
renamelist_NC2h <- '/storage/cbq/cbq_backup2/disk1_old/NGS_data/20230606_2h_newdata/renamelist_NC2h_lcd.csv'

NC1h_yeast_all_counts <- "/storage/cbq/cbq_backup2/disk3_old/NGS_data/20231005_integrate_data/data/time_point/NC1h/NC1h_yeast_all_counts.csv"
NC2h_yeast_all_counts <- "/storage/cbq/cbq_backup2/disk3_old/NGS_data/20231005_integrate_data/data/time_point/NC2h_new2/NC2h_yeast_all_counts.csv" # This sample use: disk1/cbq/NGS_data/20230606_2h_newdata/NC2h_lcd_*
NC30min_yeast_all_counts <- "/storage/cbq/cbq_backup2/disk3_old/NGS_data/20231005_integrate_data/data/time_point/NC30min/NC30min_yeast_all_counts.csv"

source("R/functions_HEK293_new_old_isoform_diff_2024.3.7.R")

all_samples_df_30min <- process_normalize_counts(renamelist_NC30min, NC30min_yeast_all_counts, '/storage/cbq/cbq_backup2/disk3_old/NGS_data/20240229_isoform_mapping/new_old_isoform/mapping_NC30min/cuffnorm_counts')
all_samples_df_NC1h <- process_normalize_counts(renamelist_NC1h, NC1h_yeast_all_counts, '/storage/cbq/cbq_backup2/disk3_old/NGS_data/20240229_isoform_mapping/new_old_isoform/mapping_NC1h/cuffnorm_counts')
all_samples_df_NC2h <- process_normalize_counts(renamelist_NC2h, NC2h_yeast_all_counts, '/storage/cbq/cbq_backup2/disk3_old/NGS_data/20240229_isoform_mapping/new_old_isoform/mapping_NC2h_lcd/cuffnorm_counts')


load("data/HEK293_transcript_identify_clear.RData")
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

source("R/polysome_seq_functions_20230703.R")
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
isoform_expression_select <- isoform_expression[isoform_expression$expression>100,]

isoform_counts_new_select <- cbind(NC30min_new, NC1h_new, NC2h_new)
isoform_counts_old_select <- cbind(NC30min_old, NC1h_old, NC2h_old)

isoform_level_MRLs <- calculate_polysome_load(isoform_counts_new_select, isoform_counts_old_select, c("NC30min", "NC1h", "NC2h"))

MRL_new <- isoform_level_MRLs$new_polysome_load
MRL_old <- isoform_level_MRLs$old_polysome_load

library(pracma)
time_points <- c(0.5, 1, 2)  # 30 min, 60 min (1h), 120 min (2h)
calculate_auc <- function(row) {
  trapz(time_points, row)
}

auc_MRL_new <- apply(MRL_new, 1, calculate_auc)
auc_MRL_old <- apply(MRL_old, 1, calculate_auc)

auc_diff <- auc_MRL_new - auc_MRL_old

isoform_level_MRL_diff <- isoform_level_MRLs$new_polysome_load - isoform_level_MRLs$old_polysome_load
isoform_level_MRL_diff$transcript_id <- NC30min_merge_reps$tracking_id
isoform_level_MRL_diff$auc_diff <- auc_diff

isoform_level <- merge(isoform_level_MRL_diff, isoform_expression_select, by = "transcript_id")

load("data/RNA_features/RNA_features_transcript_level.RData")
RNA_features_transcript_level <- RNA_features_transcript_level %>%
  dplyr::rename(transcript_id = transcriptID)
RNA_features_transcript_level_select_2 <- merge(RNA_features_transcript_level, isoform_level, by = "transcript_id")



compare_feature_to_MRL_diff <- function(data, feature_of_interest) {
  
  # Step 1: Filter for genes with multiple isoforms and handle NA values in numeric columns
  multi_iso_genes <- data
  
  # Step 2: Select the top 1 (highest) and bottom 1 (lowest) isoforms based on NC2h values for each gene
  high_low_data <- multi_iso_genes %>%
    group_by(gene_id) %>%
    arrange(desc(NC2h), .by_group = TRUE) %>%
    dplyr::slice(c(1, n())) %>%  # Select top 1 (highest) and bottom 1 (lowest) NC2h values
    mutate(
      kinetics_level = ifelse(row_number() == 1, "High", "Low")  # Label highest as High, lowest as Low
    ) %>%
    ungroup() %>%
    dplyr::select(gene_id, kinetics_level, !!sym(feature_of_interest))  # Select relevant columns
  
  # Calculate the count of High and Low values
  count_data <- high_low_data %>%
    group_by(kinetics_level) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Create modified labels for High and Low with counts
  high_label <- paste0("High (n=", count_data$count[count_data$kinetics_level == "High"], ")")
  low_label <- paste0("Low (n=", count_data$count[count_data$kinetics_level == "Low"], ")")
  
  # Step 3: Perform an unpaired t-test between High and Low groups
  high_values <- high_low_data %>% filter(kinetics_level == "High") %>% pull(.data[[feature_of_interest]])
  low_values <- high_low_data %>% filter(kinetics_level == "Low") %>% pull(.data[[feature_of_interest]])
  
  # Ensure there are enough values for the test
  if (length(high_values) > 1 && length(low_values) > 1) {
    t_test_result <- t.test(high_values, low_values, paired = TRUE)
    p_value <- t_test_result$p.value
  } else {
    p_value <- NA
  }
  
  # Step 4: Plot the data with feature levels on the y-axis and high/low MRL_diff_NC2h on the x-axis
  plot <- ggplot(high_low_data, aes(x = kinetics_level, y = .data[[feature_of_interest]], fill = kinetics_level)) +
    geom_boxplot() +
    scale_x_discrete(labels = c("High","Low")) +  # Apply custom labels with counts
    scale_fill_manual(values = c("High" = "#8382b1", "Low" = "#90b987")) +
    labs(
      x = paste0(high_label,',',low_label),
      y = feature_of_interest,
      title = feature_of_interest,
      subtitle = paste("Paired p-value:", ifelse(!is.na(p_value), signif(p_value, digits = 3), "N/A"))
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
  
  # Return the plot and the p-value
  return(list(feature = feature_of_interest, p_value = p_value, plot = plot))
}



# List of features to compare
features <- c(
  "transcript_length", "5utr_gc_content", "5utr_length", "5utr_MFE_min", "5utrcap_MFE", 
  "5utr_utr5_cap_gc_content", "start_kozak_score", "start_uORF_count", "start_uORF_overlap",
  "cds_gc_content", "cds_length", "cds_MFE_min", "cds_avg_codon_freq", "cds_min_codon_freq",
  "3utr_gc_content", "3utr_length", "3utr_MFE_min", "3utr_au_element_count", 
  "3utr_au_element_frac", "3utr_max_au_length", "exons_gc_content", "exons_length", 
  "exons_exonct", "introns_gc_content", "introns_length", "m6A_5utr_site_num", 
  "m6A_5utr_average_level", "m6A_cds_site_num", "m6A_cds_average_level", "m6A_3utr_site_num", 
  "m6A_3utr_average_level", "m6A_exons_site_num", "m6A_exons_average_level", 
  "m6A_introns_site_num", "m6A_introns_average_level", "RNA_half_life"
)

RNA_features_transcript_level_select_2 <- RNA_features_transcript_level_select_2 %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))


RNA_features_transcript_level_filtered <- RNA_features_transcript_level_select_2 %>%
  group_by(gene_id) %>%
  filter(n() > 1) %>%  # Keep only genes with more than 2 transcripts
  ungroup()

# Step 2: Compute the median auc_diff across all transcripts
median_auc_diff <- median(RNA_features_transcript_level_filtered$auc_diff, na.rm = TRUE)

# Step 3: Classify transcripts as 'high' or 'low' based on the median auc_diff
RNA_features_transcript_level_filtered <- RNA_features_transcript_level_filtered %>%
  mutate(
    class = ifelse(auc_diff > median_auc_diff, "high", "low")  # Classify as 'high' or 'low'
  )


# Load necessary libraries
library(dplyr)
library(caret)
library(xgboost)
library(ranger)
library(gbm)
library(rpart)
library(pROC)
library(ggplot2)

# Data preparation: Start from RNA_features_transcript_level_multi_iso
RNA_features <- RNA_features_transcript_level_filtered %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%  # Replace NA values with 0
  dplyr::select(-c(transcript_id, gene_name, gene_id, NC30min, NC1h, NC2h, auc_diff, expression)) %>%
  mutate(class = factor(class, levels = c("low", "high")))  # Ensure class levels are correct

# Normalize numeric features (excluding the class column)
feature_columns <- setdiff(names(RNA_features), "class")
RNA_features[, feature_columns] <- scale(RNA_features[, feature_columns])

# Sanitize column names
names(RNA_features) <- make.names(names(RNA_features), unique = TRUE)

# Split into training and testing sets
set.seed(123)
training_indices <- sample(1:nrow(RNA_features), 0.7 * nrow(RNA_features))
train_data <- RNA_features[training_indices, ]
test_data <- RNA_features[-training_indices, ]

# Convert 'class' to numeric for models that require binary encoding
train_data_gbm <- train_data
train_data_gbm$class <- as.numeric(train_data_gbm$class) - 1  # 0 = "low", 1 = "high"

test_data_gbm <- test_data
test_data_gbm$class <- as.numeric(test_data_gbm$class) - 1

# Prepare data for XGBoost
xgb_data <- xgb.DMatrix(
  data = as.matrix(train_data[, -which(names(train_data) == "class")]),
  label = as.numeric(train_data$class) - 1  # Convert factor levels to binary (0, 1)
)

# Train models
xgb_model <- xgboost(data = xgb_data, objective = "binary:logistic", nrounds = 500, verbose = 0)
glm_model <- glm(class ~ ., data = train_data, family = binomial(link = "logit"))
rf_model <- ranger(class ~ ., data = train_data, probability = TRUE)
gbm_model <- gbm(class ~ ., data = train_data_gbm, distribution = "bernoulli", n.trees = 500)
rpart_model <- rpart(class ~ ., data = train_data, method = "class")

# Initialize ROC data frame
roc_data <- data.frame()

# Evaluate models and compute ROC
models <- list("RF" = rf_model, "GLM" = glm_model, "GBM" = gbm_model, "XGB" = xgb_model, "DT" = rpart_model)

for (model_name in names(models)) {
  model <- models[[model_name]]
  
  # Predict
  if (model_name == "XGB") {
    pred <- predict(model, as.matrix(test_data[, -which(names(test_data) == "class")]))
  } else if (model_name == "RF") {
    pred <- predict(model, data = test_data, type = "response")$predictions[, 2]
  } else if (model_name == "GBM") {
    pred <- predict(model, newdata = test_data_gbm, n.trees = 500, type = "response")
  } else if (model_name == "DT") {
    pred <- predict(model, newdata = test_data, type = "prob")[, 2]
  } else {
    pred <- predict(model, newdata = test_data, type = "response")
  }
  
  # Compute ROC
  roc_res <- roc(response = test_data$class, predictor = pred)
  auc_res <- auc(roc_res)
  
  # Store ROC data
  roc_data <- rbind(roc_data, data.frame(
    model = model_name,
    sensitivity = roc_res$sensitivities,
    specificity = 1 - roc_res$specificities,
    auc = rep(auc_res, length(roc_res$sensitivities))
  ))
}

library(ggplot2)
library(pROC)
library(dplyr)

# Prepare AUC labels dynamically
auc_labels <- roc_data %>%
  group_by(model) %>%
  summarize(
    avg_specificity = mean(specificity),
    avg_sensitivity = mean(sensitivity),
    auc = unique(auc)
  )


FigS3B <- ggplot(data = roc_data, aes(x = specificity, y = sensitivity, color = model)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  # Add dynamic AUC labels at appropriate positions
  geom_text(data = auc_labels, 
            aes(label = paste0(model, " (AUC = ", round(auc, 3), ")"), 
                x = avg_specificity, y = avg_sensitivity), 
            hjust = -0.1, size = 4) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  ggtitle("ROC Curves") +
  scale_color_manual(values = c("#6086ff", "#4cb7a5", "#f29e8e", "#ff0051", "#062ede")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    legend.position = "none"
  )



# get the main isoform for each genes
library(readr)
main_transcript_isoform <- read_csv("data/main_transcript_isoform.csv")
GRCh38_gene_name_id <- read_csv("data/GRCh38_gene_name_id.csv")
load("data/MRL_STM2457_diff.RData")

# Load dplyr for data manipulation
library(dplyr)

# Clean up the `gene_id` column by removing everything after the point
GRCh38_gene_name_id <- GRCh38_gene_name_id %>%
  mutate(gene_id_clean = sub("\\..*", "", gene_id))

# Join the two dataframes based on `transcript_id` and `tracking_id` (removing version number from `tracking_id`)
main_transcript_with_gene_id <- main_transcript_isoform %>%
  mutate(transcript_id_clean = sub("\\..*", "", transcript_id)) %>%
  left_join(GRCh38_gene_name_id %>% 
              mutate(tracking_id_clean = sub("\\..*", "", tracking_id)),
            by = c("transcript_id_clean" = "tracking_id_clean")) %>%
  dplyr::select(transcript_id, UTR5, coding, UTR3, gene_id_clean) # Keep desired columns

# Rename `gene_id_clean` to `gene_id` for clarity
main_transcript_with_gene_id <- main_transcript_with_gene_id %>%
  dplyr::rename(gene_id = gene_id_clean)

MRL_STM2457_diff <- MRL_STM2457_diff %>%
  mutate(gene_id = sub("\\..*", "", gene_id))

main_transcript_with_gene_id <- main_transcript_with_gene_id %>%
  mutate(gene_id_clean = sub("\\..*", "", gene_id))

# Join the dataframes to get the sequences and transcript_id
MRL_STM2457_with_sequences <- merge(MRL_STM2457_diff, main_transcript_with_gene_id, by = "gene_id")


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
           dplyr::select(gene_id_clean, transcript_id, motif_count, sequence_length, ratio, group, part))
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
sequences_upper <- MRL_STM2457_with_sequences %>% filter(highlight_upper1)
sequences_lower <- MRL_STM2457_with_sequences %>% filter(highlight_lower1)


#YTHDF1
# Process both groups
result_upper <- process_group(sequences_upper, YTHDF1, "Upper")
result_lower <- process_group(sequences_lower, YTHDF1, "Lower")

# Combine results for all parts and groups
combined_results <- bind_rows(result_upper, result_lower)
combined_results <- combined_results %>%
  mutate(part = factor(part, levels = c("UTR5", "CDS", "UTR3")))
# Create subfigures for each sequence part with custom colors
plot_df1 <- ggplot(combined_results, aes(x = group, y = ratio, fill = group)) +
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


#YTHDF2

# Process both groups
result_upper <- process_group(sequences_upper, YTHDF2, "Upper")
result_lower <- process_group(sequences_lower, YTHDF2, "Lower")

# Combine results for all parts and groups
combined_results <- bind_rows(result_upper, result_lower)
combined_results <- combined_results %>%
  mutate(part = factor(part, levels = c("UTR5", "CDS", "UTR3")))

# Create subfigures for each sequence part with custom colors
plot_df2 <- ggplot(combined_results, aes(x = group, y = ratio, fill = group)) +
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



# Process both groups
result_upper <- process_group(sequences_upper, YTHDF3, "Upper")
result_lower <- process_group(sequences_lower, YTHDF3, "Lower")

# Combine results for all parts and groups
combined_results <- bind_rows(result_upper, result_lower)
combined_results <- combined_results %>%
  mutate(part = factor(part, levels = c("UTR5", "CDS", "UTR3")))

# Create subfigures for each sequence part with custom colors
plot_df3 <- ggplot(combined_results, aes(x = group, y = ratio, fill = group)) +
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






