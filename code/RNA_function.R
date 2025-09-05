summary_gene_function <- function(dataframe) {
  # Split the string by ";" or ":"
  split_strings <- strsplit(dataframe$`Protein Function (Protein Atlas)`, ";|:")
  
  # Take the first element of each split string
  first_elements <- sapply(split_strings, function(x) if(length(x) > 0) x[1] else NA)
  
  # Find unique elements and count them
  element_counts <- table(first_elements)
  
  # Convert to a data frame
  result_df <- data.frame("Element" = names(element_counts), "Count" = as.integer(element_counts))
  
  return(result_df)
}


summary_subcellular_location <- function(dataframe) {
  # Split the string by ";" or ":"
  split_strings <- strsplit(dataframe$`Subcellular Location (Protein Atlas)`, ";")
  
  # Take the first element of each split string
  first_elements <- sapply(split_strings, function(x) if(length(x) > 0) x[1] else NA)
  
  # Find unique elements and count them
  element_counts <- table(first_elements)
  
  # Convert to a data frame
  result_df <- data.frame("Element" = names(element_counts), "Count" = as.integer(element_counts))
  
  return(result_df)
}
