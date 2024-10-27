################################################################################
# The script generates synthetic data for multiple assays, creates a MultiAssayExperiment
# object, uses this object as input to iClusterPlus for feature selection, and visualizes
# the results using ComplexHeatmap.
################################################################################

# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("MultiAssayExperiment", "iClusterPlus", "ComplexHeatmap", "magick"))

library(MultiAssayExperiment)
library(iClusterPlus)
library(ComplexHeatmap)
library(magick)
library(MultiModalGraphics)

# Parameters
set.seed(123)  # For reproducibility
num_features <- 1000  # Number of features (e.g., genes, miRNAs, proteins, etc.)
num_samples <- 50  # Number of samples

# User-specified list of assay names (can be 2, 3, 4, 5, or any number)
assay_names <- c("Methylation", "Expression", "miRNA", "Proteomics")  # Example list; replace with user input

# Function to generate synthetic data
generate_synthetic_data <- function(num_features, num_samples) {
  matrix(rnorm(num_features * num_samples, mean = 0, sd = 1), nrow = num_features, ncol = num_samples)
}

# Generate synthetic data for each assay
assay_data_list <- lapply(assay_names, function(assay) {
  data <- generate_synthetic_data(num_features, num_samples)
  rownames(data) <- paste0(assay, 1:num_features)
  colnames(data) <- paste0("Sample", 1:num_samples)
  return(data)
})
names(assay_data_list) <- assay_names

# Sample metadata
sample_metadata <- data.frame(
  SampleID = paste0("Sample", 1:num_samples),
  Group = rep(c("Control", "Treatment"), each = num_samples / 2),
  stringsAsFactors = FALSE
)
rownames(sample_metadata) <- sample_metadata$SampleID

# Create the MultiAssayExperiment object
mae <- MultiAssayExperiment(experiments = assay_data_list, colData = sample_metadata)

# Convert MultiAssayExperiment to a list for iClusterPlus
data_list <- lapply(assays(mae), function(assay) {
  assay_matrix <- as.matrix(assay)
  mode(assay_matrix) <- "numeric"
  return(t(assay_matrix))  # Transpose so that samples are rows and features are columns
})

# Run iClusterPlus for feature selection
set.seed(123)
#fit <- do.call(iClusterPlus, c(data_list, list(type = rep("gaussian", length(assay_names)),
#                                               K = 2, lambda = rep(0.03, length(assay_names)))))
# Create a list of named arguments dt1, dt2, dt3, etc.
dt_list <- setNames(data_list, paste0("dt", seq_along(data_list)))
# Pass the named list along with other required arguments to iClusterPlus
fit <- do.call(iClusterPlus, c(dt_list, list(type = rep("gaussian", length(data_list)),
                                             K = 2, lambda = rep(0.03, length(data_list)))))

# Save the iClusterPlus output
saveRDS(fit, file = "icluster_result.rds")

# Define a threshold for selecting important features based on the magnitude of the coefficients
threshold <- 0.1

# Extract selected features based on the threshold
selected_features_list <- lapply(seq_along(assay_names), function(i) {
  colnames(data_list[[i]])[apply(fit$beta[[i]], 1, function(x) any(abs(x) > threshold))]
})
names(selected_features_list) <- assay_names

# Combine selected features
selected_features <- unique(unlist(selected_features_list))

# Prepare combined data matrix with selected features
combined_selected_data <- do.call(rbind, lapply(seq_along(assay_names), function(i) {
  assay_data_list[[i]][selected_features_list[[i]], , drop = FALSE]
}))

# Ensure the column names match the sample metadata
colnames(combined_selected_data) <- sample_metadata$SampleID

# Transpose combined data matrix for heatmap
heatmap_data <- t(combined_selected_data)

# Ensure the sample_info matches the heatmap data columns
sample_info <- data.frame(Sample = colnames(heatmap_data),
                          Group = sample_metadata$Group[match(colnames(heatmap_data),
                                                              sample_metadata$SampleID)])

# Generate heatmap annotations
ha <- HeatmapAnnotation(Group = sample_info$Group)

# Set options to handle large matrices and use raster
ht_opt$message = FALSE

# Generate the heatmap and save it to a file
png("heatmap.png", width = 2400, height = 1800, res = 300)
heatmap <- Heatmap(heatmap_data, name = "Selected Features", top_annotation = ha,
                   show_row_names = FALSE, show_column_names = TRUE,
                   clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
                   use_raster = TRUE)
draw(heatmap)
dev.off()

# Display the heatmap
img <- image_read("heatmap.png")
print(img)
