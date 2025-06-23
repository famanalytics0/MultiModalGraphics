library(testthat)
library(MultiModalGraphics)

# Test Initialization with Correct Data
test_that("CompositeFeatureHeatmap initializes correctly with valid data", {
  data <- data.frame(
    tissue = factor(rep(c("Tissue1", "Tissue2"), each = 4)),
    signaling = factor(rep(c("Pathway1", "Pathway2", "Pathway3", "Pathway4"), 2)),
    Activation_z_score = runif(8, -2, 2),
    p = runif(8, 0.05),
    number_of_genes = sample(1:100, 8)
  )
  mg <- CompositeFeatureHeatmap(data)
  expect_s4_class(mg, "CompositeFeatureHeatmap")
})

# Test Initialization with Invalid Data
test_that("CompositeFeatureHeatmap fails with invalid data", {
  data <- "not a data frame"
  expect_error(CompositeFeatureHeatmap(data))
})

# Test Plotting Functionality
test_that("plot_heatmap handles parameters correctly", {
  data <- data.frame(
    tissue = factor(rep(c("Tissue1", "Tissue2"), each = 4)),
    signaling = factor(rep(c("Pathway1", "Pathway2", "Pathway3", "Pathway4"), 2)),
    Activation_z_score = runif(8, -2, 2),
    p = runif(8, 0.05),
    number_of_genes = sample(1:100, 8),
    timePoint = rep(c("Time1", "Time2"), 4)
  )
  mg <- CompositeFeatureHeatmap(data)
  expect_silent(plot_heatmap(mg, independantVariable = "timePoint")) # Expect no errors or warnings
})

# Ensure that the function fails gracefully with missing data columns
test_that("plot_heatmap fails with missing data columns", {
  data <- data.frame(
    tissue = factor(rep(c("Tissue1", "Tissue2"), each = 4)),
    signaling = factor(rep(c("Pathway1", "Pathway2", "Pathway3", "Pathway4"), 2)),
    Activation_z_score = runif(8, -2, 2), # Missing p and number_of_genes
    timePoint = rep(c("Time1", "Time2"), 4)
  )
  mg <- CompositeFeatureHeatmap(data)
  expect_error(plot_heatmap(mg, independantVariable = "timePoint"))
})

