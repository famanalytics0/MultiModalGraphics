library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)
library(MultiAssayExperiment)

test_that("MAE constructor works", {
  # build a toy SE
  expr <- matrix(rpois(30, lambda = 10), nrow = 5)
  coldata <- DataFrame(
    Group = rep(c("A","B"), each = 3),
    SampleType = paste0("S", 1:6)
  )
  se <- SummarizedExperiment(assays = list(counts = expr), colData = coldata)
  mae <- MultiAssayExperiment(list(mod = se))
  cs <- ThresholdedScatterplot_MAE(
    mae         = mae,
    assayName   = "mod",
    groupColumn = "Group",
    sampleType  = "SampleType"
  )
  expect_s4_class(cs, "ThresholdedScatterplot")
})
