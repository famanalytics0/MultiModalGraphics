# tests/testthat/test_ThresholdedScatterplot_edge_cases.R
library(testthat)
library(MultiModalGraphics)
library(MultiAssayExperiment)
library(SummarizedExperiment)

test_that("core errors when no DE results to plot", {
  # build expr/meta so every cell returns NULL
  expr <- matrix(1, nrow=5, ncol=6,
                 dimnames=list(letters[1:5], paste0("S",1:6)))
  coldata <- DataFrame(
    Group      = rep("A",6),
    SampleType = rep("X",6),
    row.names  = paste0("S",1:6)
  )
  se  <- SummarizedExperiment(assays=list(counts=expr), colData=coldata)
  mae <- MultiAssayExperiment(list(se=se))
  expect_error(
    ThresholdedScatterplot_MAE(mae, "se", groupColumn="Group", sampleType="SampleType"),
    "No DE results to plot"
  )
})

test_that("duplicate sampleType/timePoint combos are handled", {
  # simulate multiple cells with same facet
  expr <- matrix(rnorm(30), nrow=5,
                 dimnames=list(letters[1:5], paste0("S",1:6)))
  meta <- data.frame(
    Group      = rep(c("A","B"), each=3),
    SampleType = rep(c("X","X","Y"), each=2),
    timePoint  = rep(c("T1","T1","T2"), each=2),
    row.names  = paste0("S",1:6),
    stringsAsFactors = FALSE
  )
  cs <- ThresholdedScatterplot_table(
    expr, meta,
    groupColumn = "Group",
    sampleType  = "SampleType",
    timepoint   = "timePoint"
  )
  expect_s4_class(cs, "ThresholdedScatterplot")
  # ensure no error on duplicate combos
})
