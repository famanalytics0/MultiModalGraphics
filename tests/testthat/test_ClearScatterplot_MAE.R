# tests/testthat/test_ClearScatterplot_MAE.R
library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)
library(MultiAssayExperiment)

test_that("ClearScatterplot_MAE works end-to-end", {
  expr <- matrix(rnorm(30), nrow=5,
                 dimnames=list(letters[1:5], paste0("S",1:6)))
  coldata <- DataFrame(
    Group      = rep(c("A","B"), each=3),
    SampleType = rep(c("X","Y"), each=3),
    row.names  = paste0("S",1:6)
  )
  se  <- SummarizedExperiment(assays=list(counts=expr), colData=coldata)
  mae <- MultiAssayExperiment(list(myassay=se))

  cs_mae <- ClearScatterplot_MAE(
    mae,
    assayName   = "myassay",
    groupColumn = "Group",
    sampleType  = "SampleType"
  )
  expect_s4_class(cs_mae, "ClearScatterplot")
})

test_that("ClearScatterplot_MAE errors on missing assay or wrong input", {
  expr <- matrix(1,1,1)
  se   <- SummarizedExperiment(assays=list(counts=expr), colData=DataFrame(Group="A", row.names="1"))
  mae  <- MultiAssayExperiment(list(foo=se))
  expect_error(
    ClearScatterplot_MAE(mae, assayName="bar"),
    "not found in this MAE"
  )
  expect_error(
    ClearScatterplot_MAE(42, assayName="foo"),
    "must be a MultiAssayExperiment"
  )
})
