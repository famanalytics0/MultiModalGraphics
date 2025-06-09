# tests/testthat/test_MultiModalPlot_mae_volcano.R
library(testthat)
library(MultiModalGraphics)
library(MultiAssayExperiment)
library(SummarizedExperiment)

test_that("MAE → volcano works for a single modality", {
  # create a toy SE and MAE
  expr <- matrix(rnorm(20), nrow=5, dimnames=list(letters[1:5], paste0("S",1:4)))
  coldat <- DataFrame(Group=rep(c("A","B"), each=2),
                      SampleType=rep(c("X","Y"),2),
                      row.names=paste0("S",1:4))
  se <- SummarizedExperiment(list(counts=expr), colData=coldat)
  mae <- MultiAssayExperiment(list(myassay=se))
  p <- MultiModalPlot(
    inputs       = list(ACC=mae),
    assayNames   = c(ACC="counts"),
    groupColumns = c(ACC="Group"),
    sampleTypes  = c(ACC="SampleType"),
    panel_type   = "volcano"
  )
  expect_s3_class(p, "ggplot")
})
