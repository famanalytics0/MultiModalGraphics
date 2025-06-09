# tests/testthat/test_MultiModalPlot_mae_heatmap.R
library(testthat)
library(MultiModalGraphics)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(ComplexHeatmap)

test_that("MAE → heatmap works for a single modality", {
  expr <- matrix(rnorm(20), nrow=5, dimnames=list(letters[1:5], paste0("S",1:4)))
  coldat <- DataFrame(Group=rep(c("A","B"), each=2),
                      SampleType=rep(c("X","Y"),2),
                      row.names=paste0("S",1:4))
  se <- SummarizedExperiment(list(counts=expr), colData=coldat)
  mae <- MultiAssayExperiment(list(myassay=se))
  hm <- MultiModalPlot(
    inputs       = list(ACC=mae),
    assayNames   = c(ACC="counts"),
    groupColumns = c(ACC="Group"),
    sampleTypes  = c(ACC="SampleType"),
    panel_type   = "heatmap",
    col           = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
    cluster_rows  = FALSE
  )
  expect_s4_class(hm, "Heatmap")
})
