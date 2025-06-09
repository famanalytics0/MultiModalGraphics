# tests/testthat/test_MultiModalPlot_multi_panel.R
library(testthat)
library(MultiModalGraphics)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(patchwork)

test_that("multiple volcano panels combine via patchwork", {
  df1 <- data.frame(
    log2fc    = rnorm(10),
    negLog10p = runif(10,1,5),
    regulation= sample(c("up","down"),10,TRUE),
    SampleType= sample(c("A","B"),10,TRUE),
    stringsAsFactors=FALSE
  )
  df2 <- df1
  combo <- MultiModalPlot(
    inputs     = list(M1=df1, M2=df2),
    panel_type = "volcano"
  )
  expect_s3_class(combo, "patchwork")
})

test_that("multiple heatmap panels combine side-by-side", {
  expr <- matrix(rnorm(20), nrow=5,
                 dimnames=list(letters[1:5], paste0("S",1:4)))
  meta <- data.frame(Group=rep(c("A","B"), each=2),
                     SampleType=rep(c("X","Y"),2),
                     row.names=paste0("S",1:4),
                     stringsAsFactors = FALSE)
  ht1 <- list(expr=expr, meta=meta)
  ht2 <- ht1
  combo <- MultiModalPlot(
    inputs     = list(H1=ht1, H2=ht2),
    groupColumns = c(H1="Group", H2="Group"),
    panel_type = "heatmap"
  )
  # this returns a HeatmapList
  expect_true(inherits(combo, "HeatmapList"))
})
