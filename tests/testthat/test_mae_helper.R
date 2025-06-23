# tests/testthat/test_mae_helper.R
library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)
library(MultiAssayExperiment)

test_that("MAE helper works end-to-end", {
  expr <- matrix(rnorm(20), nrow=5,
                 dimnames=list(letters[1:5], paste0("S",1:4)))
  meta <- DataFrame(treatment=rep(c("A","B"), each=2),
                    row.names = paste0("S",1:4))
  se  <- SummarizedExperiment(assays = list(counts = expr), colData = meta)
  mae <- MultiAssayExperiment(list(myassay = se))

  hm <- AnnotatedHeatmapFromMAE(mae, assayName="myassay", groupColumn="treatment")
  expect_s4_class(hm, "AnnotatedHeatmap")
})

test_that("MAE helper errors on non-MAE input", {
  expect_error(
    AnnotatedHeatmapFromMAE(42),
    "must be a MultiAssayExperiment"
  )
})
