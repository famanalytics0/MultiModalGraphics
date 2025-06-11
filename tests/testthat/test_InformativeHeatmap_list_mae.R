library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)
library(MultiAssayExperiment)

test_that("InformativeHeatmap(list of MAEs) works", {
  # helper to build a toy MAE
  make_mae <- function(n) {
    expr <- matrix(rpois(n*3, 10), nrow=3)
    coldata <- DataFrame(
      Group = rep(c("A","B"), length.out=n),
      row.names=paste0("s",1:n)
    )
    se <- SummarizedExperiment(assays=list(counts=expr), colData=coldata)
    MultiAssayExperiment(list(X=se))
  }
  mae1 <- make_mae(4)
  mae2 <- make_mae(6)
  hm <- InformativeHeatmap(
    data          = list(first=mae1, second=mae2),
    assayNames    = c(first="X", second="X"),
    groupColumns  = c(first="Group", second="Group")
  )
  expect_s4_class(hm, "InformativeHeatmap")
  # multi-MAE should not have single flag
  expect_false(isTRUE(hm@params$single))
})
