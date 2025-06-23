library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)
library(MultiAssayExperiment)

test_that("Error: wrong input type", {
  expect_error(AnnotatedHeatmap(123), "`data`.*not recognized")
})

test_that("Error: mismatched names in list-of-matrices", {
  mat1 <- matrix(1:4,2)
  meta1 <- data.frame(Group=1:2, row.names=paste0("s",1:2))
  expect_error(
    AnnotatedHeatmap(data=list(A=mat1), meta=list(B=meta1), groupColumn="Group"),
    "Names of data and meta must match"
  )
})

test_that("Error: missing assayNames/groupColumns for list-of-MAEs", {
  expr <- matrix(1:6, nrow=3)
  coldata <- DataFrame(Group=rep("A",3))
  se <- SummarizedExperiment(assays=list(x=expr), colData=coldata)
  mae <- MultiAssayExperiment(list(X=se))
  expect_error(
    AnnotatedHeatmap(data=list(m=mae)),
    "must supply both `assayNames` and `groupColumns`"
  )
})

test_that("Error: MAE with multiple assays requires assayName", {
  expr <- matrix(1:6,nrow=3)
  se1 <- SummarizedExperiment(assays=list(a=expr), colData=DataFrame(Group=1:3))
  se2 <- SummarizedExperiment(assays=list(b=expr), colData=DataFrame(Group=1:3))
  mae <- MultiAssayExperiment(list(A=se1, B=se2))
  expect_error(
    AnnotatedHeatmapFromMAE(mae),
    "MAE contains multiple assays"
  )
})
