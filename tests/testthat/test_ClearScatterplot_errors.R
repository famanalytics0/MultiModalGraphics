library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)
library(MultiAssayExperiment)

test_that("error when missing DE‐data columns", {
  df <- data.frame(a=1:3, b=2:4)
  expect_error(ThresholdedScatterplot(df), "Missing columns")
})

test_that("error on mismatched names in list‐of‐matrices", {
  mat <- matrix(rnorm(4),2)
  meta <- data.frame(Group=1:2, SampleType=c("A","B"), row.names = c("X","Y"))
  expect_error(ThresholdedScatterplot(data = list(A=mat), meta = list(B=meta)),
               "Names of data & meta must match")
})

test_that("error on missing assayNames/groupColumns in list‐of‐MAEs", {
  expr <- matrix(rnorm(6),nrow=3)
  coldata <- DataFrame(Group=c("A","A","B"), SampleType=paste0("S",1:3))
  se <- SummarizedExperiment(assays=list(counts=expr), colData=coldata)
  mae <- MultiAssayExperiment(list(m=se))
  expect_error(ThresholdedScatterplot(data = list(m=mae)),
               "must supply both `assayNames` and `groupColumns`")
})

test_that("error on too few samples", {
  expr <- matrix(rnorm(4),2)
  meta <- data.frame(Group=rep("A",2), SampleType=paste0("S",1:2), row.names=paste0("S",1:2))
  expect_error(ThresholdedScatterplot_table(expr, meta), "Each level of 'Group' must have ≥ 3 samples")
})
