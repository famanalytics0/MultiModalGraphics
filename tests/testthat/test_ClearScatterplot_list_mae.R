library(testthat)
library(MultiModalGraphics)
library(SummarizedExperiment)
library(MultiAssayExperiment)

test_that("list‐of‐MAEs constructor works", {
  make_mae <- function(n) {
    expr <- matrix(rnorm(6), nrow=2)
    coldata <- DataFrame(Group=rep(c("A","B"), each= n/2), SampleType=paste0("S",1:n))
    se <- SummarizedExperiment(assays = list(counts = expr), colData = coldata)
    MultiAssayExperiment(list(mod = se))
  }
  mae1 <- make_mae(4); mae2 <- make_mae(4)
  cs <- ClearScatterplot(
    data         = list(First = mae1, Second = mae2),
    assayName    = c(First="mod", Second="mod"),
    groupColumn  = c(First="Group", Second="Group"),
    sampleType   = c(First="SampleType", Second="SampleType")
  )
  expect_s4_class(cs, "ClearScatterplot")
})
