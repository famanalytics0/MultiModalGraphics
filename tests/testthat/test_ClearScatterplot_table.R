# tests/testthat/test_ClearScatterplot_table.R
library(testthat)
library(MultiModalGraphics)

test_that("ClearScatterplot_table works with valid input", {
  expr <- matrix(rpois(20,5), nrow=5)
  rownames(expr) <- letters[1:5]
  colnames(expr) <- paste0("S",1:4)
  meta <- data.frame(
    Group      = rep(c("A","B"), each=2),
    SampleType = c("X","X","Y","Y"),
    row.names  = colnames(expr),
    stringsAsFactors = FALSE
  )
  cs_tab <- ClearScatterplot_table(
    expr,
    meta,
    groupColumn = "Group",
    sampleType  = "SampleType"
  )
  expect_s4_class(cs_tab, "ClearScatterplot")
})

test_that("ClearScatterplot_table errors when too few samples per group", {
  expr <- matrix(rnorm(10), nrow=5)
  rownames(expr) <- letters[1:5]
  colnames(expr) <- paste0("S",1:2)
  meta <- data.frame(
    Group      = c("A","A"),
    SampleType = c("X","Y"),
    row.names  = colnames(expr),
    stringsAsFactors = FALSE
  )
  expect_error(
    ClearScatterplot_table(expr, meta),
    "≥ 3 samples"
  )
})

test_that("ClearScatterplot_table errors when no features remain after filtering", {
  expr <- matrix(rnorm(10), nrow=5)
  rownames(expr) <- letters[1:5]
  colnames(expr) <- paste0("S",1:4)
  meta <- data.frame(
    Group      = rep(c("A","B"), each=2),
    SampleType = c("X","X","Y","Y"),
    row.names  = colnames(expr),
    stringsAsFactors = FALSE
  )
  expect_error(
    ClearScatterplot_table(expr, meta, var_quantile=1),
    "No features remain"
  )
})
