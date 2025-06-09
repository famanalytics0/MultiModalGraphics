# tests/testthat/test_MultiModalPlot_exprmeta.R
library(testthat)
library(MultiModalGraphics)

test_that("expr+meta list → volcano and heatmap", {
  # simulate
  expr <- matrix(rnorm(20), nrow=5,
                 dimnames=list(letters[1:5], paste0("S",1:4)))
  meta <- data.frame(Group=rep(c("A","B"), each=2),
                     SampleType=rep(c("X","Y"),2),
                     row.names=paste0("S",1:4),
                     stringsAsFactors = FALSE)
  # volcano
  p <- MultiModalPlot(
    inputs       = list(M1=list(expr=expr, meta=meta)),
    groupColumns = c(M1="Group"),
    panel_type   = "volcano"
  )
  expect_s3_class(p, "ggplot")
  # heatmap
  hm <- MultiModalPlot(
    inputs       = list(M1=list(expr=expr, meta=meta)),
    groupColumns = c(M1="Group"),
    panel_type   = "heatmap"
  )
  expect_s4_class(hm, "Heatmap")
})
