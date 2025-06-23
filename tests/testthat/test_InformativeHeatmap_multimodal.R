library(testthat)
library(MultiModalGraphics)

test_that("AnnotatedHeatmap(multimodal FC/pval lists) works", {
  fc1 <- matrix(c(1, -2), nrow=2, dimnames=list(c("g1","g2"), "A"))
  pv1 <- matrix(c(0.01, 0.2), nrow=2, dimnames=list(rownames(fc1), "A"))
  fc2 <- matrix(c(0.5, -1), nrow=2, dimnames=list(c("g1","g2"), "B"))
  pv2 <- matrix(c(0.05, 0.3), nrow=2, dimnames=list(rownames(fc2), "B"))

  hm <- AnnotatedHeatmap(
    data      = list(A=fc1, B=fc2),
    pval_list = list(A=pv1, B=pv2)
  )
  expect_s4_class(hm, "AnnotatedHeatmap")
  expect_true(hm@params$multimodal)
})
