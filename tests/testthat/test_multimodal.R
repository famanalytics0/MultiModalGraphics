# tests/testthat/test_multimodal.R
library(testthat)
library(MultiModalGraphics)

test_that("multimodal list path returns AnnotatedHeatmap", {
  fc1 <- matrix(rnorm(12), nrow=3,
                dimnames=list(letters[1:3], paste0("S",1:4)))
  pv1 <- matrix(runif(12), nrow=3,
                dimnames=list(letters[1:3], paste0("S",1:4)))
  fc2 <- matrix(rnorm(12), nrow=3,
                dimnames=list(LETTERS[1:3], paste0("S",1:4)))
  pv2 <- matrix(runif(12), nrow=3,
                dimnames=list(LETTERS[1:3], paste0("S",1:4)))

  hm <- AnnotatedHeatmap(
    data      = list(mod1=fc1, mod2=fc2),
    pval_list = list(mod1=pv1, mod2=pv2)
  )
  expect_s4_class(hm, "AnnotatedHeatmap")
})

test_that("multimodal errors when names mismatch", {
  fc <- matrix(1,1,1)
  pv <- matrix(1,1,1)
  expect_error(
    AnnotatedHeatmap(list(a=fc), list(b=pv)),
    "must match"
  )
})
