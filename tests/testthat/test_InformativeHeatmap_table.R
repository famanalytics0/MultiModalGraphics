library(testthat)
library(MultiModalGraphics)

test_that("AnnotatedHeatmap_table works", {
  # toy FC & p-value matrices
  fc <- matrix(c(-1,2,0,3), nrow=2, dimnames=list(c("g1","g2"), c("M1","M2")))
  pv <- matrix(c(0.01,0.2,0.05,0.001), nrow=2, dimnames=list(rownames(fc), colnames(fc)))
  hm <- AnnotatedHeatmap_table(
    fc_matrix      = fc,
    pval_matrix    = pv,
    pvalue_cutoff  = 0.05,
    trending_cutoff= 0.2,
    fc_cutoff      = 0,
    max_features   = NULL
  )
  expect_s4_class(hm, "AnnotatedHeatmap")
  ht_obj <- getHeatmapObject(hm)
  expect_true(inherits(ht_obj, "Heatmap"))
  # check that rows match features
  expect_equal(rownames(hm@params$fc_matrix), rownames(fc))
})
