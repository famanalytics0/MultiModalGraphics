#' @rdname AnnotatedHeatmap
#' @export
setMethod(
  "AnnotatedHeatmap",
  signature = c(data = "matrix", meta = "data.frame", pval_list = "ANY"),
  function(data, meta, pval_list = NULL, ...) {
    # delegate to the existing matrix,ANY,ANY method
    callGeneric(data = data, meta = meta, pval_list = pval_list, ...)
  }
)
