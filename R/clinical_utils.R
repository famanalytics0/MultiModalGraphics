#' Copy global (clinical) colData to all assays in a MultiAssayExperiment
#' @param mae MultiAssayExperiment
#' @param clinical_cols Character vector of global colData column names to copy (e.g., c("vital_status", "radiation_therapy"))
#' @return MAE with updated colData for each assay
add_clinical_to_all_assays <- function(mae, clinical_cols) {
  top_coldata <- as.data.frame(colData(mae))
  sm <- sampleMap(mae)
  exps <- experiments(mae)
  for (assay_name in names(exps)) {
    se <- exps[[assay_name]]
    sel <- sm$assay == assay_name
    global_samples <- sm$primary[sel]
    assay_samples <- sm$colname[sel]
    # Defensive: subset to available samples and columns
    to_add <- top_coldata[global_samples, clinical_cols, drop = FALSE]
    rownames(to_add) <- assay_samples
    cd <- as.data.frame(colData(se))
    # Copy over columns (overwrite if present)
    for (nm in clinical_cols) {
      cd[[nm]] <- to_add[rownames(cd), nm]
    }
    colData(se) <- S4Vectors::DataFrame(cd)
    exps[[assay_name]] <- se
  }
  experiments(mae) <- exps
  mae
}



# in R/clinical_utils.R

#’ @rdname AnnotatedHeatmap
#’ @export
setMethod(
  "AnnotatedHeatmap",
  signature = c(data = "matrix", meta = "data.frame", pval_list = "ANY"),
  function(data, meta, pval_list = NULL, ...) {
    # look up the original default method for matrix,ANY,ANY
    default <- getMethod("AnnotatedHeatmap", signature = c("matrix","ANY","ANY"))
    default(data, meta, pval_list, ...)
  }
)

#’ @rdname AnnotatedHeatmap
#’ @export
setMethod(
  "AnnotatedHeatmap",
  signature = c(data = "matrix", meta = "data.frame", pval_list = "missing"),
  function(data, meta, ...) {
    # forward to the same default, supplying an explicit NULL
    default <- getMethod("AnnotatedHeatmap", signature = c("matrix","ANY","ANY"))
    default(data, meta, NULL, ...)
  }
)







