# MultiModalGraphics User Manual

**Version:** 0.99.5

**Authors:** Foziya Ahmed Mohammed, Malick Fall, Rasha Hammieh, Seid Muhie

**License:** MIT

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation & Setup](#installation)
3. [Overview of Package Functionality](#overview)
4. [Section 1: ThresholdedScatterplot (Volcano Plots)](#clear_scatterplot)

   1. [Background & Purpose](#cs_background)
   2. [API Reference](#cs_api)

      * [ThresholdedScatterplot (Table-only)](#cs_table)
      * [ThresholdedScatterplot\_table (Matrix + Metadata)](#cs_matrix)
      * [ThresholdedScatterplot\_MAE (MultiAssayExperiment-based)](#cs_mae)
      * [createPlot](#cs_createplot)
      * [show](#cs_show)
   3. [Practical Examples](#cs_examples)

      * [Example 1A: Built-in `miniACC` (MAE)](#cs_example1a)
      * [Example 1B: `curatedPCaData::Taylor` (Matrix + Metadata)](#cs_example1b)
      * [Example 1C: `curatedTCGAData::BRCA` (Precomputed DE Table)](#cs_example1c)
   4. [Advanced Customization](#cs_advanced)
   5. [Troubleshooting & Tips](#cs_troubleshoot)
5. [Section 2: AnnotatedHeatmap (Custom Heatmaps)](#informative_heatmap)

   1. [Background & Purpose](#ih_background)
   2. [API Reference](#ih_api)

      * [AnnotatedHeatmap (Matrix- or MAE-based DE)](#ih_matrix_mae)
      * [AnnotatedHeatmapFromMAT (Fold-Change + P-Value Matrices)](#ih_matrices)
      * [AnnotatedHeatmapFromMAE (iClusterPlus + DE)](#ih_mae)
      * [updateLayerFun](#ih_update)
      * [getHeatmapObject](#ih_get)
   3. [Practical Examples](#ih_examples)

      * [Example 2A: Subset of `miniACC` for Heatmap (Matrix + Metadata)](#ih_example2a)
      * [Example 2B: Small Simulated Fold-Change + P-Value Matrix](#ih_example2b)
      * [Example 2C: Custom MAE with `AnnotatedHeatmapFromMAE`](#ih_example2c)
   4. [Advanced Customization](#ih_advanced)
   5. [Troubleshooting & Tips](#ih_troubleshoot)
6. [Section 3: CompositeFeatureHeatmap (2D Tile-and-Point Heatmaps)](#multifeature_grid)

   1. [Background & Purpose](#mg_background)
   2. [API Reference](#mg_api)

      * [CompositeFeatureHeatmap (Constructor)](#mg_constructor)
      * [plot\_heatmap](#mg_plot)
   3. [Practical Examples](#mg_examples)

      * [Example 3A: Built-in Demo Data via `get_multifeature_grid_df()`](#mg_example3a)
      * [Example 3B: Custom Table from Pathway Enrichment (Real-World)](#mg_example3b)
   4. [Advanced Customization](#mg_advanced)
   5. [Troubleshooting & Tips](#mg_troubleshoot)
7. [Section 4: Unified Multimodal API (`MultiModalPlot`)](#multimodal_plot)

   1. [Background & Purpose](#mm_background)
   2. [API Reference](#mm_api)
   3. [Practical Examples](#mm_examples)

      * [Example 4A: Single MAE Volcano (Volcano Mode)](#mm_example4a)
      * [Example 4B: Two Modalities, Combined Heatmap (Heatmap Mode)](#mm_example4b)
      * [Example 4C: Three Modalities, Mixed Volcano + Heatmap](#mm_example4c)
   4. [Advanced Customization](#mm_advanced)
   5. [Troubleshooting & Tips](#mm_troubleshoot)
8. [Session Information & Citation](#session)

---

<a name="introduction"></a>

## 1. Introduction

**MultiModalGraphics** is an R/Bioconductor package that streamlines creation of publication-quality visualizations—specifically faceted volcano plots and custom heatmaps—from heterogeneous biological data. It was designed to accommodate real-world workflows that start with:

1. **MultiAssayExperiment (MAE)** objects containing multi-omic assays,
2. **Separate expression/count matrices + metadata** tables,
3. **Precomputed DE result tables**,
4. **Raw fold-change + p-value matrices**,
5. **Feature-level summary tables** for tile-and-point heatmaps.

In addition to class-specific constructors and plot methods, the package offers a **unified API (`MultiModalPlot`)** that automatically dispatches inputs to the right routines and produces combined “multimodal” displays. All functions provide *sensible defaults* yet expose *every threshold, palette, facet, font, or clustering parameter* for power users.

This manual covers three S4 classes:

1. **ThresholdedScatterplot**: Faceted volcano‐plot pipelines.
2. **AnnotatedHeatmap**: Custom ComplexHeatmap pipelines.
3. **CompositeFeatureHeatmap**: 2D tile‐and‐point heatmaps.

Each section provides API documentation, novice‐friendly step-by-step examples (using lightweight real-world datasets), and advanced usage tips.

---

<a name="installation"></a>

## 2. Installation & Setup

### 2.1. Prerequisites

* **R ≥ 4.1.0**
* **Bioconductor ≥ 3.15** (for Bioconductor dependencies and data packages)

### 2.2. Install Bioconductor (if not already)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()  # install/update Bioconductor itself
```

### 2.3. Install MultiModalGraphics & Dependencies

When you install **MultiModalGraphics** from Bioconductor or GitHub, *Imports* are installed automatically. You do not need to install them one-by-one.

#### 2.3.1. From Bioconductor

```r
BiocManager::install("MultiModalGraphics")
```

#### 2.3.2. From GitHub (Development Version)

```r
if (!requireNamespace("git2r", quietly=TRUE)) install.packages("git2r")
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
remotes::install_git(
  "https://github.com/famanalytics0/MultiModalGraphics.git",
  ref         = "master",     # or specific branch/tag
  dependencies= TRUE,
  upgrade     = "never"
)
```

> **Tip:** To verify installed dependencies:
>
> ```r
> tools::package_dependencies("MultiModalGraphics", 
>                             dependencies = c("Imports","LinkingTo"))
> ```

### 2.4. Install Optional Data Packages (Examples)

The following *Suggests* packages are only needed for example code blocks. You can skip them unless you want to run those examples.

```r
BiocManager::install(c("curatedPCaData","curatedTCGAData","curatedTCGADataData"))
# For vignette building / shiny gadget if desired:
install.packages(c("knitr","rmarkdown","testthat","GetoptLong","ggthemes","gridtext","paletteer","reshape","seriation","wesanderson","circlize","shiny"))
```

### 2.5. Load Core Packages

In a fresh R session, load:

```r
library(MultiModalGraphics)
library(ggplot2)            # for custom themes
library(dplyr)              # for data manipulation (group_by, tally)
library(magrittr)           # for `%>%`
library(matrixStats)        # for rowVars (variance filtering)
library(BiocParallel)       # for parallel DE
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(limma)
library(RColorBrewer)       # for color palettes
library(ComplexHeatmap)     # for heatmap functions
```

> **Note:** Load only what you need. For volcano‐only workflows, you need `ggplot2`, `dplyr`, `matrixStats`, `BiocParallel`, `MultiAssayExperiment`, `SummarizedExperiment`, `limma`. Other imports can remain unloaded until used.

---

<a name="overview"></a>

## 3. Overview of Package Functionality

MultiModalGraphics provides three core S4 classes (with associated constructors and plot methods):

1. **ThresholdedScatterplot**

   * Faceted volcano plots of log2 fold-change vs. –log10(p-value).
   * Constructors:

     * `ThresholdedScatterplot(data.frame_of_DE_results, …)`
     * `ThresholdedScatterplot_table(expr_matrix, meta_data, …)`
     * `ThresholdedScatterplot_MAE(mae_object, assayName, …)`
   * Plotting:

     * `createPlot(ThresholdedScatterplot_object, …)`
     * `show(ThresholdedScatterplot_object)`

2. **AnnotatedHeatmap**

   * Customized ComplexHeatmap displays with “layer\_fun” overlays of significant/trending points.
   * Constructors:

     * `AnnotatedHeatmap(expr_matrix, meta_data, …)` (does DE under the hood)
     * `AnnotatedHeatmapFromMAT(logFC_matrix, pval_matrix, …)`
     * `AnnotatedHeatmapFromMAE(mae_object, …)` (optional iClusterPlus + DE).
   * Updates / Retrieval:

     * `updateLayerFun(AnnotatedHeatmap_object, new_layer_fun)`
     * `getHeatmapObject(AnnotatedHeatmap_object)`

3. **CompositeFeatureHeatmap**

   * 2D heatmap + point overlay where tile fill = z‐score, point size/color = –log10(p) & feature count.
   * Constructor:

     * `CompositeFeatureHeatmap(data_frame, …)`
   * Plotting:

     * `plot_heatmap(CompositeFeatureHeatmap_object, …)`

In addition, a **unified high-level wrapper**:

* `MultiModalPlot(inputs_list, assayNames, groupColumns, sampleTypes, timepoints, panel_type = c("volcano","heatmap"), …)`
  automatically dispatches each named input to the correct constructor and stitches the resulting panels side-by-side.

All constructors perform:

* **Input validation**: ensure required columns/assays exist; stop early with clear error if mismatched.
* **NA removal**: drop samples or rows with missing grouping/faceting values.
* **Variance filtering** (optional): drop low-variance features by quantile.
* **DE analysis** (for ThresholdedScatterplot and AnnotatedHeatmap) using `limma` (continuous) or `limma::voom` (counts).
* **Sensible defaults** for thresholds (fold-change, p-value), colors (red/blue/grey), facets, themes.
* **Full customization**: every color, theme, facet formula, clustering parameter, font, text size, etc. is exposed via arguments.
* **Parallel / vectorized** execution for scalability.

Novice users can get a complete faceted volcano or heatmap with two lines of code; advanced users can tune *every* parameter or combine multiple modalities in one figure.

---

<a name="clear_scatterplot"></a>

## 4. Section 1: ThresholdedScatterplot (Volcano Plots)

<a name="cs_background"></a>

### 4.1. Background & Purpose

A **volcano plot** displays log<sub>2</sub> fold-change (x‐axis) versus –log<sub>10</sub>(p‐value) (y‐axis), highlighting significant “up” and “down” features. In multi-omic or multi-condition studies, you may want **faceted volcanoes** per sample type, time point, or other metadata strata. The **ThresholdedScatterplot** class encapsulates:

* **Data slot**: A `data.frame` of DE results for all facets (columns:
  `log2fc`, `negLog10p`, `regulation`, `SampleType`, optional `timePoint`, plus internally computed `color_flag` and `category`).
* **Plot slot**: A `ggplot` object (initially `NULL`) built via `createPlot()`.

Three constructors allow you to start from:

1. **A DE‐table** you already computed (e.g., from DESeq2 or limma).
2. **An expression/count matrix + metadata**; the constructor will run `limma` or `limma::voom + limma` for each “cell” (combination of sampleType/timePoint).
3. **A MultiAssayExperiment (MAE)**; the constructor extracts an assay, merges metadata via `sampleMap` if needed, runs DE per cell, and aggregates results.

After construction, call `createPlot()` to build the faceted `ggplot`. Then call `show()` (or simply print the object) to render it.

---

<a name="cs_api"></a>

### 4.2. API Reference

#### 4.2.1. ThresholdedScatterplot (Table-only Constructor) <a name="cs_table"></a>

```r
ThresholdedScatterplot(
  data,                  # data.frame with DE results
  highLog2fc     =  0.585,   # |log2FC| > 0.585 (≃1.5×) → "up"/"down"
  lowLog2fc      = -0.585,
  negLog10pValue =  1.301,   # –log10(p) > 1.301 → p < 0.05
  dropNA         = TRUE,     # drop rows with NA in key columns
  …
)
```

* **Required columns in** `data`:

  * `log2fc` (numeric)
  * `negLog10p` (numeric = –log<sub>10</sub>(p‐value))
  * `regulation` (character: `"up"` or `"down"`)
  * `SampleType` (character or factor)
  * Optionally, `timePoint` (character or factor)

* **Checks & Sanity**:

  1. `stopifnot(is.data.frame(data))`.
  2. Ensure `c("log2fc","negLog10p","regulation","SampleType") %in% colnames(data)`.
  3. If `dropNA=TRUE`: drop rows with `NA` in `log2fc` or `negLog10p` (with a warning). Stop if zero rows remain.
  4. Confirm `data$log2fc` and `data$negLog10p` are numeric; else `stop()`.
  5. Compute `color_flag`:

     ```r
     data$color_flag <- with(data,
        ifelse(log2fc >  highLog2fc  & negLog10p > negLog10pValue,  1,
        ifelse(log2fc <  lowLog2fc   & negLog10p > negLog10pValue, -1, 0)))
     data$category <- factor(data$color_flag,
        levels = c(-1,0,1),
        labels = c("down","neutral","up"))
     ```
  6. Return:

     ```r
     methods::new("ThresholdedScatterplot", data = data, plot = NULL)
     ```

* **Exports**:

  ```r
  @exportClass ThresholdedScatterplot
  @exportMethod ThresholdedScatterplot
  ```

---

#### 4.2.2. ThresholdedScatterplot\_table (Expression Matrix + Metadata) <a name="cs_matrix"></a>

```r
ThresholdedScatterplot_table(
  expr,                   # numeric matrix: features × samples
  meta,                   # data.frame: samples × metadata; rownames(meta) = colnames(expr)
  groupColumn    = "Group",     # metadata column for DE contrast
  sampleType     = "SampleType",# metadata column for X‐facet
  timepoint      = NULL,        # metadata column for Y‐facet (optional)
  dataType       = c("auto","continuous","count"),
  vectorized     = c("auto","perCell","vectorized"),
  parallel       = TRUE,
  BPPARAM        = BiocParallel::bpparam(),
  var_quantile   = 0.75,        # drop features below 75th percentile of variance
  pvalue_cutoff  = 0.05,
  fc_cutoff      = 0.585,
  min_samples    = 3,           # at least min_samples per group in each cell
  max_features   = NULL,        # optionally cap to top n features by variance
  dropNA         = TRUE,        # drop samples with NA in grouping/faceting columns
  …
)
```

* **Checks & Sanity**:

  1. `stopifnot(is.matrix(expr), is.data.frame(meta))`.
  2. Ensure `colnames(expr)` and `rownames(meta)` match exactly (case-sensitive). If not, try reordering `meta <- meta[colnames(expr), , drop=FALSE]`; else `stop()`.
  3. Verify `groupColumn %in% colnames(meta)` and `sampleType %in% colnames(meta)`. If `timepoint` non-NULL, also check `timepoint %in% colnames(meta)`. Else, `stop()`.
  4. If `dropNA=TRUE`:

     ```r
     needed_cols <- c(groupColumn, sampleType, timepoint)[!is.null(c(groupColumn, sampleType, timepoint))]
     keep_idx <- rowSums(is.na(meta[, needed_cols, drop=FALSE])) == 0
     if (!all(keep_idx)) {
       warning(sum(!keep_idx), " sample(s) removed due to NA in group/facet cols.")
       expr <- expr[, keep_idx, drop=FALSE]
       meta <- meta[keep_idx, , drop=FALSE]
     }
     ```
  5. **Variance filter** (if `var_quantile ∈ (0,1)`):

     ```r
     rv <- matrixStats::rowVars(expr, na.rm=TRUE)
     thr <- quantile(rv, var_quantile, na.rm=TRUE)
     keep_genes <- rv >= thr
     if (sum(keep_genes)==0) stop("No features pass the variance filter.")
     expr <- expr[keep_genes, , drop=FALSE]
     meta  <- meta                 # metadata unaffected
     ```
  6. If `max_features` is numeric < `nrow(expr)`, keep top `max_features` by variance.
  7. **Determine `dataType`**:

     ```r
     dataType <- match.arg(dataType)
     if (dataType == "auto") {
       dataType <- if (all(expr == floor(expr), na.rm=TRUE) && max(expr, na.rm=TRUE) > 30)
                     "count" else "continuous"
     }
     ```
  8. **Group size checks**: For each combination of `(SampleType, timePoint)`, ensure at least `min_samples` per level of `groupColumn`.
  9. **Split into “cells”**: Unique combinations of `sampleType × timePoint` (or just `sampleType` if `timepoint=NULL`).
  10. **Per-cell DE** (vectorized or parallel via `bplapply`):

      * Subset `expr[:, idx]` and `meta[idx, ]` for that cell.
      * If `dataType=="count"`, run `voom()`, else if `"continuous"` and `max(ce)>50`, log2‐transform.
      * `lmFit()`, `eBayes()`, `topTable()` → `data.frame(log2fc, negLog10p, regulation, SampleType, timePoint)`.
      * If zero results, warning and skip.
      * Combine all cell results with `do.call(rbind, lst)`.
      * If combined DE table is empty, `stop("No DE results to plot.")`.
  11. **Call** `ThresholdedScatterplot(de_table, highLog2fc=fc_cutoff, lowLog2fc=-fc_cutoff, negLog10pValue=-log10(pvalue_cutoff), dropNA=TRUE)`; return that object.

* **Exports**:

  ```r
  @exportClass ThresholdedScatterplot_table
  @exportMethod ThresholdedScatterplot_table
  ```

---

#### 4.2.3. ThresholdedScatterplot\_MAE (MultiAssayExperiment-based Constructor) <a name="cs_mae"></a>

```r
ThresholdedScatterplot_MAE(
  mae,                    # MultiAssayExperiment object
  assayName,              # character: name of assay to extract
  groupColumn     = "Group",     # metadata column (in assay-level or top-level colData) for DE
  sampleType      = NULL,        # metadata column for X-facet
  timepoint       = NULL,        # metadata column for Y-facet
  dataType        = c("auto","continuous","count"),
  vectorized      = c("auto","perCell","vectorized"),
  parallel        = TRUE,
  BPPARAM         = BiocParallel::bpparam(),
  var_quantile    = 0.75,
  pvalue_cutoff   = 0.05,
  fc_cutoff       = 0.585,
  min_samples     = 3,
  max_features    = NULL,
  dropNA          = TRUE,
  …
)
```

* **Checks & Sanity**:

  1. `stopifnot(inherits(mae, "MultiAssayExperiment"))`.
  2. `if (!(assayName %in% names(experiments(mae)))) stop("Assay '", assayName, "' not found in MAE.")`.
  3. Extract assay:

     ```r
     se <- experiments(mae)[[assayName]]
     expr_full <- SummarizedExperiment::assay(se)
     stopifnot(is.matrix(expr_full))
     mode(expr_full) <- "numeric"
     ```
  4. **Metadata merging**:

     * `meta_assay <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors=FALSE)`
     * If any of `groupColumn, sampleType, timepoint` not in `colnames(meta_assay)`, attempt to merge with top-level `colData(mae)` via `sampleMap(mae)`.
     * After merging, ensure `rownames(meta)` correspond to sample IDs that match `colnames(expr_full)`.
  5. **Intersect samples**:

     ```r
     shared <- intersect(colnames(expr_full), rownames(meta))
     if (length(shared)==0) stop("No overlapping samples between assay and metadata.")
     expr <- expr_full[, shared, drop=FALSE]
     meta <- meta[shared, , drop=FALSE]
     ```
  6. **Drop NA** in `groupColumn, sampleType, timepoint` (if `dropNA=TRUE`), with warning.
  7. **Group size checks**: Ensure each level of `groupColumn` ≥ `min_samples`; total samples ≥ `2*min_samples`.
  8. **Variance filter** (if `var_quantile ∈ (0,1)`): Same as in table version.
  9. **Cap features** (`max_features`) by variance if specified.
  10. **Determine `dataType`**: as in table version.
  11. **Split into cells**: combinations of `(sampleType, timepoint)` (or just `sampleType` if `timepoint=NULL`).
  12. **Per-cell DE**: identical logic as in `ThresholdedScatterplot_table`, using `voom()` for `count` or direct `lmFit()` for `continuous`. Use `bplapply()` if `(vectorized=="vectorized" || (vectorized=="auto" && n.cells > bpworkers(BPPARAM))) && parallel==TRUE`.
  13. Combine DE tables; if empty, `stop("No DE results to plot.")`.
  14. **Call** `ThresholdedScatterplot(de_table, highLog2fc=fc_cutoff, lowLog2fc=-fc_cutoff, negLog10pValue=-log10(pvalue_cutoff), dropNA=TRUE)`; return that object.

* **Exports**:

  ```r
  @exportClass ThresholdedScatterplot_MAE
  @exportMethod ThresholdedScatterplot_MAE
  ```

---

#### 4.2.4. createPlot (Faceted Volcano Builder) <a name="cs_createplot"></a>

```r
setMethod(
  "createPlot",
  signature(object = "ThresholdedScatterplot"),
  function(
    object,
    color_up          = "indianred",         # color for “up” (category == "up")
    color_neutral     = "grey",              # color for “neutral” (category == "neutral")
    color_down        = "cornflowerblue",    # color for “down” (category == "down")
    xlab              = expression(log2~fold~change),
    ylab              = expression(-log10~p),
    legend_position   = "bottom",            # "bottom","top","left","right"
    legend_title      = NULL,                # title of color legend
    legend_labels     = NULL,                # c("down","neutral","up") or custom
    text_family       = "sans",              # font family
    text_size         = 10,                  # base font size
    point_size        = 1.75,                # point size
    point_alpha       = 0.5,                 # point transparency
    custom_theme      = NULL,
    …
  ) {
    df <- object@data
    has_tp <- "timePoint" %in% names(df) && !all(is.na(df$timePoint))
    facet <- if (has_tp && length(unique(df$timePoint)) > 1) {
               "timePoint ~ SampleType"
             } else {
               ". ~ SampleType"
             }

    # Base ggplot
    p <- ggplot2::ggplot(
           df,
           ggplot2::aes(x = log2fc, y = negLog10p, color = factor(color_flag))
         ) +
         ggplot2::geom_point(alpha = point_alpha, size = point_size) +
         ggplot2::labs(
           x = xlab,
           y = ylab,
           color = legend_title,
           …
         ) +
         ggplot2::scale_color_manual(
           values = c(color_down, color_neutral, color_up),
           labels = legend_labels
         ) +
         ggplot2::facet_grid(stats::as.formula(facet), space = "free") +
         ggplot2::theme_bw() +
         ggplot2::theme(
           text            = ggplot2::element_text(family = text_family, size = text_size),
           strip.background= ggplot2::element_rect(fill = "white", color = "black"),
           strip.text      = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
           axis.title      = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
           axis.text       = ggplot2::element_text(size = text_size, family = text_family),
           legend.position = legend_position,
           legend.title    = ggplot2::element_text(size = text_size + 2, face = "bold", family = text_family),
           legend.text     = ggplot2::element_text(size = text_size, family = text_family),
           panel.grid.major= ggplot2::element_line(color = "grey80"),
           panel.grid.minor= ggplot2::element_blank()
         )

    if (!is.null(custom_theme)) {
      p <- p + custom_theme
    }

    # Compute counts per facet
    if (has_tp && length(unique(df$timePoint)) > 1) {
      up_df <- df %>%
        dplyr::filter(color_flag == 1) %>%
        dplyr::group_by(timePoint, SampleType) %>%
        dplyr::tally(name = "n_up")
      dn_df <- df %>%
        dplyr::filter(color_flag == -1) %>%
        dplyr::group_by(timePoint, SampleType) %>%
        dplyr::tally(name = "n_down")
    } else {
      up_df <- df %>%
        dplyr::filter(color_flag == 1) %>%
        dplyr::group_by(SampleType) %>%
        dplyr::tally(name = "n_up")
      dn_df <- df %>%
        dplyr::filter(color_flag == -1) %>%
        dplyr::group_by(SampleType) %>%
        dplyr::tally(name = "n_down")
    }

    # Add up/down counts
    if (nrow(up_df) > 0) {
      p <- p + ggplot2::geom_text(
             data   = up_df,
             ggplot2::aes(label = n_up),
             x      = Inf, y = Inf,
             hjust  = 1.1, vjust = 1.1,
             color  = color_up,
             family = text_family,
             size   = text_size / ggplot2::.pt
           )
    }
    if (nrow(dn_df) > 0) {
      p <- p + ggplot2::geom_text(
             data   = dn_df,
             ggplot2::aes(label = n_down),
             x      = -Inf, y = Inf,
             hjust  = -0.1, vjust = 1.1,
             color  = color_down,
             family = text_family,
             size   = text_size / ggplot2::.pt
           )
    }

    object@plot <- p
    invisible(object)
  }
)
```

* **Parameters**:

  * `color_up`, `color_neutral`, `color_down`: Colors for “up”, “neutral”, “down.”
  * `xlab`, `ylab`: Axis labels (strings or expressions).
  * `legend_position`: `"bottom"`, `"top"`, `"left"`, or `"right"`.
  * `legend_title`: Legend title for color keys.
  * `legend_labels`: Character vector of length 3 for `c("down","neutral","up")`.
  * `text_family`, `text_size`: Font family & base size.
  * `point_size`, `point_alpha`: Dot size & transparency.
  * `custom_theme`: Extra `ggplot2` theme (e.g., `theme_minimal()`).
  * `…`: Additional `ggplot2::labs()` arguments (e.g., `title`, `subtitle`, `caption`).

* **Exports**:

  ```r
  @export createPlot
  ```

---

#### 4.2.5. show (Render Volcano Plot) <a name="cs_show"></a>

```r
setMethod(
  "show",
  signature(object = "ThresholdedScatterplot"),
  function(object) {
    if (is.null(object@plot)) {
      stop("Plot not yet built. Run `object <- createPlot(object, …)` first.")
    }
    print(object@plot)
    invisible(NULL)
  }
)
```

* **Purpose**: Print the internal `ggplot` stored in `object@plot`.
* **Exports**:

  ```r
  @exportMethods(show)
  ```

---

<a name="cs_examples"></a>

### 4.3. Practical Examples

Below are three lightweight, **real-world-style** examples—one for each constructor. Each example can be run in a fresh R session; we assume required packages are installed. We deliberately choose small subsets (or built-in “mini” datasets) to keep runtime & memory low.

---

<a name="cs_example1a"></a>

#### Example 1A: Built-in `miniACC` (MAE)

**Goal:** Faceted volcano of `miniACC`, DE between `"C1A"` vs `"C1B"`, stratified by `"pathologic_stage"` (columns) and `"MethyLevel"` (rows).

1. **Load Packages & Data**

   ```r
   library(MultiModalGraphics)
   library(MultiAssayExperiment)
   library(SummarizedExperiment)
   library(matrixStats)
   library(BiocParallel)
   library(limma)
   library(dplyr)
   library(ggplot2)

   data("miniACC", package = "MultiAssayExperiment")
   miniACC
   ```

   *Expected output:*

   ```
   class: MultiAssayExperiment 
   dim: 21683 79 
   experiments: RNASeq2GeneNorm(1) 
   colData: 79 samples × ~30 columns
   sampleMap: 79 rows
   ```

2. **Define Parallel Backend (Optional)**

   Use 2 cores for speed:

   ```r
   BPPARAM <- MulticoreParam(workers = 2)
   ```

3. **Construct `ThresholdedScatterplot_MAE` Object**

   ```r
   cs_mae <- ThresholdedScatterplot_MAE(
     mae         = miniACC,
     assayName   = "RNASeq2GeneNorm",
     groupColumn = "C1A.C1B",
     sampleType  = "pathologic_stage",
     timepoint   = "MethyLevel",
     dataType    = "auto",       
     vectorized  = "auto",       
     parallel    = TRUE,         
     BPPARAM     = BPPARAM,
     var_quantile= 0.75,         
     pvalue_cutoff = 0.05,       
     fc_cutoff     = 0.585,      
     min_samples   = 3,
     max_features  = 5000        # keep top 5k by variance; speeds things up
   )
   ```

   * Internally:

     1. Extract assay matrix (`expr`) and metadata (`meta`).
     2. Drop low-variance genes (bottom 75%).
     3. Align samples, drop any with `NA` in `"C1A.C1B"`, `"pathologic_stage"`, `"MethyLevel"`.
     4. For each `(pathologic_stage, MethyLevel)` cell, run `voom()` → `limma::lmFit()` → `topTable()`.
     5. Combine results into `cs_mae@data`.

4. **Build & Render Volcano Plot**

   ```r
   cs_mae <- createPlot(
     cs_mae,
     color_up        = "tomato",        # “up” = tomato
     color_neutral   = "lightgrey",     # “neutral” = light grey
     color_down      = "steelblue",     # “down” = steelblue
     xlab            = "Log2 Fold Change",
     ylab            = "-Log10 P-Value",
     text_family     = "Times New Roman",
     text_size       = 12,
     point_size      = 2,
     point_alpha     = 0.7,
     legend_position = "right",
     legend_title    = "Regulation",
     legend_labels   = c("Up", "Neutral", "Down"),
     custom_theme    = theme_minimal(),
     title           = "miniACC: C1A vs C1B Volcano",
     subtitle        = "Faceted by Pathologic Stage (X) & Methylation Level (Y)",
     caption         = "Data: miniACC (MultiAssayExperiment)"
   )

   show(cs_mae)
   ```

   **What you see:**

   * A multi-panel volcano:

     * Columns = “Stage I”, “Stage II”, … (pathologic\_stage)
     * Rows = “high”, “low” (MethyLevel)
   * Points colored by `category`:

     * “Up” (tomato), “Down” (steelblue), “Neutral” (lightgrey)
   * Top-corner counts of up/down per facet.
   * Title, subtitle, caption as specified.

---

<a name="cs_example1b"></a>

#### Example 1B: `curatedPCaData::Taylor` (Matrix + Metadata)

**Goal:** Faceted volcano of prostate cancer “Taylor” dataset: DE Tumor vs Normal, faceted by Gleason Score (X) and Race (Y).

> **Note:** This example requires installing `curatedPCaData`. The “Taylor” dataset has \~12,600 genes × 106 samples. We filter to top 25% variance to speed up.

1. **Load Packages & Data**

   ```r
   BiocManager::install("curatedPCaData")    # if not already installed
   library(MultiModalGraphics)
   library(curatedPCaData)
   library(SummarizedExperiment)
   library(matrixStats)
   library(BiocParallel)
   library(limma)
   library(dplyr)
   library(ggplot2)

   # Fetch “Taylor” prostate cancer data
   pcad_list <- getPCa("Taylor")
   se_taylor <- pcad_list$Taylor
   ```

   *Quick check:*

   ```r
   se_taylor
   # SummarizedExperiment: 12625 genes × 106 samples
   # assay: counts; colData has “DiseaseStatus”, “GleasonScore”, “Race”, etc.
   ```

2. **Extract Expression Matrix & Metadata**

   ```r
   expr_mat <- assay(se_taylor, "counts")                  # 12625 × 106
   meta_df   <- as.data.frame(colData(se_taylor), stringsAsFactors = FALSE)  # 106 × ~8

   # Ensure alignment (rownames(meta_df) must match colnames(expr_mat)):
   if (!all(colnames(expr_mat) == rownames(meta_df))) {
     meta_df <- meta_df[colnames(expr_mat), , drop=FALSE]
   }
   ```

3. **Define Parallel Backend (Optional)**

   ```r
   BPPARAM <- MulticoreParam(workers = 2)
   ```

4. **Construct `ThresholdedScatterplot_table` Object**

   ```r
   cs_table <- ThresholdedScatterplot_table(
     expr          = expr_mat,
     meta          = meta_df,
     groupColumn   = "DiseaseStatus",   # “Tumor” vs “Normal”
     sampleType    = "GleasonScore",    # e.g. “6”, “7”, “8+”
     timepoint     = "Race",            # e.g. “white”, “black”, “other”
     dataType      = "auto",            # auto-detect “count” → use voom
     vectorized    = "auto",
     parallel      = TRUE,
     BPPARAM       = BPPARAM,
     var_quantile  = 0.75,              # keep top 25% by variance
     pvalue_cutoff = 0.05,
     fc_cutoff     = 1,                 # |log2FC| > 1
     min_samples   = 3,
     max_features  = 3000               # keep top 3k genes by variance
   )
   ```

5. **Build & Render Volcano Plot**

   ```r
   cs_table <- createPlot(
     cs_table,
     color_up        = "darkgreen",
     color_neutral   = "lightgrey",
     color_down      = "maroon",
     xlab            = "Log2 Fold Change (Tumor vs Normal)",
     ylab            = "-Log10 P-Value",
     text_family     = "Helvetica",
     text_size       = 11,
     point_size      = 2,
     point_alpha     = 0.6,
     legend_position = "bottom",
     legend_title    = "Regulation",
     legend_labels   = c("Up", "Neutral", "Down"),
     custom_theme    = theme_classic(),
     title           = "Prostate Cancer (Taylor): Tumor vs Normal",
     subtitle        = "Faceted by Gleason Score (X) & Race (Y)",
     caption         = "Data: Taylor (curatedPCaData)"
   )

   show(cs_table)
   ```

   **What you see:**

   * Columns = Gleason Score categories (“6”, “7”, “8+”).
   * Rows = Race categories (“white”, “black”, “other”).
   * Points colored:

     * “Up” (darkgreen), “Down” (maroon), “Neutral” (lightgrey).
   * Counts per facet in corners.
   * Title, subtitle, caption per specification.

---

<a name="cs_example1c"></a>

#### Example 1C: `curatedTCGAData::BRCA` (Precomputed DE Table)

**Goal:** Run DE manually on two subtypes (“Basal-like” vs “Luminal A”) of TCGA-BRCA, then plot with `ThresholdedScatterplot(data.frame_DE)`.

> **Note:** We only demonstrate building a small DE table and passing it to `ThresholdedScatterplot`. We do not facet by `timePoint` (set to `NA`).

1. **Load Packages & Data**

   ```r
   BiocManager::install("curatedTCGAData")    # if not installed
   library(MultiModalGraphics)
   library(curatedTCGAData)
   library(SummarizedExperiment)
   library(limma)
   library(matrixStats)
   library(ggplot2)
   library(dplyr)

   # Get TCGA-BRCA (use default version or specify one)
   brca_list <- curatedTCGAData("BRCA", dry.run = FALSE)
   se_brca <- experiments(brca_list$BRCA_RNASeq2GeneNorm)[[1]]  
   # or brca_list[[1]] if only one element
   ```

   *Check colData:*

   ```r
   meta_brca <- as.data.frame(colData(se_brca), stringsAsFactors=FALSE)
   table(meta_brca$PAM50_subtype_PAM50)
   # “Basal-like”: ~189; “Luminal A”: ~565; others & NA
   ```

2. **Subset to “Basal-like” vs “Luminal A”**

   ```r
   keep <- meta_brca$PAM50_subtype_PAM50 %in% c("Basal-like","Luminal A")
   expr2 <- assay(se_brca, "RNASeq2GeneNorm")[, keep]
   meta2 <- meta_brca[keep, , drop=FALSE]

   # Reorder metadata to match columns of expr2
   if (!all(colnames(expr2) == rownames(meta2))) {
     meta2 <- meta2[colnames(expr2), , drop=FALSE]
   }
   ```

3. **Run DE with limma (`continuous` workflow)**

   ```r
   group_factor <- factor(meta2$PAM50_subtype_PAM50, levels = c("Luminal A","Basal-like"))
   design       <- model.matrix(~ group_factor)

   fit  <- lmFit(expr2, design)
   fit  <- eBayes(fit)
   tt   <- topTable(fit, coef = "group_factorBasal-like", number = Inf)

   df_de <- data.frame(
     log2fc     = tt$logFC,
     negLog10p  = -log10(tt$P.Value),
     regulation = ifelse(tt$logFC > 0, "up", "down"),
     SampleType = "Basal-vs-LuminalA",
     timePoint  = NA_character_,
     stringsAsFactors = FALSE,
     row.names  = rownames(tt)
   )

   # Drop any NAs or infinite
   bad <- is.na(df_de$log2fc) | is.na(df_de$negLog10p) | is.infinite(df_de$negLog10p)
   if (any(bad)) {
     warning(sum(bad), " row(s) removed due to NA/Inf in log2fc or negLog10p.")
     df_de <- df_de[!bad, ]
   }
   ```

4. **Construct & Plot Volcano with ThresholdedScatterplot**

   ```r
   cs_de <- ThresholdedScatterplot(
     data            = df_de,
     highLog2fc      = 1.0,   # |log2FC| > 1
     lowLog2fc       = -1.0,
     negLog10pValue  = 1.301  # –log10(p) > 1.301 (p < 0.05)
   )

   cs_de <- createPlot(
     cs_de,
     color_up        = "purple",
     color_neutral   = "grey",
     color_down      = "darkred",
     xlab            = "Log2 Fold Change (Basal-like vs Luminal A)",
     ylab            = "-Log10 P-Value",
     text_family     = "Arial",
     text_size       = 12,
     point_size      = 2,
     point_alpha     = 0.6,
     legend_position = "top",
     legend_title    = "Regulation",
     legend_labels   = c("Down","Neutral","Up"),
     custom_theme    = theme_light(),
     title           = "TCGA-BRCA: Basal-like vs Luminal A Volcano",
     subtitle        = "limma DE (continuous data)",
     caption         = "Data: TCGA BRCA (curatedTCGAData)"
   )

   show(cs_de)
   ```

   **What you see:** A single volcano plot (no facets) colored by regulation. Top-corner labels show counts of up/down genes.

---

<a name="cs_advanced"></a>

### 4.4. Advanced Customization

1. **Axis Expressions**
   Use R expressions for scientific notation:

   ```r
   createPlot(cs_mae,
              xlab = expression(log[2]~fold~change),
              ylab = expression(-log[10]~p),
              …)
   ```

2. **Custom Color Palettes**
   Supply named color vectors:

   ```r
   createPlot(cs_table,
              color_up      = "#1B9E77",   # teal 
              color_neutral = "#D3D3D3",   # light grey
              color_down    = "#D95F02",   # orange 
              …)
   ```

3. **Changing Legend Labels/Order**

   ```r
   createPlot(cs_de,
              legend_labels   = c("Down-regulated","No Change","Up-regulated"),
              color_down      = "navy",
              color_neutral   = "grey50",
              color_up        = "firebrick",
              …)
   ```

4. **Themes**
   Append any `ggplot2` theme:

   ```r
   createPlot(cs_mae, custom_theme = theme_classic() + theme(plot.background = element_rect(fill="lightyellow")))
   ```

5. **Faceting Without `timePoint`**
   If your data has no `timePoint` column (all `NA`), the code falls back to `. ~ SampleType` (single-row layout).

6. **Plotting in Two Steps**
   Store then print:

   ```r
   cs <- ThresholdedScatterplot(df_de, highLog2fc=1, lowLog2fc=-1, negLog10pValue=1.301)
   cs <- createPlot(cs, …)
   # Do other data prep, then:
   show(cs)
   ```

---

<a name="cs_troubleshoot"></a>

### 4.5. Troubleshooting & Tips

1. **“Column not found” Errors**

   * *Symptom*:

     ```
     Error in ThresholdedScatterplot_table(...): Missing columns: X, Y
     ```
   * *Fix*: Check `colnames(meta)` and ensure exact, case-sensitive match for `groupColumn`, `sampleType`, and `timepoint`.

2. **“No overlapping samples”**

   * *Symptom*:

     ```
     Error in ThresholdedScatterplot_MAE(...): No overlapping samples between expression and metadata
     ```
   * *Fix*: Verify that `colnames(expr_matrix)` match `rownames(meta)` exactly. If using an MAE, confirm `sampleMap(mae)` is correct or merge top-level `colData(mae)` properly.

3. **“<3 samples per group”**

   * *Symptom*:

     ```
     Error in ThresholdedScatterplot_table(...): Each level of 'DiseaseStatus' must have ≥3 samples
     ```
   * *Fix*: Check `table(meta$DiseaseStatus, meta$SampleType, meta$timePoint)`. Choose different grouping/faceting columns or adjust `min_samples` argument.

4. **Zero DE Results**

   * *Symptom*:

     ```
     Warning: No DE results for SampleType=“X” / timePoint=“Y”
     Error: No DE results to plot.
     ```
   * *Cause*: All p-values > 1 (→ no –log10(p) > 1.301) or fewer than 2 samples in a cell.
   * *Fix*: Lower thresholds (`pvalue_cutoff`, `fc_cutoff`) or remove one facet (set `timepoint=NULL`) to reduce fragmentation.

5. **`createPlot` Not Called Before `show`**

   * *Symptom*:

     ```
     Error in show(cs): Plot not yet built. Run `object <- createPlot(object, …)` first.
     ```
   * *Fix*: Always run `cs <- createPlot(cs, …)` before `show(cs)`.

6. **Missing ggplot2/dplyr Functions**

   * *Symptom*:

     ```
     Error: could not find function "ggplot"
     ```
   * *Fix*: `library(ggplot2); library(dplyr); library(magrittr)`

7. **Performance Tips**

   * Increase `BPPARAM <- MulticoreParam(workers = n)` for more cores.
   * Lower `var_quantile` (e.g., 0.5) to filter more genes early.
   * Use `max_features` to cap DE on only top N features by variance.

---

<a name="informative_heatmap"></a>

## 5. Section 2: AnnotatedHeatmap (Custom Heatmaps)

<a name="ih_background"></a>

### 5.1. Background & Purpose

A **ComplexHeatmap** is a highly customizable heatmap that can integrate multiple annotations and “layer functions” (custom drawing of grid elements). The **AnnotatedHeatmap** class wraps around a `ComplexHeatmap::Heatmap` object, adding:

* **Built-in differential expression pipelines** (for matrix + metadata or MAE) to compute fold-changes and p-values.
* **Optional clustering via iClusterPlus** (for multi-omic MAE), combining multiple assay DE results.
* **Custom `layer_fun`** overlay that plots points on each cell (tile) whose color indicates statistical significance (p < threshold) or “trending” (p between two thresholds).

Core slots:

* `@heatmap`: A `Heatmap` object (from ComplexHeatmap).
* `@params`: A `list` of parameters used in constructor (e.g., thresholds, palette, etc.), to allow updating.

Three entry points:

1. **`AnnotatedHeatmap(expr_matrix, meta_data, …)`** — runs DE per facet and builds an expression/fold-change matrix.
2. **`AnnotatedHeatmapFromMAT(logFC_matrix, pvalue_matrix, …)`** — uses precomputed fold-change + p-value matrices.
3. **`AnnotatedHeatmapFromMAE(mae_object, …)`** — runs iClusterPlus (optional) → DE on clusters in each assay → assembles combined fold-change & p-value matrices → builds heatmap.

Additional methods:

* `updateLayerFun(AnnotatedHeatmap_object, new_layer_fun)`: Replace the overlay function, re-draw heatmap.
* `getHeatmapObject(AnnotatedHeatmap_object)`: Extract the underlying `Heatmap` (ComplexHeatmap) for further manipulation.

---

<a name="ih_api"></a>

### 5.2. API Reference

#### 5.2.1. AnnotatedHeatmap (Matrix + Metadata or MAE) <a name="ih_matrix_mae"></a>

```r
AnnotatedHeatmap(
  data,                    # Either:
                           #  1) numeric matrix (features × samples) + meta (if !MAE), OR
                           #  2) a MultiAssayExperiment (MAE) object
  meta                 = NULL,      # data.frame if data is matrix; rownames(meta)=colnames(data)
  assayName            = NULL,      # character, only if data is an MAE
  groupColumn          = NULL,      # grouping column for DE
  sampleType           = NULL,      # X-facet column for DE
  timepoint            = NULL,      # Y-facet column for DE (optional)
  dataType             = c("auto","continuous","count"), # counts vs continuous
  var_quantile         = 0.75,      # drop low-variance features
  pvalue_cutoff        = 0.05,      # “significant” threshold
  trending_cutoff      = 0.1,       # “trending” threshold
  fc_cutoff            = 0.585,     # |log2FC| threshold for including in matrix
  min_samples          = 3,         # min samples per group per facet
  max_features         = NULL,      # cap on features by variance or |logFC|
  runClustering        = FALSE,     # only used if data is MAE; run iClusterPlus?
  K                    = 3,         # # of clusters for iClusterPlus
  lambda               = 0.2,       # regularization for iClusterPlus
  coef                 = 2,         # contrast coefficient in limma
  pch_val              = 16,        # dot shape 
  unit_val             = 4,         # dot size in “mm”
  significant_color    = "black",   # color for p < pvalue_cutoff
  trending_color       = "yellow",  # color for p ∈ [pvalue_cutoff, trending_cutoff)
  heatmap_scale        = c("expression","logFC","foldChange"),
                             # if “expression”: tile fill = raw expression
                             # if “logFC” or “foldChange”: tile fill = logFC
  max_users            = NULL,      # future expansion
  …
)
```

* **Behavior**:

  1. **If** `inherits(data, "MultiAssayExperiment")`:

     * Extract assays → list of matrices (genes × samples).
     * **Optionally** run `iClusterPlus` on transposed assays to cluster samples.
     * For each assay, run `limma` to compare clusters → produce fold-change & p-value vectors.
     * **Combine** fold-change vectors across assays → `combined_logFC` (features × assays).
     * **Combine** p-value vectors → `combined_pvalues` (features × assays).
     * Decide `heatmap_data <- t(combined_logFC)` (samples/assays as rows/columns).
  2. **Else if** `is.matrix(data)` & `!is.null(meta)`:

     * **DE per facet** identical to `ThresholdedScatterplot_table`, but store final `logFC` & `pvalue` for each gene/sampleType/timePoint cell.
     * **Assemble** a matrix of fold-change: rows = features, columns = “cell” names (e.g., “StageI\_high”).
     * **Assemble** matching p-value matrix.
     * **Set** `heatmap_data <- logFC_matrix` (features × cells).
  3. **Else if** `is.matrix(data)` & `…` contains `pvalues = <matrix>`:

     * **Skip** DE; directly use provided `logFC_matrix = data` & `pvalue_matrix = pvalues`.
     * **Optionally** filter by `var_quantile` on `logFC_matrix` or by top `max_features` by absolute `logFC`.
  4. **Variance filtering** (if `var_quantile ∈ (0,1)`):

     ```r
     rv <- matrixStats::rowVars(heatmap_data, na.rm=TRUE)
     thr <- quantile(rv, var_quantile, na.rm=TRUE)
     keep <- rv >= thr
     heatmap_data <- heatmap_data[keep, , drop=FALSE]
     pvalue_matrix <- pvalue_matrix[keep, , drop=FALSE]
     ```
  5. **Max features** (if `max_features` numeric < `nrow(heatmap_data)`):

     ```r
     vv <- matrixStats::rowVars(heatmap_data, na.rm=TRUE)
     top_idx <- order(vv, decreasing=TRUE)[seq_len(max_features)]
     heatmap_data <- heatmap_data[top_idx, , drop=FALSE]
     pvalue_matrix <- pvalue_matrix[top_idx, , drop=FALSE]
     ```
  6. **Define `layer_fun`** (vectorized, no for-loops):

     ```r
     layer_fun <- function(j, i, x, y, w, h, fill) {
       ind_mat      <- ComplexHeatmap::restore_matrix(j, i, x, y)
       feature_idx  <- as.vector(row(ind_mat))
       cell_idx     <- as.vector(col(ind_mat))
       all_inds     <- as.vector(ind_mat)
       pvals_vec    <- pvalue_matrix[cbind(feature_idx, cell_idx)]
       color_vec    <- ifelse(
                         pvals_vec < pvalue_cutoff, significant_color,
                         ifelse(pvals_vec < trending_cutoff, trending_color, NA)
                       )
       keep_plots   <- !is.na(color_vec)
       grid::grid.points(
         x[all_inds][keep_plots],
         y[all_inds][keep_plots],
         pch  = pch_val,
         gp   = grid::gpar(col = color_vec[keep_plots]),
         size = grid::unit(unit_val, "mm")
       )
     }
     ```
  7. **Call** `ComplexHeatmap::Heatmap()`:

     ```r
     hm <- ComplexHeatmap::Heatmap(
             heatmap_data,
             name            = if (heatmap_scale[1] == "expression") "expr" else "logFC",
             layer_fun       = layer_fun,
             col             = circlize::colorRamp2(c(-2, 0, 2), c("blue","white","red")),
             cluster_rows    = TRUE,
             cluster_columns = TRUE,
             …
           )
     ```
  8. **Wrap**:

     ```r
     methods::new("AnnotatedHeatmap", heatmap = hm, params = as.list(environment()))
     ```

* **Exports**:

  ```r
  @exportClass AnnotatedHeatmap
  @exportMethod AnnotatedHeatmap
  ```

---

#### 5.2.2. AnnotatedHeatmapFromMAT (Fold-Change + P-Value Matrices) <a name="ih_matrices"></a>

```r
AnnotatedHeatmapFromMAT(
  logFC_matrix,      # numeric matrix: features × conditions
  pvalue_matrix,     # numeric matrix: same dims as logFC_matrix
  pvalue_cutoff     = 0.05,
  trending_cutoff   = 0.1,
  pch_val           = 16,
  unit_val          = 4,
  significant_color = "black",
  trending_color    = "yellow",
  col               = circlize::colorRamp2(c(-2, 0, 2), c("blue","white","red")),
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  show_row_names    = FALSE,
  show_column_names = FALSE,
  …
)
```

* **Checks & Sanity**:

  1. `stopifnot(is.matrix(logFC_matrix), is.matrix(pvalue_matrix))`
  2. `if (!all(dim(logFC_matrix) == dim(pvalue_matrix))) stop("logFC and pvalue matrices must have identical dimensions.")`
  3. **(Optional)** filter by `max_features` or user can subset before calling.
  4. **Define vectorized `layer_fun`** (as in previous method).
  5. **Call** `Heatmap(logFC_matrix, …, layer_fun = layer_fun, col = col, cluster_rows = cluster_rows, cluster_columns = cluster_columns, show_row_names=show_row_names, show_column_names=show_column_names, …)`.
  6. **Wrap**:

     ```r
     methods::new("AnnotatedHeatmap", heatmap = hm, params = as.list(environment()))
     ```

* **Exports**:

  ```r
  @export AnnotatedHeatmapFromMAT
  ```

---

#### 5.2.3. AnnotatedHeatmapFromMAE (MAE + iClusterPlus) <a name="ih_mae"></a>

```r
AnnotatedHeatmapFromMAE(
  mae,                    # MultiAssayExperiment object
  significant_pvalue   = 0.05,
  trending_pvalue      = 0.1,
  pch_val              = 16,
  unit_val             = 4,
  significant_color    = "black",
  trending_color       = "yellow",
  K                    = 3,      # # clusters for iClusterPlus
  lambda               = 0.2,    # regularization
  coef                 = 2,      # contrast index for limma
  heatmap_scale        = c("logFC","expression"),
  pvalue_cutoff        = 0.05,
  …
)
```

* **Workflow**:

  1. `assays(mae)` → list of assays (each a `SummarizedExperiment`).
  2. Convert each assay to a numeric matrix (genes × samples) and transpose → list of `sample × feature` matrices for iClusterPlus.
  3. Run `iClusterPlus` on that list:

     ```r
     fit <- iClusterPlus::iClusterPlus(
              data_list[["dt1"]], data_list[["dt2"]], …,
              type   = rep("gaussian", length(data_list)),
              K      = K,
              lambda = rep(lambda, length(data_list)),
              maxiter= 20
            )
     clusters <- factor(fit$clusters)  
     ```
  4. Build `design <- model.matrix(~ 0 + clusters)` and run `limma` on each assay:

     ```r
     limma_results <- lapply(names(assays(mae)), function(aname) {
       assay_data <- assay(mae, aname)  # genes × samples
       fit_limma  <- limma::lmFit(assay_data, design)
       fit_limma  <- limma::eBayes(fit_limma)
       topTab     <- limma::topTable(fit_limma, coef = coef, number = Inf)
       list(logFC   = topTab$logFC,
            p_values= topTab$P.Value)
     })
     names(limma_results) <- names(assays(mae))
     ```
  5. **Combine** logFCs & p-values across assays:

     ```r
     combined_logFC    <- do.call(rbind, lapply(limma_results, `[[`, "logFC"))
     combined_pvalues  <- do.call(rbind, lapply(limma_results, `[[`, "p_values"))
     if (nrow(combined_logFC) != nrow(combined_pvalues)) stop("Row mismatch between logFC & p-values.")
     heatmap_data     <- t(combined_logFC)  # samples/clusters × features
     pvalue_matrix    <- t(combined_pvalues)
     ```
  6. **Variance filter** & **max\_features** (optional).
  7. **Define vectorized `layer_fun`** as before.
  8. **Call** `Heatmap(heatmap_data, layer_fun=layer_fun, col=…, cluster_rows=TRUE, cluster_columns=TRUE, …)`.
  9. **Wrap**:

     ```r
     methods::new("AnnotatedHeatmap", heatmap = hm, params = as.list(environment()))
     ```

* **Dependencies**:

  * `iClusterPlus` (Bioconductor) for clustering.
  * `limma` for DE.
  * `matrixStats` for variance.

* **Exports**:

  ```r
  @export AnnotatedHeatmapFromMAE
  ```

---

#### 5.2.4. updateLayerFun (Replace Overlay Function) <a name="ih_update"></a>

```r
setMethod(
  "updateLayerFun",
  signature(x = "AnnotatedHeatmap"),
  function(x, layer_fun) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
      stop("ComplexHeatmap required to update layer_fun. Install via BiocManager::install('ComplexHeatmap').")
    params       <- x@params
    params$layer_fun <- layer_fun
    mat_data     <- x@heatmap@matrix  # original data used for Heatmap
    args         <- c(list(mat_data), params)
    new_ht       <- do.call(ComplexHeatmap::Heatmap, args)
    x@heatmap    <- new_ht
    x@params     <- params
    return(x)
  }
)
```

* **Purpose**: Take an existing `AnnotatedHeatmap` object, supply a brand-new `layer_fun` (e.g., plot different symbols or threshold criteria), and re-draw the heatmap.
* **Exports**:

  ```r
  @export updateLayerFun
  ```

---

#### 5.2.5. getHeatmapObject (Retrieve Underlying `Heatmap`) <a name="ih_get"></a>

```r
setMethod(
  "getHeatmapObject",
  signature(x = "AnnotatedHeatmap"),
  function(x) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
      stop("ComplexHeatmap required to retrieve Heatmap object.")
    return(x@heatmap)
  }
)
```

* **Purpose**: Extract the raw `Heatmap` object so you can:

  1. Add further annotations (e.g., row/column annotations),
  2. Use `ComplexHeatmap::draw()` or `ComplexHeatmap::ht_list` operations,
  3. Inspect heatmap internals.
* **Exports**:

  ```r
  @export getHeatmapObject
  ```

---

<a name="ih_examples"></a>

### 5.3. Practical Examples

Below are three lightweight examples showcasing use of **AnnotatedHeatmap** in different modes. Each example can be run in a fresh R session. We deliberately work with small subsets to minimize computation time.

---

<a name="ih_example2a"></a>

#### Example 2A: Subset of `miniACC` for Heatmap (Matrix + Metadata)

**Goal:** Build a heatmap of top DE genes (top 100 by variance) from a subset of the `miniACC` “RNASeq2GeneNorm” assay, faceted by “pathologic\_stage”.

> **Concept:** We will:
>
> 1. Extract a small subset of samples and features from `miniACC`.
> 2. Run `AnnotatedHeatmap` in “matrix + metadata” mode to perform DE within each `pathologic_stage` (two groups: “C1A” vs “C1B”).
> 3. Visualize a heatmap of logFC per gene per stage, overlaying points where p < 0.05.

1. **Load Packages & Data**

   ```r
   library(MultiModalGraphics)
   library(MultiAssayExperiment)
   library(SummarizedExperiment)
   library(matrixStats)
   library(BiocParallel)
   library(limma)
   library(ComplexHeatmap)
   library(dplyr)
   ```

2. **Extract `miniACC` & Subset Samples/Genes**

   ```r
   data("miniACC", package = "MultiAssayExperiment")
   se_rna <- experiments(miniACC)[["RNASeq2GeneNorm"]]
   expr_full <- assay(se_rna)             # ~21k × 79
   meta_full <- as.data.frame(colData(se_rna), stringsAsFactors = FALSE)

   # For speed, keep just 30 samples randomly
   set.seed(42)
   samp_keep <- sample(colnames(expr_full), 30)
   expr_small <- expr_full[, samp_keep]
   meta_small <- meta_full[samp_keep, , drop=FALSE]

   # Keep top 200 genes by variance across these 30 samples
   rv <- matrixStats::rowVars(expr_small, na.rm=TRUE)
   top_genes <- order(rv, decreasing = TRUE)[1:200]
   expr_small <- expr_small[top_genes, ]
   ```

3. **Run `AnnotatedHeatmap` in Matrix+Metadata Mode**

   ```r
   # Define parallel backend (2 cores)
   BPPARAM <- MulticoreParam(workers = 2)

   ih_matrix <- AnnotatedHeatmap(
     data               = expr_small,        # 200 × 30 matrix
     meta               = meta_small,       
     groupColumn        = "C1A.C1B",         # “C1A” vs “C1B”
     sampleType         = "pathologic_stage",# X-facet 
     timepoint          = NULL,              # single row 
     dataType           = "auto",            # assembly sees counts? Likely “continuous” 
     var_quantile       = 0.5,               # keep top 50% by variance (100 → 50 genes)
     pvalue_cutoff      = 0.05,
     trending_cutoff    = 0.1,
     fc_cutoff          = 0.585,
     min_samples        = 3,
     max_features       = 50,                # cap to top 50 genes 
     runClustering      = FALSE,             # no iClusterPlus here
     BPPARAM            = BPPARAM,
     pch_val            = 16,
     unit_val           = 3,
     significant_color  = "red",
     trending_color     = "orange",
     heatmap_scale      = "logFC",
     col                = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
     cluster_rows       = TRUE,
     cluster_columns    = FALSE,
     show_row_names     = TRUE,
     show_column_names  = TRUE,
     name               = "logFC"
   )
   ```

4. **Extract & Draw Heatmap**

   ```r
   ht_matrix <- getHeatmapObject(ih_matrix)
   ComplexHeatmap::draw(ht_matrix, heatmap_legend_side = "right")
   ```

   **What you see:**

   * A heatmap (50 genes × “pathologic\_stage” facets).
   * Tile fill = logFC for C1A vs C1B for each pathologic\_stage (one column per stage).
   * Overlay points (red/orange) on cells whose p‐value is < 0.05 or < 0.1.
   * Row dendrogram clustering genes; no column clustering.

---

<a name="ih_example2b"></a>

#### Example 2B: Small Simulated Fold-Change + P-Value Matrix

**Goal:** Create a heatmap from a precomputed 15×4 fold-change matrix and matching p-value matrix.

1. **Simulate Data**

   ```r
   set.seed(123)
   # 15 features (rows) × 4 conditions (columns)
   logFC_mat    <- matrix(rnorm(60, sd = 1), nrow = 15, ncol = 4)
   pval_mat     <- matrix(runif(60, 0, 0.2), nrow = 15, ncol = 4)
   rownames(logFC_mat) <- paste0("Gene", 1:15)
   colnames(logFC_mat) <- paste0("Cond", 1:4)
   rownames(pval_mat) <- rownames(logFC_mat)
   colnames(pval_mat) <- colnames(logFC_mat)
   ```

2. **Construct `AnnotatedHeatmapFromMAT` Object**

   ```r
   ih_mat <- AnnotatedHeatmapFromMAT(
     logFC_matrix     = logFC_mat,
     pvalue_matrix    = pval_mat,
     pvalue_cutoff    = 0.05,
     trending_cutoff  = 0.1,
     pch_val          = 17,
     unit_val         = 3,
     significant_color= "darkred",
     trending_color   = "darkorange",
     col              = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
     cluster_rows     = TRUE,
     cluster_columns  = TRUE,
     show_row_names   = TRUE,
     show_column_names= TRUE
   )
   ```

3. **Draw Heatmap**

   ```r
   ht_mat <- getHeatmapObject(ih_mat)
   ComplexHeatmap::draw(ht_mat, heatmap_legend_side = "right")
   ```

   **What you see:**

   * A heatmap (15 features × 4 conditions).
   * Tile fill = `logFC`.
   * Red dots for p < 0.05; orange dots for 0.05 ≤ p < 0.1 overlaid on corresponding cells.

---

<a name="ih_example2c"></a>

#### Example 2C: Custom MAE with `AnnotatedHeatmapFromMAE`

**Goal:** Simulate a tiny MAE with two assays (gene expression & methylation), run iClusterPlus to identify 2 clusters, run DE on each assay comparing clusters, and build a combined heatmap.

> **Note:** iClusterPlus can be slow. We simulate a very small dataset:  50 features × 8 samples for each assay and 2 clusters.

1. **Simulate Two Assays & Build Tiny MAE**

   ```r
   library(MultiAssayExperiment)
   library(SummarizedExperiment)

   set.seed(456)
   # Simulate 50 genes × 8 samples expression
   expr1 <- matrix(rnorm(400, mean=5, sd=2), nrow=50, ncol=8)
   rownames(expr1) <- paste0("Gene", 1:50)
   colnames(expr1) <- paste0("S", 1:8)
   se_expr <- SummarizedExperiment(assays = list(expr = expr1),
                                   colData = DataFrame(row.names = paste0("S", 1:8)))

   # Simulate 50 “methylation” features × same 8 samples
   meth1 <- matrix(rnorm(400, mean=0.5, sd=0.1), nrow=50, ncol=8)
   rownames(meth1) <- paste0("Meth", 1:50)
   colnames(meth1) <- paste0("S", 1:8)
   se_meth <- SummarizedExperiment(assays = list(beta = meth1),
                                   colData = DataFrame(row.names = paste0("S", 1:8)))

   # Build sampleMap and MAE
   sample_map <- DataFrame(
     assay   = rep(c("Expr","Meth"), each=8),
     primary = rep(paste0("S", 1:8), 2),
     colname = rep(paste0("S", 1:8), 2)
   )
   tiny_mae <- MultiAssayExperiment(
     experiments = SimpleList(Expr = se_expr, Meth = se_meth),
     colData     = DataFrame(row.names = paste0("S", 1:8)),
     sampleMap   = sample_map
   )
   ```

2. **Run `AnnotatedHeatmapFromMAE`**

   ```r
   # Use 2 clusters, lambda=0.1
   ih_mae <- AnnotatedHeatmapFromMAE(
     mae                 = tiny_mae,
     significant_pvalue  = 0.05,
     trending_pvalue     = 0.1,
     pch_val             = 16,
     unit_val            = 2,
     significant_color   = "darkred",
     trending_color      = "darkorange",
     K                   = 2,
     lambda              = 0.1,
     coef                = 2,
     heatmap_scale       = "logFC",
     col                 = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
     cluster_rows        = TRUE,
     cluster_columns     = TRUE,
     show_row_names      = TRUE,
     show_column_names   = TRUE,
     name                = "logFC"
   )
   ```

3. **Draw Heatmap**

   ```r
   ht_mae <- getHeatmapObject(ih_mae)
   ComplexHeatmap::draw(ht_mae, heatmap_legend_side = "right")
   ```

   **What you see:**

   * Heatmap of combined logFC (features × 2 clusters) for both assays stacked (iClusterPlus).
   * Dots where p < 0.05 (dark red) or 0.05 ≤ p < 0.1 (dark orange).

---

<a name="ih_advanced"></a>

### 5.4. Advanced Customization

1. **Custom Color Ramps**
   Provide any `circlize::colorRamp2()` call:

   ```r
   ih_matrix <- AnnotatedHeatmap(
     …,
     col = circlize::colorRamp2(c(-3,0,3), c("#1B9E77","white","#D95F02")),
     …
   )
   ```

2. **Turn Off Clustering**

   ```r
   ih_mat <- AnnotatedHeatmapFromMAT(
     …,
     cluster_rows       = FALSE,
     cluster_columns   = FALSE,
     …
   )
   ```

3. **Adjust Dot Shapes & Sizes**

   ```r
   ih_mat <- AnnotatedHeatmapFromMAT(
     …,
     pch_val = 15,  # square dots
     unit_val= 5,   # bigger
     …
   )
   ```

4. **Change “Significant” vs “Trending” Thresholds**

   ```r
   ih_matrix <- AnnotatedHeatmap(
     …,
     pvalue_cutoff   = 0.01,
     trending_cutoff = 0.05,
     …
   )
   ```

5. **Update Layer Function** (e.g., highlight only p < 0.001)

   ```r
   new_layer_fun <- function(j, i, x, y, w, h, fill) {
     ind_mat      <- ComplexHeatmap::restore_matrix(j, i, x, y)
     feature_idx  <- as.vector(row(ind_mat))
     cell_idx     <- as.vector(col(ind_mat))
     all_inds     <- as.vector(ind_mat)
     # pvalue_matrix is stored in ih_mat@params$pvalue_matrix
     pvals_vec    <- ih_mat@params$pvalue_matrix[cbind(feature_idx, cell_idx)]
     color_vec    <- ifelse(pvals_vec < 0.001, "magenta", NA)
     keep_plots   <- !is.na(color_vec)
     grid::grid.points(
       x[all_inds][keep_plots],
       y[all_inds][keep_plots],
       pch = 16,
       gp  = grid::gpar(col = color_vec[keep_plots]),
       size= grid::unit(4, "mm")
     )
   }
   ih_matrix <- updateLayerFun(ih_matrix, new_layer_fun)
   ComplexHeatmap::draw(getHeatmapObject(ih_matrix))
   ```

---

<a name="ih_troubleshoot"></a>

### 5.5. Troubleshooting & Tips

1. **“Data must be a matrix”**

   * *Symptom*:

     ```
     Error in AnnotatedHeatmap(data = df, …): Data must be a matrix.
     ```
   * *Fix*: If you intended to run the “table-only” or “precomputed” workflow, use `AnnotatedHeatmapFromMAT`. To use `AnnotatedHeatmap` with a data.frame, first convert to matrix.

2. **“Missing pvalues”**

   * *Symptom*:

     ```
     Error: object 'pvalue_matrix' not found
     ```
   * *Fix*: For `AnnotatedHeatmapFromMAT`, ensure you passed `pvalue_matrix=…` exactly matching dims of `logFC_matrix`.

3. **“iClusterPlus took too long”**

   * *Symptom*:
     Long wait or “NA/NaN/Inf in foreign function call” if parameters are not sensible.
   * *Fix*:

     * Use fewer features or fewer dimensions (e.g., set `var_quantile=0.9`, `max_features=1000`).
     * Check `K` and `lambda` values are appropriate.
     * For small pilot runs, use a tiny MAE (50–100 features, 8–12 samples).

4. **“Heatmap intermittently blank”**

   * *Symptom*:
     No colors appear or no dots appear.
   * *Fix*:

     * Verify `heatmap_data` and `pvalue_matrix` are non-empty.
     * Ensure thresholds are not too strict (e.g., `pvalue_cutoff=0.001` might yield zero points).
     * Check `cluster_rows`/`cluster_columns` arguments.

5. **Missing `ComplexHeatmap` or `circlize` functions**

   * *Symptom*:

     ```
     Error: could not find function "Heatmap"
     ```
   * *Fix*:

     ```r
     library(ComplexHeatmap)
     library(circlize)
     ```

---

<a name="multifeature_grid"></a>

## 6. Section 3: CompositeFeatureHeatmap (2D Tile-and-Point Heatmaps)

<a name="mg_background"></a>

### 6.1. Background & Purpose

A **CompositeFeatureHeatmap** plot displays a grid where:

* **X-axis**: one categorical variable (e.g., tissue type).
* **Y-axis**: another categorical variable (e.g., signaling pathway).
* **Tile fill**: an activation z-score (or other continuous metric).
* **Point overlay**: colored by –log<sub>10</sub>(p-value) and sized by number of items (e.g., genes) associated with that tile.
* **Faceting**: optionally by a third categorical variable (e.g., timePoint).

This plot is especially useful for summarizing pathway-level activation (z-scores), with p-value significance and size indicating gene counts, across tissues or conditions.

Core slots:

* `@data`: A `data.frame` with at least:

  * `x` (categorical, e.g. “tissue”)
  * `y` (categorical, e.g. “signaling”)
  * `Activation_z_score` (numeric)
  * `p` (numeric p-value)
  * `number_of_genes` (numeric)
  * `timePoint` (optional, for facet)

* `@title`, `@x_label`, `@y_label`, `@logpval_label`, `@zscore_label`, `@numitems_label`, `@color_palette`, `@breaks`: plot configuration.

Two key methods:

* Constructor: `CompositeFeatureHeatmap(data_frame, …)`.
* Plot: `plot_heatmap(CompositeFeatureHeatmap_object, pValueColumn, lowColor, highColor, borderColor, columnForNumber, independantVariable)`.

---

<a name="mg_api"></a>

### 6.2. API Reference

#### 6.2.1. CompositeFeatureHeatmap (Constructor) <a name="mg_constructor"></a>

```r
CompositeFeatureHeatmap(
  data,                # data.frame with required columns
  title            = "Heatmap",
  x_label          = "X Label",
  y_label          = "Y Label",
  logpval_label    = "-log10(p-value)",
  zscore_label     = "Activation z-score",
  numitems_label   = "Number of Genes",
  color_palette    = "RdYlBu",   # valid RColorBrewer palette
  breaks           = seq(-1, 1, 0.5)
)
```

* **Checks & Sanity**:

  1. `stopifnot(is.data.frame(data))`.
  2. Ensure `c("tissue","signaling","Activation_z_score","p","number_of_genes") %in% colnames(data)`.
  3. If `timePoint` (or user-specified faceting column) is used, ensure `timePoint %in% colnames(data)`.
  4. Confirm `color_palette` is valid:

     ```r
     if (!color_palette %in% rownames(RColorBrewer::brewer.pal.info)) {
       stop("Invalid RColorBrewer palette: ", color_palette)
     }
     ```
  5. Confirm `breaks` is a numeric vector of length ≥ 2.
  6. Return:

     ```r
     methods::new(
       "CompositeFeatureHeatmap",
       data            = data,
       title           = title,
       x_label         = x_label,
       y_label         = y_label,
       logpval_label   = logpval_label,
       zscore_label    = zscore_label,
       numitems_label  = numitems_label,
       color_palette   = color_palette,
       breaks          = breaks
     )
     ```

* **Exports**:

  ```r
  @exportClass CompositeFeatureHeatmap
  @exportMethod CompositeFeatureHeatmap
  ```

---

#### 6.2.2. plot\_heatmap (Tile-and-Point Plot) <a name="mg_plot"></a>

```r
setMethod(
  "plot_heatmap",
  signature(object = "CompositeFeatureHeatmap"),
  function(
    object,
    pValueColumn      = "p",                   # column for p-values
    lowColor          = "yellow",              # color for low –log10(p)
    highColor         = "red",                 # color for high –log10(p)
    borderColor       = "grey60",              # tile border color
    columnForNumber   = "number_of_genes",     # column for point size 
    independantVariable= "timePoint"           # column for facet rows
  ) {
    data <- object@data
    title <- object@title
    x_label <- object@x_label
    y_label <- object@y_label
    logpval_label <- object@logpval_label
    zscore_label <- object@zscore_label
    numitems_label <- object@numitems_label
    color_palette <- object@color_palette
    breaks <- object@breaks

    # Check columns
    required_cols <- c("tissue","signaling","Activation_z_score", pValueColumn, columnForNumber)
    missing_cols  <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns in data: ", paste(missing_cols, collapse=", "))
    }
    if (!is.factor(data$tissue)) {
      data$tissue <- factor(data$tissue)
    }
    if (!is.factor(data$signaling)) {
      data$signaling <- factor(data$signaling)
    }
    if (!is.numeric(data$Activation_z_score)) {
      stop("`Activation_z_score` must be numeric.")
    }
    if (!is.numeric(data[[pValueColumn]])) {
      stop("`", pValueColumn, "` must be numeric.")
    }
    if (!is.numeric(data[[columnForNumber]])) {
      stop("`", columnForNumber, "` must be numeric.")
    }
    if (!independantVariable %in% colnames(data)) {
      stop("`independantVariable` not found in data: ", independantVariable)
    }
    if (!is.factor(data[[independantVariable]])) {
      data[[independantVariable]] <- factor(data[[independantVariable]])
    }

    # Calculate –log10(p)
    data$neglog10p <- -log10(data[[pValueColumn]])

    # Define color ramp for tile fill
    pal_colors <- rev(RColorBrewer::brewer.pal(n = max(3, length(breaks)-1), name = color_palette))
    fill_colors <- grDevices::colorRampPalette(pal_colors)(100)

    # Build ggplot
    p <- ggplot2::ggplot(data, ggplot2::aes(x = tissue, y = signaling)) +
         ggplot2::geom_tile(ggplot2::aes(fill = Activation_z_score), colour = borderColor) +
         ggplot2::scale_fill_gradientn(
           colours = fill_colors,
           breaks  = breaks,
           labels  = scales::comma_format()
         ) +
         ggplot2::geom_point(
           ggplot2::aes(
             colour = neglog10p,
             size   = .data[[columnForNumber]]
           )
         ) +
         ggplot2::scale_color_gradient(low = lowColor, high = highColor) +
         ggplot2::scale_size(range = c(1, 10)) +
         ggplot2::labs(
           x      = x_label,
           y      = y_label,
           title  = title,
           fill   = zscore_label,
           colour = logpval_label,
           size   = numitems_label
         ) +
         ggplot2::facet_grid(
           cols = vars(.data[[independantVariable]]),
           scales = "free_x",
           space  = "free"
         ) +
         ggplot2::theme_bw() +
         ggplot2::theme(
           axis.text.x       = ggplot2::element_text(angle = 90, hjust = 1),
           panel.border      = ggplot2::element_rect(fill = NA, colour = "grey80", linewidth = 0.6),
           axis.text         = ggplot2::element_text(size = 14, face = "bold"),
           axis.title        = ggplot2::element_text(size = 18, face = "bold"),
           title             = ggplot2::element_text(size = 18),
           strip.text.x      = ggplot2::element_text(size = 14, face = "bold", colour = "black")
         )

    print(p)
  }
)
```

* **Parameters**:

  * `pValueColumn`: Column name in `data` for p-values (numeric).
  * `lowColor`, `highColor`: Colors for low/high –log<sub>10</sub>(p) gradient.
  * `borderColor`: Color for tile borders.
  * `columnForNumber`: Column name for “number\_of\_genes” (size mapping).
  * `independantVariable`: Column name for facet (rows).

* **Exports**:

  ```r
  @export plot_heatmap
  ```

---

<a name="mg_examples"></a>

### 6.3. Practical Examples

Below are two lightweight, **real-world-style** examples demonstrating `CompositeFeatureHeatmap`. Both use a small “demo” dataset to minimize memory/time.

---

<a name="mg_example3a"></a>

#### Example 3A: Built-in Demo Data via `get_multifeature_grid_df()`

**Goal:** Show a tile-and-point heatmap of pathway activation data (z-scores, p-values, gene counts) faceted by time point.

> **MultiModalGraphics** includes a helper `get_multifeature_grid_df()` that returns a small example data.frame (≈ 16 rows). We use that.

1. **Load Packages & Demo Data**

   ```r
   library(MultiModalGraphics)

   # Assume get_multifeature_grid_df() is exported
   data_mg <- get_multifeature_grid_df()  
   str(data_mg)
   # 'data.frame': e.g. 16 observations × 6 columns:
   #  $ tissue            : Factor w/ 4 levels “Tissue1”, “Tissue2”, …
   #  $ signaling         : Factor w/ 4 levels “Pathway1”, …
   #  $ Activation_z_score: num
   #  $ p                 : num
   #  $ number_of_genes   : int
   #  $ timePoint         : Factor w/ 2 levels “Time1”, “Time2”
   ```

2. **Construct `CompositeFeatureHeatmap` Object**

   ```r
   mg <- CompositeFeatureHeatmap(
     data           = data_mg,
     title          = "Demo Multifeature Grid",
     x_label        = "Tissue Type",
     y_label        = "Signaling Pathway",
     logpval_label  = "-Log10 P-Value",
     zscore_label   = "Activation Z-Score",
     numitems_label = "Gene Count",
     color_palette  = "RdYlBu",
     breaks         = seq(-2, 2, 0.5)
   )
   ```

3. **Plot the Heatmap**

   ```r
   plot_heatmap(
     mg,
     pValueColumn      = "p",
     lowColor          = "lightblue",
     highColor         = "darkred",
     borderColor       = "black",
     columnForNumber   = "number_of_genes",
     independantVariable= "timePoint"
   )
   ```

   **What you see:**

   * Tiles colored by z-score (blue to red via “RdYlBu”).
   * Points sized by gene count (1–10 range) and colored by –log10(p) (lightblue→darkred).
   * Facets (columns) for `Time1` and `Time2`.
   * X-axis = tissue (“Tissue1”, “Tissue2”, …); Y-axis = signaling pathways.

---

<a name="mg_example3b"></a>

#### Example 3B: Custom Table from Pathway Enrichment (Real-World)

**Goal:** Suppose you ran pathway enrichment on Differential Expression results from the “Taylor” dataset and obtained a small table of 12 pathways. We demonstrate building a `CompositeFeatureHeatmap` from that real-world output.

1. **(Simulate) Pathway Enrichment**

   > In practice, you would run e.g. `clusterProfiler::enrichKEGG()` or `fgsea::fgsea()`. Here, we simulate a simplified output:

   ```r
   set.seed(789)
   pathways <- paste0("Pathway", 1:12)
   tissues  <- rep(c("Normal","Tumor"), each = 12)
   # Each pathway appears in both tissues
   df_pe <- data.frame(
     tissue            = factor(rep(c("Normal","Tumor"), each = 12)),
     signaling         = factor(rep(pathways, 2)),
     Activation_z_score= runif(24, -3, 3),
     p                 = runif(24, 1e-5, 0.05),
     number_of_genes   = sample(10:200, 24, TRUE),
     timePoint         = factor(sample(c("Baseline","FollowUp"), 24, replace = TRUE)),
     stringsAsFactors  = FALSE
   )
   df_pe$tissue    <- factor(df_pe$tissue, levels = c("Normal","Tumor"))
   df_pe$signaling <- factor(df_pe$signaling, levels = pathways)
   ```

2. **Construct & Plot**

   ```r
   mg2 <- CompositeFeatureHeatmap(
     data           = df_pe,
     title          = "Pathway Enrichment Grid: Taylor PCa",
     x_label        = "Tissue",
     y_label        = "Pathway",
     logpval_label  = "-Log10 P-Value",
     zscore_label   = "Z-Score",
     numitems_label = "Gene Count",
     color_palette  = "PuOr",         # diverging
     breaks         = seq(-3, 3, 1)
   )

   plot_heatmap(
     mg2,
     pValueColumn      = "p",
     lowColor          = "lightgreen",
     highColor         = "darkgreen",
     borderColor       = "grey50",
     columnForNumber   = "number_of_genes",
     independantVariable= "timePoint"
   )
   ```

   **What you see:**

   * Tiles colored by `Activation_z_score` (green→white→purple).
   * Points sized by gene count and colored by –log10(p).
   * Facets for “Baseline” vs “FollowUp.”

---

<a name="mg_advanced"></a>

### 6.4. Advanced Customization

1. **Custom Breakpoints & Palettes**

   ```r
   mg <- CompositeFeatureHeatmap(
     data           = data_mg,
     color_palette  = "Spectral",
     breaks         = seq(-4,4,0.5)
   )
   ```

2. **Removing Facets**
   If your data has no `timePoint`, set `independantVariable = NULL` (or a column that is constant).

   ```r
   plot_heatmap(mg, independantVariable = NULL)
   ```

3. **Adjust Tile Borders**

   ```r
   plot_heatmap(mg, borderColor = "white", lowColor = "cyan", highColor = "magenta")
   ```

4. **Change Point Aesthetics**
   Control point range:

   ```r
   p <- plot_heatmap(
          mg,
          lowColor        = "blue",
          highColor       = "red",
          columnForNumber = "number_of_genes"
        )
   # To manually adjust point range, modify scale_size:
   p + scale_size(range = c(2, 15))
   ```

5. **Axis & Font Styling**
   After `plot_heatmap()`, you can add further `ggplot2` modifications:

   ```r
   p <- plot_heatmap(mg)
   p + theme(
     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
     axis.text.y = element_text(size = 12),
     plot.title  = element_text(size = 20, face = "bold")
   )
   ```

---

<a name="mg_troubleshoot"></a>

### 6.5. Troubleshooting & Tips

1. **“Missing Required Columns”**

   * *Symptom*:

     ```
     Error in plot_heatmap(mg, …): Missing required columns in data: …
     ```
   * *Fix*: Rename your data columns to include at least `tissue`, `signaling`, `Activation_z_score`, your `pValueColumn` (default `"p"`), and `columnForNumber` (default `"number_of_genes"`).

2. **Invalid Palette Name**

   * *Symptom*:

     ```
     Error: Invalid RColorBrewer palette: "XYZ"
     ```
   * *Fix*: Choose a palette from `rownames(RColorBrewer::brewer.pal.info)` (e.g., “RdYlBu”, “PuOr”, “Spectral”, etc.)

3. **Axis Text Overlapping**

   * If your `tissue` or `signaling` has long labels, add:

     ```r
     + theme(axis.text.x = element_text(angle = 45, hjust = 1))
     ```

4. **Large Dataframes**

   * If you have many pathways & tissues (e.g., 100 × 20), faceting can become cluttered. Consider:

     * Reducing levels of `timePoint`.
     * Splitting into multiple panels.
     * Decreasing tile size or removing tile borders.

---

<a name="multimodal_plot"></a>

## 7. Section 4: Unified Multimodal API (`MultiModalPlot`)

<a name="mm_background"></a>

### 7.1. Background & Purpose

Often, you want to display **multiple visualizations** (e.g., several volcano plots or several heatmaps) side-by-side for different modalities (RNA, methylation, proteomics). Instead of manually invoking each constructor and combining with `patchwork` or `ComplexHeatmap::+`, **MultiModalPlot** automatically:

1. Detects each input’s type (MAE, matrix+meta, DE table).
2. Dispatches to the appropriate builder (`ThresholdedScatterplot_MAE`, `ThresholdedScatterplot_table`, `ThresholdedScatterplot`, or `AnnotatedHeatmap`).
3. Returns a **combined plot**:

   * For `panel_type = "volcano"`: a `ggplot` (if one modality) or a `patchwork` object (if >1).
   * For `panel_type = "heatmap"`: a `Heatmap` or `HeatmapList` (if >1), ready for `draw()`.

Thus, `MultiModalPlot` is a **single entry point** for complex, multimodal displays.

---

<a name="mm_api"></a>

### 7.2. API Reference

```r
MultiModalPlot(
  inputs,                # named list of inputs: each element is either:
                         #   • a MultiAssayExperiment (MAE), or
                         #   • a list(expr = <matrix>, meta = <data.frame>), or
                         #   • a precomputed DE data.frame
  assayNames     = NULL, # named character vector: for each MAE input, which assay to use
  groupColumns   = NULL, # named character vector: grouping column per modality
  sampleTypes    = NULL, # named character vector: X-facet column per modality
  timepoints     = NULL, # named character vector: Y-facet column per modality (optional)
  dataType       = "auto",    # “auto”/“continuous”/“count” per modality
  var_quantile   = 0.75,      # variance filter threshold per modality
  pvalue_cutoff  = 0.05,      # –log10(p) threshold per modality
  trending_cutoff= 0.1,       # trending threshold for heatmaps
  fc_cutoff      = 0.585,     # fold-change threshold per modality
  max_features   = NULL,      # cap on features per modality
  parallel       = TRUE,      # whether to parallelize per modality
  BPPARAM        = BiocParallel::bpparam(),
  panel_type     = c("volcano","heatmap"),
  …
)
```

* **Inputs**: Must be a **named** list. Each element’s name is used as the panel label.

  * **Case A**: `inherits(x, "MultiAssayExperiment")` → calls `ThresholdedScatterplot_MAE()` or `AnnotatedHeatmap()` depending on `panel_type`.
  * **Case B**: `is.list(x) && is.matrix(x$expr) && is.data.frame(x$meta)` → calls `ThresholdedScatterplot_table()` or `AnnotatedHeatmap()`.
  * **Case C**: `is.data.frame(x)` with columns `c("log2fc","negLog10p","regulation","SampleType")` → calls `ThresholdedScatterplot()`. Only valid if `panel_type = "volcano"`.
  * Else: `stop("Input ‘…’ not recognized. Must be MAE, (expr, meta), or DE table.")`.

* **Per-modality arguments**: If you pass vectors (e.g., `dataType = c(RNA="auto", Meth="count")`), `MultiModalPlot` picks the element matching the modality name. Otherwise, it recycles or defaults.

* **panel\_type**:

  * `"volcano"` → returns a `ggplot` if only one modality, or a `patchwork`‐combined plot if >1.
  * `"heatmap"` → returns a single `Heatmap` (if one modality) or a `HeatmapList` combined with `+` (side-by-side) if >1.

* **Additional arguments** (`…`):

  * For volcano mode: passed to `createPlot()` (e.g., `color_up`, `color_down`, `title`, etc.).
  * For heatmap mode: passed to the respective `AnnotatedHeatmap*` constructor (e.g., `col`, `cluster_rows`, `pch_val`, etc.).

* **Example skeleton**:

  ```r
  if (!is.list(inputs) || is.null(names(inputs))) {
    stop("`inputs` must be a named list.")
  }
  for (mod in names(inputs)) {
    x <- inputs[[mod]]
    if (inherits(x, "MultiAssayExperiment")) {
      # check assayNames[[mod]], groupColumns[[mod]], sampleTypes[[mod]]
      if (panel_type=="volcano") {
        panel_list[[mod]] <- ThresholdedScatterplot_MAE(
          mae          = x,
          assayName    = assayNames[[mod]],
          groupColumn  = groupColumns[[mod]],
          sampleType   = sampleTypes[[mod]],
          timepoint    = timepoints[[mod]],
          dataType     = dataType[[mod]],
          var_quantile = var_quantile,
          pvalue_cutoff= pvalue_cutoff,
          fc_cutoff    = fc_cutoff,
          max_features = max_features,
          parallel     = parallel,
          BPPARAM      = BPPARAM,
          …
        )
      } else {
        panel_list[[mod]] <- AnnotatedHeatmap(
          data            = x,
          assayName       = assayNames[[mod]],
          groupColumn     = groupColumns[[mod]],
          sampleType      = sampleTypes[[mod]],
          timepoint       = timepoints[[mod]],
          dataType        = dataType[[mod]],
          var_quantile    = var_quantile,
          pvalue_cutoff   = pvalue_cutoff,
          trending_cutoff = trending_cutoff,
          fc_cutoff       = fc_cutoff,
          max_features    = max_features,
          runClustering   = TRUE,
          K               = 3,
          lambda          = 0.2,
          BPPARAM         = BPPARAM,
          …
        )
      }
    } else if (is.list(x) && is.matrix(x$expr) && is.data.frame(x$meta)) {
      # call ThresholdedScatterplot_table() or AnnotatedHeatmap()
      …
    } else if (is.data.frame(x) &&
               all(c("log2fc","negLog10p","regulation","SampleType") %in% colnames(x))) {
      if (panel_type != "volcano") {
        stop("DE tables only support `panel_type='volcano'`.")
      }
      panel_list[[mod]] <- ThresholdedScatterplot(
        data           = x,
        highLog2fc     = fc_cutoff,
        lowLog2fc      = -fc_cutoff,
        negLog10pValue = -log10(pvalue_cutoff),
        …
      )
    } else {
      stop("Input ‘", mod, "’ not recognized.")
    }
  }
  ```

* **Combining Panels**:

  * **Volcano**:

    ```r
    library(patchwork)
    gglist     <- lapply(panel_list, function(cs) {
                   cs2 <- createPlot(cs, …)
                   cs2@plot
                 })
    combined_gg <- Reduce(`+`, gglist)
    return(combined_gg)
    ```
  * **Heatmap**:

    ```r
    ht_list <- lapply(panel_list, getHeatmapObject)
    if (length(ht_list)==1) {
      return(ht_list[[1]])
    } else {
      combined_ht <- Reduce(`+`, ht_list)
      return(combined_ht)
    }
    ```

* **Exports**:

  ```r
  @export MultiModalPlot
  ```

---

<a name="mm_examples"></a>

### 7.3. Practical Examples

Below are three examples illustrating `MultiModalPlot` for multimodal integration.

---

<a name="mm_example4a"></a>

#### Example 4A: Single MAE Volcano (Volcano Mode)

**Goal:** One-line volcano plot from `miniACC`, same as Example 1A but via `MultiModalPlot`.

```r
library(MultiModalGraphics)
library(MultiAssayExperiment)
library(BiocParallel)

data("miniACC", package="MultiAssayExperiment")
BPPARAM <- MulticoreParam(workers = 2)

volcano_plot <- MultiModalPlot(
  inputs        = list(ACC = miniACC),
  assayNames    = c(ACC = "RNASeq2GeneNorm"),
  groupColumns  = c(ACC = "C1A.C1B"),
  sampleTypes   = c(ACC = "pathologic_stage"),
  timepoints    = c(ACC = "MethyLevel"),
  panel_type    = "volcano",
  dataType      = "auto",
  var_quantile  = 0.75,
  pvalue_cutoff = 0.05,
  fc_cutoff     = 0.585,
  max_features  = 5000,
  parallel      = TRUE,
  BPPARAM       = BPPARAM,
  # Volcano plot custom args:
  color_up        = "tomato",
  color_neutral   = "lightgrey",
  color_down      = "steelblue",
  xlab            = "Log2 FC",
  ylab            = "-Log10 P",
  text_family     = "Times",
  text_size       = 12,
  point_size      = 2,
  point_alpha     = 0.7,
  legend_position = "right",
  legend_title    = "Regulation",
  legend_labels   = c("Up","Neutral","Down"),
  custom_theme    = theme_minimal(),
  title           = "miniACC Volcano (MultiModalPlot)",
  subtitle        = "C1A vs C1B by stage & methylation",
  caption         = "Data: miniACC"
)

print(volcano_plot)
```

**What happens internally**:

* Detects `inputs$ACC` is an MAE → calls `ThresholdedScatterplot_MAE(...)` → builds `ThresholdedScatterplot` object → `createPlot(...)` → `ggplot`.
* Returns a `ggplot` object ready to print.

---

<a name="mm_example4b"></a>

#### Example 4B: Two Modalities, Combined Heatmap (Heatmap Mode)

**Goal:** Suppose you have RNA expression and DNA methylation MAEs for the same 40 samples. We build two heatmap panels side-by-side comparing Tumor vs Normal DE for each modality.

> **Note:** For demonstration, we subset `miniACC` to RNA (already present) and simulate a matching methylation MAE with the same samples.

1. **Load & Prepare Data**

   ```r
   library(MultiModalGraphics)
   library(MultiAssayExperiment)
   library(SummarizedExperiment)
   library(matrixStats)
   library(BiocParallel)
   library(limma)
   library(ComplexHeatmap)
   library(circlize)

   # 1) Use miniACC RNA assay (30 samples)
   data("miniACC", package="MultiAssayExperiment")
   se_rna_full <- experiments(miniACC)[["RNASeq2GeneNorm"]]
   samp30 <- sample(colnames(assay(se_rna_full)), 30)
   expr_rna30 <- assay(se_rna_full)[, samp30]
   meta_rna30 <- as.data.frame(colData(se_rna_full)[samp30, ], stringsAsFactors=FALSE)
   se_rna30   <- SummarizedExperiment(assays = list(expr = expr_rna30),
                                      colData = DataFrame(meta_rna30, row.names = samp30))

   # 2) Simulate a methylation assay for same 30 samples
   set.seed(101)
   expr_meth30 <- matrix(runif(50*30, 0, 1), nrow = 50, ncol = 30)
   rownames(expr_meth30) <- paste0("CpG", 1:50)
   colnames(expr_meth30) <- samp30
   se_meth30 <- SummarizedExperiment(assays = list(beta = expr_meth30),
                                     colData = DataFrame(meta_rna30, row.names = samp30))

   # 3) Build sampleMap & MAEs
   smap <- DataFrame(
     assay   = rep(c("RNA","Meth"), each = 30),
     primary = rep(samp30, 2),
     colname = rep(samp30, 2)
   )
   mae_rna_meth <- MultiAssayExperiment(
     experiments = SimpleList(RNA = se_rna30, Meth = se_meth30),
     colData     = DataFrame(meta_rna30, row.names = samp30),
     sampleMap   = smap
   )
   ```

2. **Run `MultiModalPlot` in Heatmap Mode**

   ```r
   BPPARAM <- MulticoreParam(workers = 2)

   ht_all <- MultiModalPlot(
     inputs        = list(RNA = mae_rna_meth, Meth = mae_rna_meth),
     assayNames    = c(RNA = "expr", Meth = "beta"),
     groupColumns  = c(RNA = "C1A.C1B", Meth = "C1A.C1B"),
     sampleTypes   = c(RNA = "pathologic_stage", Meth = "pathologic_stage"),
     panel_type    = "heatmap",
     dataType      = c(RNA="continuous", Meth="continuous"),
     var_quantile  = 0.75,
     pvalue_cutoff = 0.05,
     trending_cutoff= 0.1,
     fc_cutoff     = 0.585,
     max_features  = 30,
     runClustering = TRUE,
     K             = 2,
     lambda        = 0.1,
     BPPARAM       = BPPARAM,
     heatmap_scale = "logFC",
     # Custom color palettes per modality
     col_RNA       = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
     col_Meth      = circlize::colorRamp2(c(-2,0,2), c("purple","white","green")),
     cluster_rows  = TRUE,
     cluster_columns = TRUE,
     show_row_names  = TRUE,
     show_column_names = FALSE
   )

   # Draw combined heatmap list
   ComplexHeatmap::draw(ht_all, heatmap_legend_side = "right")
   ```

   **What happens internally**:

   * **RNA**: `AnnotatedHeatmapFromMAE(mae_rna_meth, assayName="expr", groupColumn="C1A.C1B", sampleType="pathologic_stage", …)` → heatmap of RNA logFC.
   * **Meth**: same with `assayName="beta"` → heatmap of methylation logFC.
   * Both `Heatmap` objects combined side-by-side (`Heatmap1 + Heatmap2`) → returned as `HeatmapList`.

---

<a name="mm_example4c"></a>

#### Example 4C: Three Modalities, Mixed Volcano + Heatmap

**Goal:**

1. Build a volcano from an MAE (e.g. miniACC).
2. Build a heatmap from a fold-change + p-value matrix.
3. Combine both into one figure: left = volcano, right = heatmap.

> **Trick:** We can call `MultiModalPlot` with a mixture of inputs if we specify `panel_type="volcano"` or `"heatmap"` for all. To mix volcano + heatmap, call them separately and combine manually with `patchwork` or `ComplexHeatmap`.

1. **Volcano Panel (miniACC MAE)**

   ```r
   library(MultiModalGraphics)
   library(MultiAssayExperiment)
   library(BiocParallel)

   data("miniACC", package = "MultiAssayExperiment")
   BPPARAM <- MulticoreParam(workers = 2)

   volcano_panel <- MultiModalPlot(
     inputs        = list(ACC = miniACC),
     assayNames    = c(ACC = "RNASeq2GeneNorm"),
     groupColumns  = c(ACC = "C1A.C1B"),
     sampleTypes   = c(ACC = "pathologic_stage"),
     timepoints    = c(ACC = "MethyLevel"),
     panel_type    = "volcano",
     var_quantile  = 0.75,
     pvalue_cutoff = 0.05,
     fc_cutoff     = 0.585,
     max_features  = 5000,
     parallel      = TRUE,
     BPPARAM       = BPPARAM,
     title         = "miniACC: Volcano"
   )
   ```

2. **Heatmap Panel (Simulated Fold-Change + P-Value)**

   ```r
   # Re-use matrices from Example 2B
   logFC_mat <- matrix(rnorm(32, sd = 1.5), nrow = 8, ncol = 4)
   pval_mat  <- matrix(runif(32, 0, 0.2),   nrow = 8, ncol = 4)
   rownames(logFC_mat) <- paste0("Gene", 1:8)
   colnames(logFC_mat) <- paste0("Cond", 1:4)
   rownames(pval_mat)  <- rownames(logFC_mat)
   colnames(pval_mat)  <- colnames(logFC_mat)

   ih_mat <- AnnotatedHeatmapFromMAT(
     logFC_matrix     = logFC_mat,
     pvalue_matrix    = pval_mat,
     pvalue_cutoff    = 0.05,
     trending_cutoff  = 0.1,
     pch_val          = 16,
     unit_val         = 3,
     significant_color= "darkred",
     trending_color   = "darkorange",
     col              = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
     cluster_rows     = TRUE,
     cluster_columns  = TRUE,
     show_row_names   = TRUE,
     show_column_names= TRUE,
     name             = "logFC"
   )
   heatmap_panel <- getHeatmapObject(ih_mat)
   ```

3. **Combine Volcano & Heatmap**

   * **Approach A (Patchwork + gridGraphics)**:
     Draw volcano `ggplot` on left, draw heatmap on right. Requires converting `ggplot` to a grob.

   ```r
   library(patchwork)
   library(ggplot2)
   library(grid)

   # Convert volcano (ggplot) to grob
   volcano_grob <- ggplotify::as.grob(volcano_panel)

   # Draw heatmap to a new page as grob
   heatmap_grob <- grid.grabExpr(ComplexHeatmap::draw(heatmap_panel, heatmap_legend_side="right"))

   # Combine side-by-side
   grid.newpage()
   pushViewport(viewport(layout = grid.layout(1, 2, widths = unit.c(unit(0.5, "npc"), unit(0.5, "npc")))))
   print(volcano_grob, vp = viewport(layout.pos.row=1, layout.pos.col=1))
   grid.draw(heatmap_grob)  # vp automatically next region
   ```

   * **Approach B (GridGraphics)**:
     Use `grid.newpage()` and custom `viewport()`s to place each.

   **Result:**

   * Left half of window: faceted volcano (miniACC).
   * Right half: custom heatmap of simulated FC.

---

<a name="mm_advanced"></a>

### 7.4. Advanced Customization

1. **Per-Modality Color Palettes**
   You can pass modality-specific aesthetics by naming them:

   ```r
   MultiModalPlot(
     inputs        = list(RNA = mae1, Meth = mae2, Prot = mae3),
     assayNames    = c(RNA = "RNASeq2GeneNorm", Meth = "MethylAssay", Prot = "ProtAssay"),
     groupColumns  = c(RNA = "TumorVsNormal", Meth = "TumorVsNormal", Prot = "TumorVsNormal"),
     sampleTypes   = c(RNA = "Stage", Meth = "Stage", Prot = "Stage"),
     panel_type    = "heatmap",
     pvalue_cutoff = 0.01,
     fc_cutoff     = 1,
     max_features  = 200,
     # Per-modality color 'col_RNA', 'col_Meth', 'col_Prot'
     col_RNA       = circlize::colorRamp2(c(-3,0,3), c("#1B9E77","white","#D95F02")),
     col_Meth      = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
     col_Prot      = circlize::colorRamp2(c(-1,0,1), c("purple","white","green")),
     cluster_rows  = list(RNA = TRUE, Meth = FALSE, Prot = TRUE),
     cluster_columns = TRUE,
     show_row_names  = FALSE,
     show_column_names = FALSE
   )
   ```

2. **Mixed-Modality Plots**
   You cannot mix volcano + heatmap within a single `MultiModalPlot` call. Instead, build them separately via two calls and combine manually with `patchwork` or `grid`.

3. **Custom DE Thresholds**
   Per modality, pass named thresholds:

   ```r
   MultiModalPlot(
     inputs        = list(RNA = mae1, Meth = mae2),
     groupColumns  = c(RNA="Contrast1", Meth="Contrast2"),
     fc_cutoff     = c(RNA=1.0, Meth=0.5),
     pvalue_cutoff = c(RNA=0.01, Meth=0.05),
     …
   )
   ```

4. **Parallel vs Sequential Execution**

   * `parallel = TRUE` uses `BiocParallel::bplapply()` if `vectorized=="auto"` and #cells > #workers. To force sequential, use `parallel = FALSE`.
   * For large multi-omic MAEs, consider `BPPARAM = MulticoreParam(workers = 4)` or `SnowParam()` on cluster.

---

<a name="mm_troubleshoot"></a>

### 7.5. Troubleshooting & Tips

1. **“`panel_type='heatmap'` but DE Table Provided”**

   * *Symptom*:

     ```
     Error: DE tables only support `panel_type='volcano'`
     ```
   * *Fix*: Use `panel_type = "volcano"` if your input is precomputed DE. To heatmap that, use `AnnotatedHeatmapFromMAT()` instead.

2. **Missing Modality Mismatch**

   * *Symptom*:

     ```
     For matrix modality 'XYZ', please supply `groupColumns[['XYZ']]`.
     ```
   * *Fix*: Provide `groupColumns` entry for each named modality that is a matrix or MAE.

3. **One Volcano + One Heatmap in Single Call**

   * *Not Supported*: You must call `MultiModalPlot` separately in “volcano” and “heatmap” modes, then combine manually.

4. **“Insufficient Samples”**

   * If any cell (SampleType × timePoint) has < `min_samples` or only one level of `groupColumn`, that cell is skipped. If *all* cells skip, you get `stop("No DE results to plot.")`.
   * *Fix*: Lower `min_samples`, remove `timepoint`, or adjust grouping.

5. **Large MultiModal Arrays**

   * For >3 modalities, heatmap combining can become unreadable. Consider splitting into separate figures or using scrolling (e.g., `ComplexHeatmap::draw(..., heatmap_height = unit(n, "cm"))`).

---

<a name="session"></a>

## 8. Session Information & Citation

For reproducibility, here is a sample `sessionInfo()`.

```r
sessionInfo()
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.2.3

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] MultiModalGraphics_0.99.5 ggplot2_3.4.0            dplyr_1.0.10             
# [4] magrittr_2.0.3            matrixStats_0.63.1        BiocParallel_1.30.0      
# [7] limma_3.52.0              SummarizedExperiment_1.26.0 MultiAssayExperiment_2.8.0
# [10] ComplexHeatmap_2.6.0      circlize_0.4.0            RColorBrewer_1.1-3       

# loaded via a namespace (and not attached):
# [1] BiocGenerics_0.36.0 XVector_0.36.0      IRanges_2.32.0      S4Vectors_0.36.0   
# [5] Biobase_2.56.0      lifecycle_1.0.1     grid_4.2.0          gtable_0.3.0       
# [9] permute_0.9-5       farver_2.1.0        usethis_2.1.6       zoo_1.8-9          
# [13] foreach_1.5.2       cluster_2.1.4       doParallel_1.0.17   codetools_0.2-18   
# [17] ...  
```

If you use **MultiModalGraphics** in a publication, please cite:

```
Mohammed FA, Fall M, Hammamieh R, Muhie S (2025). MultiModalGraphics: Multi‐Modal Visualization & Plotting Utilities. R package version 0.99.5.
```

Visit the GitHub repository for more examples, bug reports, and feature requests:

```
https://github.com/famanalytics0/MultiModalGraphics
```

---

**End of MultiModalGraphics User Manual**
