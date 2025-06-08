# MultiModalGraphics

**Version:** 0.99.5  
**Authors:** Foziya Ahmed Mohammed, Malick Fall, Rasha Hammamieh, Seid Muhie  
**License:** Artistic-2.0  

A unified R/Bioconductor package for creating publication-quality volcano plots and custom heatmaps from diverse biological data sources. Designed for maximum flexibility and ease of use, the package supports:

- **Faceted volcano plots** (log₂ fold-change vs. –log₁₀(p-value)) from:
  - Precomputed DE tables
  - Expression/count matrices + metadata
  - `MultiAssayExperiment` (MAE) objects
- **Informative heatmaps** with overlaid significance markers from:
  - Fold-change & p-value matrices
  - Matrix + metadata DE pipelines
  - MAE + optional iClusterPlus DE pipelines
- **Multifeature grids**: 2D tile-and-point heatmaps (z-scores as tile fills, p-values as point color, counts as point size)

This README guides you through installation, per-class usage (ClearScatterplot, InformativeHeatmap, MultifeatureGrid), and advanced tips. Examples range from “one-line” commands for novices to full customization for power users.

---


## Quick Start

Welcome to **MultiModalGraphics**. This section gives you a rapid, hands-on introduction to installing the package and producing your first plots. For full details, see the [User Manual](MANUAL.md).

---

### 1. Installation

```r
# (1) Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# (2) Install MultiModalGraphics from Bioconductor
BiocManager::install("MultiModalGraphics")

# (3) Load the package
library(MultiModalGraphics)

# (4) (Optional) Install suggested data packages for examples:
BiocManager::install(c("curatedPCaData", "curatedTCGAData", "curatedTCGADataData"))
````

---

### 2. Core Dependencies

```r
library(ggplot2)             # for plotting
library(dplyr)               # for data handling
library(matrixStats)         # for variance filtering
library(BiocParallel)        # for parallel DE
library(MultiAssayExperiment) 
library(SummarizedExperiment)
library(limma)               # for differential expression
library(RColorBrewer)        # for palettes
library(ComplexHeatmap)      # for heatmaps
```

---

### 3. Volcano Plots with **ClearScatterplot**

#### 3.1. From a Precomputed DE Table

```r
# (a) Simulate a small DE table
set.seed(101)
df_de <- data.frame(
  log2fc     = rnorm(200, sd=1),
  negLog10p  = -log10(runif(200, 0, 0.1)),
  regulation = sample(c("up","down"), 200, replace=TRUE),
  SampleType = rep(c("Group1","Group2"), each=100),
  timePoint  = rep(c("T1","T2"), times=100),
  stringsAsFactors = FALSE
)

# (b) Construct a ClearScatterplot object
cs1 <- ClearScatterplot(
  data           = df_de,
  highLog2fc     = 1.0,      # log2FC > 1 → “up”
  lowLog2fc      = -1.0,     # log2FC < -1 → “down”
  negLog10pValue = 1.301     # –log10(p) > 1.301 → p < 0.05
)

# (c) Build and display the ggplot
cs1 <- createPlot(
  cs1,
  color1 = "tomato",      # color for “down”
  color2 = "grey50",         # color for “neutral”
  color3 = "steelblue",         # color for “up”
  title  = "Example Volcano (Precomputed DE)",
  subtitle = "Faceted by SampleType (X) and timePoint (Y)"
)
show(cs1)
```

# Tiny smoke‐test dataset

```r
set.seed(42)
df <- data.frame(
  log2fc     = rnorm(100),
  negLog10p  = runif(100, 0, 5),
  regulation = sample(c("up","down"), 100, TRUE),
  SampleType = sample(c("A","B"), 100, TRUE),
  stringsAsFactors = FALSE
)

# This line must appear first, or you’ll see “object 'cs' not found”
cs <- ClearScatterplot(df)

# Now build + draw in one step
show(
  createPlot(
    cs,
    color1          = "tomato",     # up = steelblue
    color2          = "lightgrey",     # neutral = lightgrey
    color3          = "steelblue",        # down = tomato
    xlab            = "Custom X-Axis",
    ylab            = "Custom Y-Axis",
    text_family     = "Times New Roman",
    text_size       = 14,
    point_size      = 3,
    point_alpha     = 0.8,
    legend_position = "right",
    legend_title    = "My Legend",
    legend_labels   = c("Down","Neutral","Up"),
    custom_theme    = theme_minimal(),
    title           = "Volcano Plot—No Hard‐Coded Numeric Flags"
  )
)
```

> **Result:** A faceted volcano plot with up/down/neutral points colored as specified, and counts of “up” and “down” annotated in each facet.

---

#### 3.2. From an Expression Matrix + Metadata

```r
# (a) Simulate a small counts matrix: 100 genes × 20 samples
set.seed(102)
expr <- matrix(rpois(100*20, lambda=10), nrow=100, ncol=20)
rownames(expr) <- paste0("Gene", 1:100)
colnames(expr) <- paste0("Sample", 1:20)

# (b) Simulate sample metadata
meta <- data.frame(
  Group      = rep(c("A","B"), each=10),
  SampleType = rep(c("X","Y"), times=10),
  timePoint  = rep(c("T1","T2"), each=5, times=2),
  row.names  = paste0("Sample", 1:20),
  stringsAsFactors = FALSE
)

# (c) Create a ClearScatterplot_table object (runs limma DE under the hood)
cs2 <- ClearScatterplot_table(
  expr         = expr,
  meta         = meta,
  groupColumn  = "Group",
  sampleType   = "SampleType",
  timepoint    = "timePoint",
  dataType     = "auto",       # automatically detects “count”
  var_quantile = 0.75,         # filters out lowest 25% variance genes
  parallel     = TRUE,
  BPPARAM      = MulticoreParam(workers = 2)
)

# (d) Build and display the volcano
cs2 <- createPlot(
  cs2,
  color1 = "navy", color2 = "grey80", color3 = "firebrick",
  title  = "Matrix + Metadata Volcano"
)
show(cs2)
```

> **Result:** The package performs variance filtering, runs voom+limma per SampleType–timePoint cell, then stitches all DE results into a single faceted volcano plot.

---

#### 3.3. From a **MultiAssayExperiment** Object

```r
# (a) Load a small example MAE (miniACC)
data("miniACC", package="MultiAssayExperiment")

# (b) Set up parallel workers
BPPARAM <- MulticoreParam(workers = 2)

# (c) Create a ClearScatterplot_MAE object
cs3 <- ClearScatterplot_MAE(
  mae            = miniACC,
  assayName      = "RNASeq2GeneNorm",
  groupColumn    = "C1A.C1B",
  sampleType     = "pathologic_stage",
  timepoint      = "MethyLevel",
  dataType       = "auto",
  vectorized     = "auto",
  parallel       = TRUE,
  BPPARAM        = BPPARAM,
  var_quantile   = 0.75
)

# (d) Build and display
cs3 <- createPlot(
  cs3,
  color1 = "darkblue",
  color2 = "grey70",
  color3 = "darkred",
  title  = "miniACC Volcano: C1A vs C1B"
)
show(cs3)
```

> **Result:** DE is run per facet (Stage × Methylation) using the chosen assay in the MAE. The final output is a publication-ready faceted volcano.

---

### 4. Custom Heatmaps with **InformativeHeatmap**

#### 4.1. From a Fold-Change + P-Value Matrix

```r
# (a) Simulate logFC and p-value matrices
set.seed(103)
logFC_mat <- matrix(rnorm(50*4, mean=0, sd=1), nrow=50, ncol=4)
pval_mat  <- matrix(runif(50*4, 0, 0.2), nrow=50, ncol=4)
rownames(logFC_mat) <- paste0("Gene", 1:50)
colnames(logFC_mat) <- paste0("Cond", 1:4)
rownames(pval_mat)  <- rownames(logFC_mat)
colnames(pval_mat)  <- colnames(logFC_mat)

# (b) Build InformativeHeatmap directly from matrices
ih1 <- InformativeHeatmapFromMAT(
  logFC_matrix     = logFC_mat,
  pvalue_matrix    = pval_mat,
  pvalue_cutoff    = 0.05,
  trending_cutoff  = 0.1,
  pch_val          = 16,
  unit_val         = 4,
  significant_color= "black",
  trending_color   = "yellow",
  col              = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = FALSE,
  show_column_names= TRUE,
  name             = "log2FC"
)

# (c) Extract and draw
ht1 <- getHeatmapObject(ih1)
ComplexHeatmap::draw(ht1, heatmap_legend_side = "right")
```

> **Result:** A 50×4 heatmap of log₂FC with colored points overlayed where p < 0.05 (black) or 0.05 ≤ p < 0.1 (yellow).

---

#### 4.2. From a Matrix + Metadata (DE‐Based Heatmap)

```r
# (a) Reuse expr (100×20) and meta from Section 3.2
# (b) Create InformativeHeatmap (performs DE per SampleType as columns)
ih2 <- InformativeHeatmap(
  data             = expr,
  meta             = meta,
  groupColumn      = "Group",
  sampleType       = "SampleType",
  timepoint        = NULL,
  dataType         = "continuous",
  var_quantile     = 0.5,
  pvalue_cutoff    = 0.05,
  trending_cutoff  = 0.1,
  fc_cutoff        = 0.585,
  min_samples      = 3,
  max_features     = 30,
  runClustering    = FALSE,
  pch_val          = 16,
  unit_val         = 3,
  significant_color= "darkred",
  trending_color   = "orange",
  heatmap_scale    = "logFC",
  col              = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = FALSE,
  show_column_names= TRUE,
  name             = "log2FC"
)

# (c) Extract and draw
ht2 <- getHeatmapObject(ih2)
ComplexHeatmap::draw(ht2, heatmap_legend_side = "right")
```

> **Result:** DE is run per SampleType, the top 30 features by variance are displayed as log₂FC, and significant/trending points are overlaid.

---

#### 4.3. From a **MultiAssayExperiment** + **iClusterPlus**

```r
# (a) Assume you have a small MAE named `tiny_mae` with two assays
#     (see the Manual for how to build tiny_mae from curatedTCGAData)

# (b) Build InformativeHeatmapFromMAE (runs iClusterPlus → limma)
ih3 <- InformativeHeatmapFromMAE(
  mae                = tiny_mae,
  significant_pvalue = 0.05,
  trending_pvalue    = 0.1,
  pch_val            = 16,
  unit_val           = 3,
  significant_color  = "darkred",
  trending_color     = "darkorange",
  K                  = 2,
  lambda             = 0.1,
  coef               = 2,
  heatmap_scale      = "logFC",
  col                = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
  cluster_rows       = TRUE,
  cluster_columns    = TRUE,
  show_row_names     = FALSE,
  show_column_names  = TRUE,
  name               = "log2FC"
)

# (c) Extract and draw
ht3 <- getHeatmapObject(ih3)
ComplexHeatmap::draw(ht3, heatmap_legend_side = "right")
```

> **Result:** Samples are clustered via iClusterPlus, DE is run per cluster vs. cluster, and the combined log₂FC matrix is shown with significant/trending overlays.

---

### 5. 2D Tile-and-Point Heatmaps with **MultifeatureGrid**

#### 5.1. Using Built-In Demo Data

```r
# (a) Fetch the demo data frame
demo_df <- get_multifeature_grid_df()
str(demo_df)
# 'data.frame': 16 rows × 6 columns:
#  $ tissue             : Factor
#  $ signaling          : Factor
#  $ Activation_z_score : num
#  $ p                  : num
#  $ number_of_genes    : int
#  $ timePoint          : Factor

# (b) Construct the MultifeatureGrid object
mg <- MultifeatureGrid(
  data           = demo_df,
  title          = "Demo Multifeature Grid",
  x_label        = "Tissue",
  y_label        = "Pathway",
  logpval_label  = "-Log10(P-Value)",
  zscore_label   = "Activation Z-Score",
  numitems_label = "Gene Count",
  color_palette  = "RdYlBu",
  breaks         = seq(-2, 2, 0.5)
)

# (c) Plot the heatmap
plot_heatmap(
  mg,
  pValueColumn       = "p",
  lowColor           = "lightblue",
  highColor          = "darkblue",
  borderColor        = "grey50",
  columnForNumber    = "number_of_genes",
  independantVariable= "timePoint"
)
```

> **Result:** A tile-and-point grid: tiles show z-scores, colored points show –log₁₀(p), and point size reflects gene count, faceted by `timePoint`.

---

### 6. Unified Interface with `MultiModalPlot`

When you have **multiple modalities** to display (e.g., both a volcano and a heatmap side by side, or two different assays), you can call a single function, **`MultiModalPlot()`**, which auto-detects input types and stitches panels together.

```r
# Example: Volcano (MAE) + Heatmap (precomputed)
volcano_panel <- ClearScatterplot_MAE(
  mae          = miniACC,
  assayName    = "RNASeq2GeneNorm",
  groupColumn  = "C1A.C1B",
  sampleType   = "pathologic_stage",
  timepoint    = "MethyLevel",
  parallel     = FALSE
)

# Simulate simple logFC + p-value matrices
set.seed(104)
lfc  <- matrix(rnorm(20*3, sd=1), nrow=20, ncol=3)
pmat <- matrix(runif(20*3, 0, 0.2), nrow=20, ncol=3)
rownames(lfc)  <- paste0("Gene", 1:20)
rownames(pmat) <- rownames(lfc)
colnames(lfc)  <- paste0("Cond", 1:3)
colnames(pmat) <- colnames(lfc)

# Wrap in a named list
inputs <- list(
  Volcano = miniACC,
  Heatmap = list(logFC = lfc, pval = pmat)
)

# Call MultiModalPlot()
mm_plot <- MultiModalPlot(
  inputs,
  assayNames    = c(Volcano = "RNASeq2GeneNorm"),
  groupColumns  = c(Volcano = "C1A.C1B"),
  sampleTypes   = c(Volcano = "pathologic_stage"),
  timepoints    = c(Volcano = "MethyLevel"),
  panel_type    = "volcano",       # “volcano” for volcano panels
  color1        = "navy",
  color3        = "firebrick"
)

print(mm_plot)

# To request a heatmap instead:
mm_heat <- MultiModalPlot(
  inputs,
  assayNames    = c(Volcano = "RNASeq2GeneNorm"),
  groupColumns  = c(Volcano = "C1A.C1B"),
  sampleTypes   = c(Volcano = "pathologic_stage"),
  timepoints    = c(Volcano = "MethyLevel"),
  panel_type    = "heatmap",
  col_Volcano   = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
  cluster_rows  = TRUE,
  cluster_columns = TRUE
)

ComplexHeatmap::draw(mm_heat, heatmap_legend_side = "right")
```

> **Note:** When `MultiModalPlot()` sees a `MultiAssayExperiment` it dispatches to `ClearScatterplot_MAE` or `InformativeHeatmapFromMAE`. If it sees a **list** with `expr` & `meta`, it dispatches to “table” constructors. If it sees precomputed matrices or DE tables, it uses the direct-from-mat constructors.

---

### 7. Linking the Manual

For **in-depth guidance**, step-by-step walkthroughs, and advanced customization, please consult the **User Manual**:

**→ [MANUAL.md](MANUAL.md) ←**

The Manual provides:

* Detailed explanations of every function and argument
* Real-world, lightweight examples (TCGA, prostate cancer, miniACC)
* Advanced tips for parallelization, custom layer functions, color theming, and troubleshooting

---

With these few lines of code, you can generate publication-ready volcano plots, custom heatmaps, or integrated multimodal displays. For more examples and explanations, see the User Manual linked above.

---
## Table of Contents

1. [Installation](#installation)  
2. [Overview](#overview)  
3. [ClearScatterplot (Volcano Plots)](#clearscatterplot)  
   1. [What It Does](#cs-what)  
   2. [Quick Start](#cs-quickstart)  
   3. [API Details](#cs-api)  
   4. [Examples](#cs-examples)  
   5. [Advanced Customization & Tips](#cs-advanced)  
4. [InformativeHeatmap (Custom Heatmaps)](#informativeheatmap)  
   1. [What It Does](#ih-what)  
   2. [Quick Start](#ih-quickstart)  
   3. [API Details](#ih-api)  
   4. [Examples](#ih-examples)  
   5. [Advanced Customization & Tips](#ih-advanced)  
5. [MultifeatureGrid (2D Tile-and-Point Heatmaps)](#multifeaturegrid)  
   1. [What It Does](#mg-what)  
   2. [Quick Start](#mg-quickstart)  
   3. [API Details](#mg-api)  
   4. [Examples](#mg-examples)  
   5. [Advanced Customization & Tips](#mg-advanced)  
6. [Unified Multimodal API: `MultiModalPlot`](#multimodal)  
7. [Session Info & Citation](#session)  

---

<a name="installation"></a>
## 1. Installation

### 1.1. Prerequisites

- **R ≥ 4.1.0**  
- **Bioconductor ≥ 3.15**  

### 1.2. Install from Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MultiModalGraphics")
````

### 1.3. Install Development Version from GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_git(
  "https://github.com/famanalytics0/MultiModalGraphics.git",
  ref         = "master",
  dependencies= TRUE,
  upgrade     = "never"
)
```

### 1.4. Load Core Dependencies

```r
library(MultiModalGraphics)
library(ggplot2)            # for plotting
library(dplyr)              # for data manipulation
library(matrixStats)        # for variance filtering
library(BiocParallel)       # for parallel DE
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(limma)
library(RColorBrewer)       # for color palettes
library(ComplexHeatmap)     # for heatmap rendering
```

> **Tip:** Data-package examples require:
>
> ```r
> BiocManager::install(c("curatedPCaData","curatedTCGAData","curatedTCGADataData"))
> ```

---

<a name="overview"></a>

## 2. Overview

MultiModalGraphics centers on **three S4 classes**:

1. **ClearScatterplot**
   Faceted volcano plots of DE results.
   Constructors:

   * `ClearScatterplot(data.frame)`: precomputed DE
   * `ClearScatterplot_table(expr_matrix, meta_data)`: DE from raw data
   * `ClearScatterplot_MAE(mae_object, assayName, …)`: DE from MAE
     Plotting:
   * `createPlot(ClearScatterplot_obj, …)` → builds ggplot
   * `show(ClearScatterplot_obj)` → prints plot

2. **InformativeHeatmap**
   Customized ComplexHeatmap with overlaid significance/trending markers.
   Constructors:

   * `InformativeHeatmap(expr_matrix, meta_data, …)`: DE on raw data
   * `InformativeHeatmapFromMAT(logFC_matrix, pval_matrix, …)`: precomputed matrices
   * `InformativeHeatmapFromMAE(mae_object, …)`: DE + iClusterPlus from MAE
     Utilities:
   * `updateLayerFun(InformativeHeatmap_obj, layer_fun)`
   * `getHeatmapObject(InformativeHeatmap_obj)`

3. **MultifeatureGrid**
   2D tile-and-point heatmap: tile fill = z-score, point color = –log₁₀(p), point size = count.
   Constructor:

   * `MultifeatureGrid(data_frame, …)`
     Plotting:
   * `plot_heatmap(MultifeatureGrid_obj, …)`

Additionally, **`MultiModalPlot`** unifies all three into a single high-level API, auto-dispatching by input type.

---

<a name="clearscatterplot"></a>

## 3. ClearScatterplot (Volcano Plots)

<a name="cs-what"></a>

### 3.1. What It Does

Creates faceted volcano plots (log₂ fold-change vs. –log₁₀(p-value)):

* **“Up”** features: log₂FC > `highLog2fc` & –log₁₀(p) > `negLog10pValue`
* **“Down”** features: log₂FC < `lowLog2fc` & –log₁₀(p) > `negLog10pValue`
* **“Neutral”**: everything else

Supports three input modes:

1. **Precomputed DE results** (`data.frame` with `log2fc`, `negLog10p`, `regulation`, `SampleType`, optional `timePoint`)
2. **Raw expression/count matrices + metadata** (runs `limma` or `limma::voom`)
3. **`MultiAssayExperiment` objects** (extract assay, optional merge with top-level metadata, runs DE per facet)

Key features:

* **Automatic variance filtering** (dropping low-variance features)
* **Optional parallelization** via `BiocParallel`
* **Faceting** by `SampleType` (columns) and optional `timePoint` (rows)
* **Count annotations**: displays number of up/down features per facet

---

<a name="cs-quickstart"></a>

### 3.2. Quick Start

```r
library(MultiModalGraphics)
library(MultiAssayExperiment)
library(BiocParallel)

# 1. Load example MAE: miniACC (Bioc data package)
data("miniACC", package = "MultiAssayExperiment")

# 2. Choose 2 cores for parallel DE
BPPARAM <- MulticoreParam(workers = 2)

# 3. Create volcano from MAE
cs_mae <- ClearScatterplot_MAE(
  mae            = miniACC,
  assayName      = "RNASeq2GeneNorm",
  groupColumn    = "C1A.C1B",         # “C1A” vs “C1B”
  sampleType     = "pathologic_stage",# facet columns
  timepoint      = "MethyLevel",      # facet rows
  dataType       = "auto",
  vectorized     = "auto",
  parallel       = TRUE,
  BPPARAM        = BPPARAM,
  var_quantile   = 0.75
)

# 4. Build & display the ggplot
cs_mae <- createPlot(
  cs_mae,
  color1 = "cornflowerblue",   # down
  color2 = "grey",             # neutral
  color3 = "indianred",        # up
  title  = "miniACC Volcano: C1A vs C1B",
  subtitle = "Faceted by Stage (X) & Methylation (Y)",
  caption  = "Data: miniACC"
)
show(cs_mae)
```

> **One-line** version using `MultiModalPlot`:
>
> ```r
> volcano_plot <- MultiModalPlot(
>   inputs        = list(ACC = miniACC),
>   assayNames    = c(ACC = "RNASeq2GeneNorm"),
>   groupColumns  = c(ACC = "C1A.C1B"),
>   sampleTypes   = c(ACC = "pathologic_stage"),
>   timepoints    = c(ACC = "MethyLevel"),
>   panel_type    = "volcano",
>   parallel      = TRUE,
>   BPPARAM       = BPPARAM,
>   title         = "miniACC Volcano"
> )
> print(volcano_plot)
> ```

---

<a name="cs-api"></a>

### 3.3. API Details

#### 3.3.1. `ClearScatterplot(data, highLog2fc=0.585, lowLog2fc=-0.585, negLog10pValue=1.301)`

* **Arguments**

  * `data`: `data.frame` with columns

    * `log2fc` (numeric)
    * `negLog10p` (numeric = –log₁₀(p))
    * `regulation` (character: `"up"`/`"down"`)
    * `SampleType` (factor/character)
    * Optional: `timePoint` (factor/character) for 2-D faceting
  * `highLog2fc`, `lowLog2fc`: thresholds for “up”/“down”
  * `negLog10pValue`: –log₁₀(p) threshold for significance

* **Returns**

  * An S4 **`ClearScatterplot`** object with:

    * `@data`: filtered `data.frame` plus computed `color_flag` (–1,0,1)
    * `@plot`: initially `NULL`

* **Validity Checks**

  * Must have required columns; stops if missing.
  * Filters out rows with `NA` in `log2fc` or `negLog10p` (warnings).

#### 3.3.2. `ClearScatterplot_table(expr, meta, groupColumn="Group", sampleType="SampleType", timepoint=NULL, dataType=c("auto","continuous","count"), vectorized=c("auto","perCell","vectorized"), parallel=TRUE, BPPARAM, var_quantile=0.75)`

* **Arguments**

  * `expr`: numeric matrix (features × samples)
  * `meta`: `data.frame` (samples × metadata), rownames = colnames(expr)
  * `groupColumn`: column in `meta` for DE contrast (must have ≥ 3 samples per group)
  * `sampleType`: column in `meta` for X-facet (will create one volcano per unique level)
  * `timepoint`: optional column in `meta` for Y-facet (optional 2D faceting)
  * `dataType`: `"auto"` (infer continuous vs count) or `"continuous"` or `"count"`
  * `vectorized`: `"auto"`, `"perCell"`, or `"vectorized"` (controls parallel logic)
  * `parallel`: `TRUE`/`FALSE`. If `FALSE`, always run sequentially.
  * `BPPARAM`: `BiocParallelParam` backend for parallel DE
  * `var_quantile`: \[0–1] quantile for filtering out low-variance features

* **Workflow**

  1. Filter out features with row-variance < `quantile(rv, var_quantile)`.
  2. Drop samples with `NA` in `groupColumn`, `sampleType`, or `timepoint`.
  3. Ensure each group has ≥ 3 samples; drop cells (SampleType × timePoint) with < 2 groups.
  4. For each cell, run:

     * If `dataType=="count"`, `voom()` → `lmFit()`→`eBayes()`.
     * Else, optionally log₂-transform (if large) → `lmFit()`→`eBayes()`.
     * `topTable()` → build `data.frame(log2fc, negLog10p, regulation, SampleType, timePoint)`.
  5. Combine all cell results; throw error if empty.
  6. Pass combined DE table to `ClearScatterplot()`.

* **Returns**

  * A **`ClearScatterplot`** object ready for `createPlot()`.

#### 3.3.3. `ClearScatterplot_MAE(mae, assayName, groupColumn="Group", sampleType="SampleType", timepoint=NULL, dataType=c("auto","continuous","count"), vectorized=c("auto","perCell","vectorized"), parallel=TRUE, BPPARAM, var_quantile=0.75)`

* **Arguments**

  * `mae`: a `MultiAssayExperiment` object
  * `assayName`: name of assay to extract (must exist in `experiments(mae)`)
  * `groupColumn`, `sampleType`, `timepoint`: metadata columns to use
  * `dataType`, `vectorized`, `parallel`, `BPPARAM`, `var_quantile`: as above

* **Workflow**

  1. Extract `SummarizedExperiment` via `experiments(mae)[[assayName]]`.
  2. Filter features by variance (top `1 − var_quantile`).
  3. Assemble assay-level `meta <- colData(se)`. If any of `groupColumn`, `sampleType`, or `timepoint` not present in `meta`, merge with top-level `colData(mae)` via `sampleMap(mae)`.
  4. Intersect sample names; drop unmatched.
  5. Drop samples with `NA` in required metadata.
  6. Enforce at least 3 samples per group, total ≥ 6 samples.
  7. Internally call the same core DE logic as `ClearScatterplot_table` (via `.ClearScatterplot_core`).

* **Returns**

  * A **`ClearScatterplot`** object ready for `createPlot()`.

#### 3.3.4. `createPlot(object, color1="cornflowerblue", color2="grey", color3="indianred", xlab=expression(log2~fold~change), ylab=expression(-log10~p), custom_theme=NULL, point_alpha=0.5, point_size=1.75, legend_position="bottom", legend_title=NULL, legend_labels=NULL, text_family="sans", text_size=10, ...)`

* **Arguments**

  * `object`: a `ClearScatterplot` instance
  * `color1`, `color2`, `color3`: colors for “down” (−1), “neutral” (0), “up” (+1) points
  * `xlab`, `ylab`: axis labels (expressions or strings)
  * `custom_theme`: extra `ggplot2` theme to overlay
  * `point_alpha`, `point_size`: aesthetics for scatter points
  * `legend_position`: `"bottom"`, `"top"`, `"left"`, or `"right"`
  * `legend_title`, `legend_labels`: legend title and labels (length 3: down, neutral, up)
  * `text_family`, `text_size`: font family & base size
  * `…`: additional `ggplot2::labs()` arguments (e.g., `title`, `subtitle`, `caption`)

* **Behavior**

  1. Checks if `@plot` is `NULL`; if so, builds a `ggplot` with:

     * Jittered + exact points (colored by `color_flag`)
     * Facets: `timePoint ~ SampleType` if `timePoint` exists; else `. ~ SampleType`
     * Custom `theme_bw()` base + text styling
     * Adds count labels (number of up/down) in facet corners

  2. Stores the ggplot in `object@plot` and returns invisibly.

* **Usage**

  ```r
  cs <- createPlot(
    cs,
    color1 = "navy", color2 = "grey50", color3 = "firebrick",
    title   = "Volcano Plot",
    point_size = 2
  )
  show(cs)  # prints the plot
  ```

#### 3.3.5. `show(object)`

* **Behavior**

  * If `object@plot` is `NULL`, calls `createPlot(object)` with defaults.
  * Prints the `ggplot`.
  * Invisibly returns `object`.

---

<a name="cs-examples"></a>

### 3.4. Examples

#### 3.4.1. Example A: Precomputed DE Table

```r
# Simulate a small DE table for “Basal-like” vs “Luminal A” in TCGA-BRCA
set.seed(202)
df_de <- data.frame(
  log2fc     = rnorm(100, sd = 1),
  negLog10p  = -log10(runif(100, 0, 0.1)),
  regulation = sample(c("up","down"), 100, replace = TRUE),
  SampleType = "Basal-vs-LuminalA",
  timePoint  = NA_character_,
  stringsAsFactors = FALSE
)

# Construct and plot
cs1 <- ClearScatterplot(df_de, highLog2fc = 1.0, lowLog2fc = -1.0, negLog10pValue = 1.301)
cs1 <- createPlot(
  cs1,
  color1 = "steelblue",     # down
  color3 = "tomato",        # up
  legend_labels = c("Down","Neutral","Up"),
  title = "Simulated TCGA-BRCA Volcano"
)
show(cs1)
```

#### 3.4.2. Example B: Matrix + Metadata

```r
library(curatedPCaData)   # for real “Taylor” example
data("Taylor", package = "curatedPCaData")

# Extract counts matrix & metadata
expr_taylor <- assay(Taylor, "counts")                       # ~12.6k × 106
meta_taylor <- as.data.frame(colData(Taylor), stringsAsFactors = FALSE)

# Limit to two groups: “Tumor” vs “Normal”
keep_idx <- meta_taylor$DiseaseStatus %in% c("Tumor","Normal")
expr2 <- expr_taylor[, keep_idx]
meta2 <- meta_taylor[keep_idx, , drop = FALSE]
rownames(meta2) <- rownames(meta_taylor)[keep_idx]

# Create ClearScatterplot_table
cs2 <- ClearScatterplot_table(
  expr          = expr2,
  meta          = meta2,
  groupColumn   = "DiseaseStatus",
  sampleType    = "Race",
  timepoint     = "GleasonScore",
  dataType      = "auto",
  var_quantile  = 0.75,
  parallel      = TRUE,
  BPPARAM       = MulticoreParam(workers = 2)
)
cs2 <- createPlot(
  cs2,
  title = "Taylor PCa: Tumor vs Normal",
  subtitle = "Faceted by Race (X) & Gleason Score (Y)",
  color1 = "navy", color3 = "darkred"
)
show(cs2)
```

#### 3.4.3. Example C: MAE-based DE

```r
library(MultiAssayExperiment)
data("miniACC", package = "MultiAssayExperiment")

cs3 <- ClearScatterplot_MAE(
  mae          = miniACC,
  assayName    = "RNASeq2GeneNorm",
  groupColumn  = "C1A.C1B",
  sampleType   = "pathologic_stage",
  timepoint    = "MethyLevel",
  var_quantile = 0.75,
  parallel     = TRUE,
  BPPARAM      = MulticoreParam(workers = 2)
)
cs3 <- createPlot(
  cs3,
  title = "miniACC: C1A vs C1B Volcano",
  color1 = "darkblue", color3 = "firebrick"
)
show(cs3)
```

---

<a name="cs-advanced"></a>

### 3.5. Advanced Customization & Tips

* **Axis Expressions**
  Use mathematical notation:

  ```r
  createPlot(cs, xlab = expression(log[2]~fold~change), ylab = expression(-log[10]~p))
  ```

* **Legend & Color Order**
  By default, `scale_color_manual(values = c(color1, color2, color3))` maps
  –1 → `color1` (“down”), 0 → `color2` (“neutral”), +1 → `color3` (“up”).
  Adjust ordering if needed via `legend_labels`.

* **Faceting Without `timePoint`**
  If your data has no `timePoint` column, faceting defaults to `. ~ SampleType` (single row).

* **Fine-Tuned Parallelization**

  * `vectorized = "auto"`: uses `bplapply()` if #cells > #workers.
  * To force no parallel: `parallel = FALSE` or `vectorized = "perCell"`.

* **Dropping Low-Variance Genes**
  Lower `var_quantile` (e.g., `0.5`) to filter more aggressively. Or set `var_quantile = 1` to skip filtering.

* **Plot in Two Steps**

  ```r
  cs <- ClearScatterplot(df_de)
  cs <- createPlot(cs, …)
  # Modify cs@plot, e.g. add theme changes
  cs@plot <- cs@plot + ggplot2::theme_minimal()
  show(cs)
  ```

* **Common Errors & Fixes**

  1. **Missing Columns**: check `colnames(meta)` for exact names (case sensitive).
  2. **No Overlap** (MAE): ensure `sampleMap(mae)` correctly links assay‐level `colData` and top-level `colData`.
  3. **Too Few Samples**: each group in each facet must have ≥ 3 samples if `min_samples=3`. Adjust group/facet columns or lower threshold.

---

<a name="informativeheatmap"></a>

## 4. InformativeHeatmap (Custom Heatmaps)

<a name="ih-what"></a>

### 4.1. What It Does

Generates a **ComplexHeatmap** of features × conditions, with:

* **Tile fills** = continuous values (expression or log₂ fold-change)
* **Overlaid points** indicating statistical significance:

  * Color = `significant_color` if p < `pvalue_cutoff`
  * Color = `trending_color` if `pvalue_cutoff` ≤ p < `trending_cutoff`
* **Optional clustering** of rows/columns
* **Sensible defaults**: blue-white-red ramp, black/yellow significance dots, etc.

Three workflows:

1. **Raw matrix + metadata**: runs DE per facet, assembles logFC & p-value matrices → heatmap.
2. **Precomputed logFC & p-value matrices**: directly builds heatmap with overlaid points.
3. **MAE + iClusterPlus**: clusters samples, runs DE per assay, combines logFC & p values → heatmap.

Built-in dot overlay is vectorized (no explicit loops), guaranteeing performance even on large matrices.

---

<a name="ih-quickstart"></a>

### 4.2. Quick Start

#### 4.2.1. Matrix + Metadata Mode

```r
library(MultiModalGraphics)
library(matrixStats)
library(limma)
library(BiocParallel)
library(ComplexHeatmap)
library(circlize)

# Simulate small expression matrix: 100 features × 20 samples
set.seed(123)
expr <- matrix(rnorm(100*20, mean=5, sd=2), nrow=100, ncol=20)
rownames(expr) <- paste0("Gene", 1:100)
colnames(expr) <- paste0("S", 1:20)

# Simulate metadata: Group (“A”/“B”), SampleType (“X”/“Y”)
meta <- data.frame(
  Group      = rep(c("A","B"), each=10),
  SampleType = rep(c("X","Y"), times=10),
  stringsAsFactors = FALSE,
  row.names = paste0("S", 1:20)
)

ih_mat <- InformativeHeatmap(
  data             = expr,
  meta             = meta,
  groupColumn      = "Group",
  sampleType       = "SampleType",
  timepoint        = NULL,
  dataType         = "continuous",
  var_quantile     = 0.5,
  pvalue_cutoff    = 0.05,
  trending_cutoff  = 0.1,
  fc_cutoff        = 0.585,
  min_samples      = 3,
  max_features     = 30,
  runClustering    = FALSE,
  pch_val          = 16,
  unit_val         = 3,
  significant_color= "darkred",
  trending_color   = "orange",
  heatmap_scale    = "logFC",
  col              = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = TRUE,
  show_column_names= TRUE,
  name             = "logFC"
)

ht_mat <- getHeatmapObject(ih_mat)
ComplexHeatmap::draw(ht_mat, heatmap_legend_side = "right")
```

#### 4.2.2. Precomputed LogFC + P-Value Matrices

```r
# Simulate logFC (8 genes × 4 comparisons) and p-values
set.seed(456)
logFC_mat <- matrix(rnorm(32, sd=1), nrow=8, ncol=4)
pval_mat  <- matrix(runif(32, 0, 0.2), nrow=8, ncol=4)
rownames(logFC_mat) <- paste0("Gene", 1:8)
colnames(logFC_mat) <- paste0("Cond", 1:4)
rownames(pval_mat)  <- rownames(logFC_mat)
colnames(pval_mat)  <- colnames(logFC_mat)

ih_pre <- InformativeHeatmapFromMAT(
  logFC_matrix     = logFC_mat,
  pvalue_matrix    = pval_mat,
  pvalue_cutoff    = 0.05,
  trending_cutoff  = 0.1,
  pch_val          = 16,
  unit_val         = 4,
  significant_color= "black",
  trending_color   = "yellow",
  col              = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = TRUE,
  show_column_names= TRUE,
  name             = "logFC"
)

ht_pre <- getHeatmapObject(ih_pre)
ComplexHeatmap::draw(ht_pre, heatmap_legend_side = "right")
```

#### 4.2.3. MAE + iClusterPlus Mode

```r
library(MultiAssayExperiment)
library(iClusterPlus)

# Assume `tiny_mae` is a small MAE with two assays (see Manual Example 2C)
ih_mae <- InformativeHeatmapFromMAE(
  mae                = tiny_mae,
  significant_pvalue = 0.05,
  trending_pvalue    = 0.1,
  pch_val            = 16,
  unit_val           = 3,
  significant_color  = "darkred",
  trending_color     = "darkorange",
  K                  = 2,
  lambda             = 0.1,
  coef               = 2,
  heatmap_scale      = "logFC",
  col                = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
  cluster_rows       = TRUE,
  cluster_columns    = TRUE,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  name               = "logFC"
)

ht_mae <- getHeatmapObject(ih_mae)
ComplexHeatmap::draw(ht_mae, heatmap_legend_side = "right")
```

---

<a name="ih-api"></a>

### 4.3. API Details

#### 4.3.1. `InformativeHeatmap(data, meta=NULL, assayName=NULL, groupColumn=NULL, sampleType=NULL, timepoint=NULL, dataType=c("auto","continuous","count"), var_quantile=0.75, pvalue_cutoff=0.05, trending_cutoff=0.1, fc_cutoff=0.585, min_samples=3, max_features=NULL, runClustering=FALSE, K=3, lambda=0.2, coef=2, pch_val=16, unit_val=4, significant_color="black", trending_color="yellow", heatmap_scale=c("expression","logFC","foldChange"), max_users=NULL, …)`

* **Arguments**

  * **Input Modes**:

    1. `inherits(data, "MultiAssayExperiment")`: use MAE pipeline (see below).
    2. `is.matrix(data)` & `!is.null(meta)`: matrix + metadata DE pipeline (like `ClearScatterplot_table`).
    3. `is.matrix(data)` & user passed `pvalue_matrix` in `…`: ready-to-plot workflow (`InformativeHeatmapFromMAT` style).
  * `meta`: required if `data` is a matrix (rows = features, columns = samples).
  * `assayName`: required if `data` is an MAE.
  * `groupColumn`, `sampleType`, `timepoint`: metadata columns for DE/contrasts/faceting.
  * `dataType`: `"auto"`, `"continuous"`, or `"count"`—only if DE pipeline.
  * `var_quantile`: \[0–1] quantile to filter low-variance features before DE or plotting.
  * `pvalue_cutoff`, `trending_cutoff`: thresholds for “significant” and “trending” points.
  * `fc_cutoff`: |log₂FC| threshold to include features in heatmap.
  * `min_samples`: minimum samples per group per facet (DE).
  * `max_features`: cap on features by variance (DE pipelines) or by |logFC| (user can pre-filter).
  * `runClustering`, `K`, `lambda`, `coef`: only for MAE pipeline (iClusterPlus clustering + limma contrast).
  * `pch_val`, `unit_val`, `significant_color`, `trending_color`: point overlay aesthetics.
  * `heatmap_scale`: `"expression"`, `"logFC"`, or `"foldChange"`.
  * `…`: additional `ComplexHeatmap::Heatmap()` arguments (e.g., `col`, `cluster_rows`, `show_row_names`, etc.).

* **Behavior**

  1. **MAE pipeline**:

     * Extract each assay → numeric matrix (features × samples).
     * Transpose & run `iClusterPlus()` across assays → sample clusters.
     * For each assay, run `limma::lmFit()` + `eBayes()` to compare clusters → `logFC` & `p-values`.
     * Combine assays’ `logFC` & `p-values` → `combined_logFC` & `combined_pvalues` matrices (features × assays/clusters).
     * Proceed to step 3.

  2. **Matrix + Metadata pipeline**:

     * Exactly as `ClearScatterplot_table`, but instead of building a DE table, assemble one `logFC_matrix` and `pvalue_matrix` (features × facets).
     * `heatmap_data <- logFC_matrix`; `pvalue_matrix <- pvalue_matrix`.
     * Proceed to step 3.

  3. **Precomputed matrices**:

     * If user passed `pvalue_matrix` in `…`, skip DE and filtering.
     * Optionally apply `var_quantile` filter on `heatmap_data`.
     * Optionally cap to `max_features` by variance.

  4. **Vectorized `layer_fun`**:

     ```r
     layer_fun <- function(j, i, x, y, w, h, fill) {
       ind_mat     <- ComplexHeatmap::restore_matrix(j, i, x, y)
       feature_idx <- as.vector(row(ind_mat))
       cell_idx    <- as.vector(col(ind_mat))
       all_inds    <- as.vector(ind_mat)
       pvals_vec   <- pvalue_matrix[cbind(feature_idx, cell_idx)]
       color_vec   <- ifelse(
                        pvals_vec < pvalue_cutoff, significant_color,
                        ifelse(pvals_vec < trending_cutoff, trending_color, NA)
                      )
       keep_plots  <- !is.na(color_vec)
       grid::grid.points(
         x[all_inds][keep_plots],
         y[all_inds][keep_plots],
         pch  = pch_val,
         gp   = grid::gpar(col = color_vec[keep_plots]),
         size = grid::unit(unit_val, "mm")
       )
     }
     ```

  5. **Call** `ComplexHeatmap::Heatmap(heatmap_data, name=…, col=…, cluster_rows=…, cluster_columns=…, layer_fun=layer_fun, …)`.

  6. **Return** an **`InformativeHeatmap`** S4 object with slots:

     * `@heatmap`: the built `Heatmap` object
     * `@params`: the list of all passed arguments, for later updates

* **MAE Pipeline Specifics**

  ```r
  assays_list <- lapply(assays(mae), function(se) {
    m <- as.matrix(assay(se))
    mode(m) <- "numeric"
    t(m)  # samples × features for iClusterPlus
  })
  fit <- iClusterPlus::iClusterPlus(
           data_list = assays_list,
           type      = rep("gaussian", length(assays_list)),
           K         = K,
           lambda    = rep(lambda, length(assays_list)),
           maxiter   = 20
         )
  clusters <- factor(fit$clusters)
  design   <- model.matrix(~ 0 + clusters)
  limma_results <- lapply(names(assays(mae)), function(aname) {
    se_assay <- assays(mae)[[aname]]
    expr_assay <- SummarizedExperiment::assay(se_assay)
    fit_limma  <- limma::lmFit(expr_assay, design)
    fit_limma  <- limma::eBayes(fit_limma)
    tt         <- limma::topTable(fit_limma, coef = coef, number = Inf)
    list(logFC   = tt$logFC, p_values = tt$P.Value)
  })
  combined_logFC   <- do.call(rbind, lapply(limma_results, `[[`, "logFC"))
  combined_pvalues <- do.call(rbind, lapply(limma_results, `[[`, "p_values"))
  heatmap_data     <- t(combined_logFC)
  pvalue_matrix    <- t(combined_pvalues)
  ```

* **Exports**

  ```r
  @exportClass InformativeHeatmap
  @exportMethod InformativeHeatmap
  @export InformativeHeatmapFromMAT
  @export InformativeHeatmapFromMAE
  @exportMethod updateLayerFun
  @exportMethod getHeatmapObject
  ```

---

<a name="ih-examples"></a>

### 4.4. Examples

#### 4.4.1. Matrix + Metadata Mode (DE Heatmap)

```r
# Simulate or load data as in Quick Start
# (Example repeated from Quick Start 4.2.1)
ih_mat <- InformativeHeatmap(
  data             = expr,
  meta             = meta,
  groupColumn      = "Group",
  sampleType       = "SampleType",
  timepoint        = NULL,
  dataType         = "continuous",
  var_quantile     = 0.5,
  pvalue_cutoff    = 0.05,
  trending_cutoff  = 0.1,
  fc_cutoff        = 0.585,
  min_samples      = 3,
  max_features     = 30,
  runClustering    = FALSE,
  pch_val          = 16,
  unit_val         = 3,
  significant_color= "darkred",
  trending_color   = "orange",
  heatmap_scale    = "logFC",
  col              = circlize::colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = TRUE,
  show_column_names= TRUE,
  name             = "logFC"
)
ht_mat <- getHeatmapObject(ih_mat)
ComplexHeatmap::draw(ht_mat, heatmap_legend_side = "right")
```

#### 4.4.2. Precomputed Matrices Mode

```r
# (Example repeated from Quick Start 4.2.2)
ih_pre <- InformativeHeatmapFromMAT(
  logFC_matrix     = logFC_mat,
  pvalue_matrix    = pval_mat,
  pvalue_cutoff    = 0.05,
  trending_cutoff  = 0.1,
  pch_val          = 16,
  unit_val         = 4,
  significant_color= "black",
  trending_color   = "yellow",
  col              = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  show_row_names   = TRUE,
  show_column_names= TRUE,
  name             = "logFC"
)
ht_pre <- getHeatmapObject(ih_pre)
ComplexHeatmap::draw(ht_pre, heatmap_legend_side = "right")
```

#### 4.4.3. MAE + iClusterPlus Mode

```r
# (Example repeated from Quick Start 4.2.3)
ih_mae <- InformativeHeatmapFromMAE(
  mae                = tiny_mae,
  significant_pvalue = 0.05,
  trending_pvalue    = 0.1,
  pch_val            = 16,
  unit_val           = 3,
  significant_color  = "darkred",
  trending_color     = "darkorange",
  K                  = 2,
  lambda             = 0.1,
  coef               = 2,
  heatmap_scale      = "logFC",
  col                = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")),
  cluster_rows       = TRUE,
  cluster_columns    = TRUE,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  name               = "logFC"
)
ht_mae <- getHeatmapObject(ih_mae)
ComplexHeatmap::draw(ht_mae, heatmap_legend_side = "right")
```

---

<a name="ih-advanced"></a>

### 4.5. Advanced Customization & Tips

* **Custom “layer\_fun”**
  Replace the default significance overlay:

  ```r
  custom_layer <- function(j, i, x, y, w, h, fill) {
    ind_mat     <- ComplexHeatmap::restore_matrix(j, i, x, y)
    feature_idx <- as.vector(row(ind_mat))
    cell_idx    <- as.vector(col(ind_mat))
    all_inds    <- as.vector(ind_mat)
    # Access original pvalue_matrix (stored in ih_mat@params)
    pvals_vec   <- ih_mat@params$pvalue_matrix[cbind(feature_idx, cell_idx)]
    color_vec   <- ifelse(pvals_vec < 0.001, "magenta", NA)
    keep_plots  <- !is.na(color_vec)
    grid::grid.points(
      x[all_inds][keep_plots],
      y[all_inds][keep_plots],
      pch  = 17,
      gp   = grid::gpar(col = color_vec[keep_plots]),
      size = grid::unit(5, "mm")
    )
  }
  ih_mat <- updateLayerFun(ih_mat, custom_layer)
  ComplexHeatmap::draw(getHeatmapObject(ih_mat))
  ```

* **Per-Modality Color Palettes**
  When using an MAE with multiple assays, specify per-modality palettes as `col_AssayName` in `…`.

* **Cluster Control**
  Set `cluster_rows = FALSE` or `cluster_columns = FALSE` to disable clustering.

* **Scaling Options**

  * `heatmap_scale = "expression"`: fill tiles with raw expression.
  * `heatmap_scale = "logFC"`: fill with log₂ fold-change.
  * `heatmap_scale = "foldChange"`: same as “logFC” (alias).

* **Performance**

  * Use `var_quantile` to reduce feature set before DE.
  * Set `max_features` to cap DE to top N features by variance.
  * Adjust `parallel` & `BPPARAM` for large data.

---

<a name="multifeaturegrid"></a>

## 5. MultifeatureGrid (2D Tile-and-Point Heatmaps)

<a name="mg-what"></a>

### 5.1. What It Does

Creates a grid plot where:

* **X-axis**: a categorical variable (e.g., “tissue”)
* **Y-axis**: a categorical variable (e.g., “signaling”)
* **Tile fill** (`fill`): continuous variable (e.g., activation z-score)
* **Point overlay**:

  * **Color** = –log₁₀(p-value) (gradient from `lowColor` to `highColor`)
  * **Size** = number of items (e.g., genes) associated with that tile
* **Optional faceting** by a third categorical variable (e.g., “timePoint”)

Ideal for summarizing pathway or module activity across tissues/conditions with significance and gene count overlays.

---

<a name="mg-quickstart"></a>

### 5.2. Quick Start

```r
library(MultiModalGraphics)

# 1. Example data (built-in helper)
demo_df <- get_multifeature_grid_df()
str(demo_df)
# 'data.frame': 16 rows × 6 columns: 
#  $ tissue            : Factor 
#  $ signaling         : Factor 
#  $ Activation_z_score: num
#  $ p                 : num
#  $ number_of_genes   : int
#  $ timePoint         : Factor 

# 2. Construct MultifeatureGrid object
mg <- MultifeatureGrid(
  data           = demo_df,
  title          = "Demo Multifeature Grid",
  x_label        = "Tissue Type",
  y_label        = "Pathway",
  logpval_label  = "-Log10 P-Value",
  zscore_label   = "Activation Z-Score",
  numitems_label = "Gene Count",
  color_palette  = "RdYlBu",
  breaks         = seq(-2, 2, 0.5)
)

# 3. Plot the heatmap
plot_heatmap(
  mg,
  pValueColumn       = "p",
  lowColor           = "lightblue",
  highColor          = "darkred",
  borderColor        = "black",
  columnForNumber    = "number_of_genes",
  independantVariable = "timePoint"
)
```

---

<a name="mg-api"></a>

### 5.3. API Details

#### 5.3.1. `MultifeatureGrid(data, title="Heatmap", x_label="X Label", y_label="Y Label", logpval_label="-log10(p-value)", zscore_label="Activation z-score", numitems_label="Number of Genes", color_palette="RdYlBu", breaks=seq(-1,1,0.5))`

* **Arguments**

  * `data`: `data.frame` with at least:

    * `tissue` (factor/character) – X-axis
    * `signaling` (factor/character) – Y-axis
    * `Activation_z_score` (numeric) – tile fill
    * `p` (numeric p-value) – point color
    * `number_of_genes` (numeric) – point size
    * `timePoint` (optional, factor/character) – facet columns
  * `title`, `x_label`, `y_label`, `logpval_label`, `zscore_label`, `numitems_label`: plot text labels
  * `color_palette`: valid RColorBrewer palette (e.g., `"RdYlBu"`, `"PuOr"`)
  * `breaks`: numeric vector of breakpoints for tile fill

* **Validity Checks**

  * Ensures required columns exist; errors if missing.
  * Converts `tissue` & `signaling` (and `timePoint`) to factors if not already.
  * Verifies `color_palette` is in `RColorBrewer::brewer.pal.info`.

* **Returns**

  * A **`MultifeatureGrid`** S4 object with slots:

    * `@data`: input `data.frame` (unchanged)
    * Text and color configuration slots

#### 5.3.2. `plot_heatmap(object, pValueColumn="p", lowColor="yellow", highColor="red", borderColor="grey60", columnForNumber="number_of_genes", independantVariable="timePoint")`

* **Arguments**

  * `object`: a `MultifeatureGrid` instance
  * `pValueColumn`: column in `object@data` for p-values (numeric)
  * `lowColor`, `highColor`: colors for low/high –log₁₀(p) gradient
  * `borderColor`: color for tile borders
  * `columnForNumber`: column in `object@data` for point size mapping (numeric)
  * `independantVariable`: column in `object@data` for facet columns (factor/character)

* **Behavior**

  1. Checks required columns exist; stops if missing.

  2. Computes `neglog10p = -log10(data[[pValueColumn]])`.

  3. Builds a `ggplot2` object:

     * `geom_tile(fill = Activation_z_score, colour = borderColor)`
     * `scale_fill_gradientn(colours = palette, breaks = breaks)`
     * `geom_point(aes(colour = neglog10p, size = number_of_genes))`
     * `scale_color_gradient(low = lowColor, high = highColor)`
     * `scale_size(range = c(1, 10))`
     * `labs(...)` using labels from `object` slots
     * `facet_grid(cols = vars(timePoint), scales="free_x", space="free")`
     * `theme_bw()` + text styling (axis text rotated, etc.)

  4. Prints the plot to the active graphics device.

* **Returns**

  * Invisibly returns `NULL`; the plot is displayed by `print()`.

---

<a name="mg-examples"></a>

### 5.4. Examples

#### 5.4.1. Example A: Built-In Demo Data

```r
# (Same as Quick Start 4.2)
demo_df <- get_multifeature_grid_df()
mg <- MultifeatureGrid(
  data         = demo_df,
  title        = "Demo Multifeature Grid",
  x_label      = "Tissue",
  y_label      = "Pathway",
  logpval_label= "-Log10 P",
  zscore_label = "Z-Score",
  numitems_label="Gene Count",
  color_palette="RdYlBu",
  breaks       = seq(-2, 2, 0.5)
)
plot_heatmap(
  mg,
  pValueColumn       = "p",
  lowColor           = "lightgreen",
  highColor          = "darkgreen",
  borderColor        = "grey30",
  columnForNumber    = "number_of_genes",
  independantVariable= "timePoint"
)
```

#### 5.4.2. Example B: Custom Pathway Enrichment Table

```r
# Simulate a real-world pathway-enrichment result
set.seed(202)
df_pe <- data.frame(
  tissue            = factor(rep(c("Normal","Tumor"), each = 12)),
  signaling         = factor(rep(paste0("Pathway", 1:12), 2)),
  Activation_z_score= runif(24, -3, 3),
  p                 = runif(24, 1e-4, 0.05),
  number_of_genes   = sample(10:200, 24, TRUE),
  timePoint         = factor(sample(c("Baseline","FollowUp"), 24, replace = TRUE)),
  stringsAsFactors  = FALSE
)
df_pe$tissue    <- factor(df_pe$tissue, levels = c("Normal","Tumor"))
df_pe$signaling <- factor(df_pe$signaling, levels = paste0("Pathway", 1:12))

mg2 <- MultifeatureGrid(
  data           = df_pe,
  title          = "Taylor PCa Pathway Enrichment",
  x_label        = "Tissue",
  y_label        = "Pathway",
  logpval_label  = "-Log10 P-Value",
  zscore_label   = "Z-Score",
  numitems_label = "Gene Count",
  color_palette  = "PuOr",
  breaks         = seq(-3,3,1)
)

plot_heatmap(
  mg2,
  pValueColumn       = "p",
  lowColor           = "lightblue",
  highColor          = "darkblue",
  borderColor        = "grey50",
  columnForNumber    = "number_of_genes",
  independantVariable= "timePoint"
)
```

---

<a name="mg-advanced"></a>

### 5.5. Advanced Customization & Tips

* **Removing Facets**
  If your data has no `timePoint`, set `independantVariable = NULL` or to a constant column.

* **Adjusting Tile Breaks**
  The `breaks` vector in the constructor defines where to place color ticks. You can supply any numeric sequence that suits your data range.

* **Customizing Point Size Scale**
  After `plot_heatmap()`, you can tweak point sizing:

  ```r
  p <- plot_heatmap(mg)
  p + scale_size(range = c(0.5, 15))
  ```

* **Changing Font & Theme**
  Since the plot is a `ggplot2` object, append any theme modifications:

  ```r
  p <- plot_heatmap(mg)
  p + theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 12)
  )
  ```

* **Large Datasets**
  For many tiles (e.g., thousands of combinations), consider:

  * Reducing factor levels (e.g., top 50 pathways only)
  * Splitting into multiple plots per facet subset
  * Adjusting `coord_fixed()` to manage aspect ratio

---

<a name="multimodal"></a>

## 6. Unified Multimodal API: `MultiModalPlot`

Rather than calling each class separately, use **`MultiModalPlot()`** to automatically dispatch inputs to their correct constructors and combine plots side-by-side. Supported input types:

* **`MultiAssayExperiment` objects** → `ClearScatterplot_MAE` or `InformativeHeatmap`
* **Lists of `expr` + `meta`** → `ClearScatterplot_table` or `InformativeHeatmap`
* **Precomputed DE tables** → `ClearScatterplot` (volcano only)
* **Fold-change + p-value matrices** (precomputed) → `InformativeHeatmapFromMAT` (heatmap only)

**Common Arguments**:

```r
MultiModalPlot(
  inputs,                # named list of inputs
  assayNames     = NULL, # per-modality assay names (for MAE inputs)
  groupColumns   = NULL, # per-modality grouping columns
  sampleTypes    = NULL, # per-modality X-facet columns
  timepoints     = NULL, # per-modality Y-facet columns
  dataType       = "auto",
  var_quantile   = 0.75,
  pvalue_cutoff  = 0.05,
  trending_cutoff= 0.1,
  fc_cutoff      = 0.585,
  max_features   = NULL,
  parallel       = TRUE,
  BPPARAM        = BiocParallel::bpparam(),
  panel_type     = c("volcano","heatmap"),
  …                  # additional class-specific args
)
```

* **`panel_type = "volcano"`**

  * Single modality → returns a `ggplot` (volcano).
  * Multiple modalities → returns a `patchwork`‐combined plot (volcano panels side-by-side).

* **`panel_type = "heatmap"`**

  * Single modality → returns a `Heatmap` object.
  * Multiple modalities → returns a `HeatmapList` (heatmap panels combined with `+`).

> **Example**
>
> ```r
> # Two modalities: MAE for RNA, and precomputed DE table
> volcano_panel <- MultiModalPlot(
>   inputs       = list(RNA = miniACC, DEtable = df_de),
>   assayNames   = c(RNA = "RNASeq2GeneNorm"),
>   groupColumns = c(RNA = "C1A.C1B"),
>   sampleTypes  = c(RNA = "pathologic_stage"),
>   timepoints   = c(RNA = "MethyLevel"),
>   panel_type   = "volcano",
>   dataType     = "auto",
>   parallel     = FALSE,
>   color1       = "navy",
>   color3       = "firebrick",
>   title        = c(RNA = "miniACC Volcano", DEtable = "Precomputed Volcano")
> )
> print(volcano_panel)
> ```

---

<a name="session"></a>

## 7. Session Info & Citation

For reproducibility:

```r
sessionInfo()
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Attached base packages: stats, graphics, grDevices, utils, datasets, methods, base
# Other attached packages: 
#   MultiModalGraphics_0.99.5 ggplot2_3.4.0 dplyr_1.0.10 matrixStats_0.63.1
#   BiocParallel_1.30.0 limma_3.52.0 SummarizedExperiment_1.26.0 MultiAssayExperiment_2.8.0
#   ComplexHeatmap_2.6.0 circlize_0.4.0 RColorBrewer_1.1-3
# Loaded via a namespace (and not attached):
#   BiocGenerics_0.36.0 IRanges_2.32.0 S4Vectors_0.36.0 Biobase_2.56.0
#   grid_4.2.0 gridGraphics_0.5-1 magrittr_2.0.3 etc.
```

Please cite:

> Mohammed FA, Fall M, Hammamieh H, Muhie S (2025). **MultiModalGraphics: Multi-Modal Visualization & Plotting Utilities**. R package version 0.99.5.

---

**Thank you for using MultiModalGraphics!**
Feedback, issues, and pull requests are welcome on [GitHub](https://github.com/famanalytics0/MultiModalGraphics).

