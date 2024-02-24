#' Provides a sample dataframe for informative heatmap visualization
#'
#' This function returns a predefined dataframe that represents data suitable for
#' creating an informative heatmap. The dataframe includes values for three groups
#' and their associated p-values, with row names indicating gene identifiers.
#'
#' @return A dataframe with sample data for an informative heatmap.
#' @examples
#' informative_heatmap_df <- get_informative_heatmap_df()
#' @export
get_informative_heatmap_df <- function() {
  informative_heatmap_df <- structure(list(group1_value = c(0.42, 0.13, -0.32, -0.42),
                                 group2_value = c(-0.36, 0.41, 0.21, 0.32),
                                 group3_value = c(-0.17, -0.15, 0.12, 0.11),
                                 pvalue1 = c(0.01, 0.12, 0.12, 0.02),
                                 pvalue2 = c(0.07, 0.03, 0.21, 0.08),
                                 pvalue3 = c(0.07, 0.3, 0.06, 0.21)),
                            class = "data.frame", row.names = c("gene1", "gene2", "gene3", "gene4"))

  return(informative_heatmap_df)
}

#' Provides a sample dataframe for clear scatterplot visualization
#'
#' This function returns a predefined dataframe that represents data suitable for
#' creating a clear scatterplot. The dataframe includes log2 fold change values, p-values,
#' and regulation direction for a set of genes across different tissues and time points.
#'
#' @return A dataframe with sample data for a clear scatterplot.
#' @examples
#' clear_scatterplot_df <- get_clear_scatterplot_df()
#' @export
get_clear_scatterplot_df <- function() {
  clear_scatterplot <- structure(list(
    log2fc = c(-3.01, 1.19, -1.39,
               -0.99, 1.99, -1.54, -1.21, -1.55, 1.85, -1.28, -1.3, 1.3,
               -1.08, -0.85, -2.62, 1.69, 1.96, -1.06, -1.51, 0.79, -1.45,
               -1.64, -1.26, 1.76, 1.48, -1.12, 1.14, 1.63, 1.21, -2.07,
               -1.49, 2.56, -0.83, 1.94, 1.8, 0.85, 1.91, 1.75, 1.36, -1.33,
               -0.8, -2.03, -1.67, 0.79, -1.34, 0.94, -0.75, 1.01, -1.09,
               0.94, -0.64, 0.98, -1.14, -0.82, -0.76, -0.69, 0.78, -0.62,
               0.32, 0.9, 0.36, -0.58, 0.38, 0.44, -0.49, 0.39, -0.31, 0.23,
               0.73, 0.52, 0.39, 0.24, -0.36, -0.37, -0.24, -0.26, -0.32,
               -0.36, 0.37, -0.58, 0.31, -0.26, 0.32, -0.24, 0.42, -0.15,
               -0.25, 0.09, -0.07, -0.17, 0.16, -0.05, 0.1, -0.04, -0.1,
               -0.05, 0.03, -0.04, 0.01, -0.01),

    p = c(1e-04, 0.0064, 0.0076,
          0.0101, 0.0101, 0.0122, 0.0131, 0.0159, 0.0182, 0.0207, 0.0226,
          0.0248, 0.0251, 0.0267, 0.0278, 0.0281, 0.03, 0.0309, 0.0314,
          0.0314, 0.0319, 0.0333, 0.0338, 0.0351, 0.0357, 0.0358, 0.0372,
          0.0385, 0.0392, 0.0395, 0.041, 0.0411, 0.0416, 0.0421, 0.0425,
          0.0436, 0.0445, 0.0459, 0.0473, 0.0496, 0.0555, 0.0702, 0.073,
          0.0971, 0.1327, 0.1566, 0.1578, 0.1673, 0.168, 0.1725, 0.1799,
          0.1812, 0.2052, 0.2276, 0.2409, 0.249, 0.2915, 0.321, 0.3332,
          0.347, 0.3613, 0.3882, 0.3913, 0.429, 0.4293, 0.445, 0.4719,
          0.4776, 0.4803, 0.491, 0.4986, 0.4987, 0.5184, 0.5317, 0.5477,
          0.555, 0.5675, 0.5763, 0.6052, 0.6241, 0.626, 0.6849, 0.6883,
          0.6892, 0.6966, 0.7444, 0.7644, 0.7883, 0.7902, 0.7919, 0.8191,
          0.833, 0.8824, 0.9051, 0.911, 0.9186, 0.9222, 0.9314, 0.9848,
          0.989),

    regulation = c("decrease", "increase", "decrease",
                   "decrease", "increase", "decrease", "decrease", "decrease",
                   "increase", "decrease", "decrease", "increase", "decrease",
                   "decrease", "decrease", "increase", "increase", "decrease",
                   "decrease", "increase", "decrease", "decrease", "decrease",
                   "increase", "increase", "decrease", "increase", "increase",
                   "increase", "decrease", "decrease", "increase", "decrease",
                   "increase", "increase", "increase", "increase", "increase",
                   "increase", "decrease", "decrease", "decrease", "decrease",
                   "increase", "decrease", "increase", "decrease", "increase",
                   "decrease", "increase", "decrease", "increase", "decrease",
                   "decrease", "decrease", "decrease", "increase", "decrease",
                   "increase", "increase", "increase", "decrease", "increase",
                   "increase", "decrease", "increase", "decrease", "increase",
                   "increase", "increase", "increase", "increase", "decrease",
                   "decrease", "decrease", "decrease", "decrease", "decrease",
                   "increase", "decrease", "increase", "decrease", "increase",
                   "decrease", "increase", "decrease", "decrease", "increase",
                   "decrease", "decrease", "increase", "decrease", "increase",
                   "decrease", "decrease", "decrease", "increase", "decrease",
                   "increase", "decrease"),

    organ = c("tissue2", "tissue2",
              "tissue2", "tissue2", "tissue2", "tissue1", "tissue2", "tissue2",
              "tissue2", "tissue2", "tissue1", "tissue2", "tissue2", "tissue2",
              "tissue1", "tissue1", "tissue2", "tissue2", "tissue2", "tissue2",
              "tissue2", "tissue1", "tissue2", "tissue1", "tissue1", "tissue2",
              "tissue2", "tissue2", "tissue2", "tissue1", "tissue2", "tissue1",
              "tissue2", "tissue1", "tissue1", "tissue2", "tissue1", "tissue1",
              "tissue1", "tissue1", "tissue2", "tissue1", "tissue1", "tissue2",
              "tissue1", "tissue1", "tissue1", "tissue1", "tissue1", "tissue1",
              "tissue1", "tissue1", "tissue1", "tissue1", "tissue1", "tissue2",
              "tissue2", "tissue1", "tissue2", "tissue1", "tissue2", "tissue1",
              "tissue2", "tissue1", "tissue1", "tissue1", "tissue1", "tissue2",
              "tissue1", "tissue1", "tissue2", "tissue2", "tissue1", "tissue1",
              "tissue1", "tissue2", "tissue2", "tissue1", "tissue1", "tissue1",
              "tissue1", "tissue2", "tissue2", "tissue1", "tissue1", "tissue1",
              "tissue1", "tissue2", "tissue2", "tissue1", "tissue2", "tissue2",
              "tissue2", "tissue2", "tissue1", "tissue1", "tissue2", "tissue1",
              "tissue2", "tissue1"),

    timePoint = structure(c(1L, 1L, 1L,
                            1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L,
                            1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L,
                            2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L,
                            2L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L,
                            1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L,
                            1L, 1L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L,
                            1L, 2L, 2L, 1L, 1L, 1L, 2L),

                          levels = c("time1", "time2"),

                          class = "factor"),

    reg_time_org = c("decrease.tissue2.time1", "increase.tissue2.time1",
                     "decrease.tissue2.time1", "decrease.tissue2.time1", "increase.tissue2.time1",
                     "decrease.tissue1.time2", "decrease.tissue2.time1", "decrease.tissue2.time1",
                     "increase.tissue2.time1", "decrease.tissue2.time1", "decrease.tissue1.time1",
                     "increase.tissue2.time1", "decrease.tissue2.time1", "decrease.tissue2.time1",
                     "decrease.tissue1.time2", "increase.tissue1.time1", "increase.tissue2.time1",
                     "decrease.tissue2.time1", "decrease.tissue2.time1", "increase.tissue2.time2",
                     "decrease.tissue2.time1", "decrease.tissue1.time1", "decrease.tissue2.time2",
                     "increase.tissue1.time2", "increase.tissue1.time1", "decrease.tissue2.time1",
                     "increase.tissue2.time1", "increase.tissue2.time1", "increase.tissue2.time1",
                     "decrease.tissue1.time2", "decrease.tissue2.time1", "increase.tissue1.time2",
                     "decrease.tissue2.time1", "increase.tissue1.time2", "increase.tissue1.time2",
                     "increase.tissue2.time1", "increase.tissue1.time1", "increase.tissue1.time2",
                     "increase.tissue1.time2", "decrease.tissue1.time2", "decrease.tissue2.time2",
                     "decrease.tissue1.time2", "decrease.tissue1.time2", "increase.tissue2.time1",
                     "decrease.tissue1.time2", "increase.tissue1.time2", "decrease.tissue1.time2",
                     "increase.tissue1.time2", "decrease.tissue1.time2", "increase.tissue1.time1",
                     "decrease.tissue1.time1", "increase.tissue1.time1", "decrease.tissue1.time2",
                     "decrease.tissue1.time1", "decrease.tissue1.time2", "decrease.tissue2.time1",
                     "increase.tissue2.time1", "decrease.tissue1.time2", "increase.tissue2.time2",
                     "increase.tissue1.time2", "increase.tissue2.time1", "decrease.tissue1.time2",
                     "increase.tissue2.time1", "increase.tissue1.time1", "decrease.tissue1.time1",
                     "increase.tissue1.time2", "decrease.tissue1.time2", "increase.tissue2.time2",
                     "increase.tissue1.time2", "increase.tissue1.time1", "increase.tissue2.time1",
                     "increase.tissue2.time2", "decrease.tissue1.time2", "decrease.tissue1.time2",
                     "decrease.tissue1.time1", "decrease.tissue2.time2", "decrease.tissue2.time1",
                     "decrease.tissue1.time1", "increase.tissue1.time1", "decrease.tissue1.time1",
                     "increase.tissue1.time2", "decrease.tissue2.time1", "increase.tissue2.time2",
                     "decrease.tissue1.time1", "increase.tissue1.time2", "decrease.tissue1.time2",
                     "decrease.tissue1.time1", "increase.tissue2.time1", "decrease.tissue2.time2",
                     "decrease.tissue1.time2", "increase.tissue2.time2", "decrease.tissue2.time2",
                     "increase.tissue2.time1", "decrease.tissue2.time1", "decrease.tissue1.time2",
                     "decrease.tissue1.time2", "increase.tissue2.time1", "decrease.tissue1.time1",
                     "increase.tissue2.time1", "decrease.tissue1.time2"),

    neglog10p = c(4, 2.19382002601611, 2.11918640771921,
                  1.99567862621736, 1.99567862621736, 1.91364016932525,
                  1.88272870434424, 1.79860287567955, 1.73992861201493,
                  1.68402965454308, 1.6458915608526, 1.60554831917378,
                  1.60032627851896, 1.57348873863542, 1.55595520408192,
                  1.55129368009492, 1.52287874528034, 1.51004152057517,
                  1.50307035192679, 1.50307035192679, 1.49620931694282,
                  1.47755576649368, 1.47108329972235, 1.45469288353418,
                  1.44733178388781, 1.44611697335613, 1.4294570601181,
                  1.4145392704915, 1.40671393297954, 1.40340290437354,
                  1.38721614328026, 1.38615817812393, 1.38090666937326,
                  1.37571790416433, 1.37161106994969, 1.36051351073141,
                  1.35163998901907, 1.33818731446274, 1.32513885926219,
                  1.3045183235098, 1.25570701687732, 1.15366288787019,
                  1.13667713987954, 1.012780770092, 0.877129077135564,
                  0.805208242278075, 0.801893001126599, 0.776504059037605,
                  0.774690718274137, 0.763210900590707, 0.744968836654449,
                  0.741841806659206, 0.687822643560221, 0.642827742276966,
                  0.618163200001657, 0.603800652904264, 0.535361440904967,
                  0.493494967595128, 0.47729500726525, 0.459670525209126,
                  0.442132038431978, 0.410944468947656, 0.40749015209932,
                  0.367542707815276, 0.367239111520561, 0.351639989019068,
                  0.326150022657051, 0.320935681878687, 0.318487413361038,
                  0.308918507877032, 0.302247725832245, 0.302160631781637,
                  0.285335007137463, 0.274333339685821, 0.261457259071215,
                  0.255707016877324, 0.24603413413484, 0.239351380418644,
                  0.218101080648851, 0.204745817419117, 0.20342566678957,
                  0.164372833790102, 0.162222230446267, 0.161654731224009,
                  0.157016529877783, 0.128193635541271, 0.116679321617025,
                  0.103308473437116, 0.102262974654573, 0.10132965703447,
                  0.0866630740673768, 0.0793549985932124, 0.0543345005678658,
                  0.0433034351053491, 0.0404816230270017, 0.0368735589180954,
                  0.0351748821166112, 0.0308637664032875, 0.0066519600767401,
                  0.0048037084028206),

    color_flag = c(-1, 1, -1, -1, 1,
                   -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1, 1,
                   -1, -1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, 1,
                   1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),

    row.names = c(NA, -100L),

    class = "data.frame")


  return(clear_scatterplot)
}

#' Provides a sample dataframe for multifeature grid visualization
#'
#' This function returns a predefined dataframe that represents data suitable for
#' creating a multifeature grid visualization. The dataframe includes signaling pathways,
#' activation z-scores, and associated p-values for different tissues and time points.
#'
#' @return A dataframe with sample data for a multifeature grid visualization.
#' @examples
#' multifeature_grid_df <- get_multifeature_grid_df()
#' @export
get_multifeature_grid_df <- function() {
  multifeature_grid <- structure(list(signaling = c("pathway1", "pathway2", "pathway3",
                                                    "pathway4", "pathway1", "pathway2", "pathway3", "pathway4", "pathway1",
                                                    "pathway2", "pathway3", "pathway4", "pathway1", "pathway2", "pathway3",
                                                    "pathway4", "pathway1", "pathway2", "pathway3", "pathway4", "pathway1",
                                                    "pathway2", "pathway3", "pathway4", "pathway1", "pathway2", "pathway3",
                                                    "pathway4", "pathway1", "pathway2", "pathway3", "pathway4", "pathway1",
                                                    "pathway1", "pathway2", "pathway2", "pathway3", "pathway3", "pathway4",
                                                    "pathway4"),
                                      Activation_z_score = c(-1.37, -0.34, -2.227, -1.642,
                                              -2.15, 1.352, 0.215, 0.013, 1.37, -2.193, -0.738, 1.433, 1.37,
                                              1.268, -0.563, 0.562, 0.033, -1.974, -1.041, 0.247, 1.37, -0.359,
                                              -0.6, -0.285, 1.37, -0.905, -1.516, -1.606, 0.813, 0.897, 0.766,
                                              0.937, NA, NA, NA, NA, NA, NA, NA, NA),
                                      timePoint = c("  acute  stress",
                                              "  acute  stress", "  acute  stress", "  acute  stress", "  acute  stress",
                                              "  acute  stress", "  acute  stress", "  acute  stress", "  acute recovery",
                                              "  acute recovery", "  acute recovery", "  acute recovery", "  acute recovery",
                                              "  acute recovery", "  acute recovery", "  acute recovery", " chronic  stress",
                                              " chronic  stress", " chronic  stress", " chronic  stress", " chronic  stress",
                                              " chronic  stress", " chronic  stress", " chronic  stress", " chronic recovery",
                                              " chronic recovery", " chronic recovery", " chronic recovery",
                                              " chronic recovery", " chronic recovery", " chronic recovery",
                                              " chronic recovery", "DEG", "DEG", "DEG", "DEG", "DEG", "DEG",
                                              "DEG", "DEG"),
                                      tissue = c("tissue1", "tissue1", "tissue1", "tissue1",
                                              "tissue2", "tissue2", "tissue2", "tissue2", "tissue1", "tissue1",
                                              "tissue1", "tissue1", "tissue2", "tissue2", "tissue2", "tissue2",
                                              "tissue1", "tissue1", "tissue1", "tissue1", "tissue2", "tissue2",
                                              "tissue2", "tissue2", "tissue1", "tissue1", "tissue1", "tissue1",
                                              "tissue2", "tissue2", "tissue2", "tissue2", "genes", "genes",
                                              "genes", "genes", "genes", "genes", "genes", "genes"),
                                      p = c(NA,
                                               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4.49e-42,
                                               1.15e-42, 1.64e-20, 1.42e-21, 6.03e-21, 1.73e-22, 7.95e-30, 9.14e-30
                                              ),
                                      number_of_genes = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                               NA, NA, NA, NA, NA, NA, 40L, 40L, 26L, 27L, 23L, 24L, 40L, 40L
                                          )),
                                 class = "data.frame", row.names = c(NA, -40L))

  return(multifeature_grid)
}
