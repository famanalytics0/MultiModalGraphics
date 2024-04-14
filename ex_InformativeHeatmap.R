library(ComplexHeatmap)
library(MultiModalGraphics)

informative_heatmap = read.csv("InformativeHeatmap.csv", row.names = 1)
informative_heatmap_matrix = as.matrix(informative_heatmap)
group_val = informative_heatmap_matrix[,1:3]
p_val = informative_heatmap_matrix[,4:6]

htmp_plus <- InformativeHeatmap(group_val,
                                 unit_val=7,
                                 pch_val=16,
                                 signicant_color="black",
                                 trending_color="turquoise",
                                 significant_pvalue=0.05,
                                 trending_pvalue=0.1,
                                 significance_level = p_val,
                                 row_title = "Genes",
                                 column_title = "Value and Significance",
                                 cluster_rows = T,
                                 show_row_names = T,
                                 row_names_side = c("left"),
                                 column_names_rot = 45,
                                 row_dend_reorder = TRUE,
                                 rect_gp = gpar(col = "white", lwd = 2))

draw(getHeatmapObject(htmp_plus))
