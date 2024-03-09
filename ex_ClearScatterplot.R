library(MultiModalGraphics)
library(magrittr)

plotdata <- get_clear_scatterplot_df()

scatterplotObject <- ClearScatterplot(data = plotdata,
                        logFoldChange = "log2fc",
                        timePointColumn = "timePoint",
                        timePointLevels = c("T10R1", "T5R1")
                        )

# Create the plot
scattered_plot <- createPlot(scatterplotObject, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                             highLog2fc = 0.585, lowLog2fc = -0.585,negLog10pValue = 1.301,
                             expressionDirection = "regulation", negativeLogPValue="negLog10p",
                              timeVariable="reg_time_org", xAxis = "organ", yAxis = "timePoint")

# Print the plot
scattered_plot  # This will call the 'show' method and display the plot

