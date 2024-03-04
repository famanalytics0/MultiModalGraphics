setwd(".")
library(MilkR)
library(magrittr)

# Load or source the ClearScatterplot.R file
#source("R/ClearScatterplot.R")

# Create an instance of the ClearScatterplot class
#data_filepath <- "scattered_data.csv"
#filepath <- "scatter_plot_data.csv"
filepath <- "ClearScatterplot_good.csv"
#filepath <- "../MilkR/mPTSD blood heart spleen brain regions DEGs for scattor plot_good_v2.csv"
plotdata <- read.csv(filepath)
View(plotdata)
# data <- get_clear_scatterplot_df()
#data
scatterplotObject <- ClearScatterplot(data = plotdata,
                        pValueColumn = "p",
                        expressionColumnName = "log2fc",
                        timePointColumn = "timePoint",
                        timePointLevels = c("T10R1", "T5R1"),
                        colorHigh = "cornflowerblue", colorNeutral = "grey", colorLow = "indianred"
                        )

# data, pValueColumn = "p", expressionColumnName = "log2fc",
# highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue = 1.301,
# timePointColumn = "timePoint", timePointLevels = NULL,
# colorHigh = "cornflowerblue", colorNeutral = "grey", colorLow = "indianred"
# scattered_plot <- new("ClearScatterplot", data = plotdata,  pValueColumn = "p",
#                       expressionColumnName = "log2fc",
#                       highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301)

# (data, pValueColumn = "p", expressionColumnName = "log2fc",
#   highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue = 1.301,
#   timePointColumn = "timePoint", timePointLevels = NULL,
#   colorHigh = "cornflowerblue", colorNeutral = "grey", colorLow = "indianred")

scatterplotObject <- createPlot(scatterplotObject,
                        highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue = 1.301,
                        color1 = "cornflowerblue", color2 = "grey", color3 = "indianred",
                        expressionDirection = "regulation", timeVariable="reg_time_org",
                        xAxis = "organ", yAxis = "timePoint")

# color1 = "cornflowerblue", color2 = "grey", color3="indianred",
# highLog2fc = 0.585, lowLog2fc = -0.585, expressionDirection = "regulation",
# timeVariable="reg_time_org", xAxis = "organ", yAxis = "timePoint"
# Call methods on the object
scattered_plot <- createPlot(scatterplotObject, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                             highLog2fc = 0.585, lowLog2fc = -0.585,negLog10pValue = 1.301,
                             expressionDirection = "regulation",
                             timeVariable="reg_time_org")  # Create the plot

# Print the plot
scattered_plot  # This will call the 'show' method and display the plot

