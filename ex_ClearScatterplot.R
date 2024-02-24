setwd(".")
library(MilkR)

# Load or source the ClearScatterplot.R file
#source("R/ClearScatterplot.R")

# Create an instance of the ClearScatterplot class
#data_filepath <- "scattered_data.csv"
#data_filepath <- "scatter_plot_data.csv"
filepath <- "ClearScatterplot.csv"
#data_filepath <- "mPTSD blood heart spleen brain regions DEGs for scattor plot_good_v2.csv"
#data <- read.csv(filepath)
data <- get_clear_scatterplot_df()

scattered_plot <- new("ClearScatterplot", data = data,
                      significanceColumn = "p", expressionColumnName = "log2fc",
                      highLog2fc = 0.585, lowLog2fc = -0.585, negLog10pValue =  1.301)

# Call methods on the object
scattered_plot <- createPlot(scattered_plot, color1 = "cornflowerblue", color2 = "grey", color3="indianred",
                             highLog2fc = 0.585, lowLog2fc = -0.585,
                             expressionDirection = "regulation",
                             timeVariable="reg_time_org")  # Create the plot

# Print the plot
scattered_plot  # This will call the 'show' method and display the plot

