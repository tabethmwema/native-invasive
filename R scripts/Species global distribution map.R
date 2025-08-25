# Clear the environment history
rm(list = ls()) 

# Set Working Directory
setwd("C:/Users/tzm0087/OneDrive/PartA")

# Install packages
install.packages("maps")

# Load packages
library(maps)
library(ggplot2)

# Load the data
cast_data <- read.csv("CAST_lat_long.csv", header = TRUE, sep = ",")
dom_data <- read.csv("DOM_lat_long.csv", header = TRUE, sep = ",")

# Set the dimensions for the plot
png("Mwema_mapA.png", res = 300, width = 190, height = 110, units = "mm")

# Plot a world map
map("world", fill = TRUE, col = "lightgrey", bg = "lightcyan1", mar = c(1, 3, 1, 3))

# Plot the coordinates
points(cast_data$long, cast_data$lat, pch=16, col="yellow2", cex = 0.5) # pch = 16 for a solid circle
points(dom_data$long, dom_data$lat, pch=17, col="tomato4", cex = 0.5) # pch = 17 for a solid triangle

# Add lakes
map('lakes', add=TRUE, fill=TRUE, col="cyan1", boundary="cyan1")

# Fill Antarctica with white color
map('world', regions='Antarctica', col='white', add=TRUE, fill = TRUE)

# Add axes around the entire map
axis(1, cex.axis = 0.7)  # x-axis at the bottom
axis(2, cex.axis = 0.7)  # y-axis on the left
axis(3, cex.axis = 0.7)  # x-axis at the top
axis(4, cex.axis = 0.7)  # y-axis on the right

# Connect the axes to form a box
box()

# Add a title
title(main = "Species Distribution Map", col.main = "black", cex.main = 1.2, line = +2.5)
title(xlab = "Latitude", ylab = "Longitude", line = 2)

# Add a legend
legend("bottomleft", legend = c("M. m. castaneus", "M. m. domesticus"), pch = c(16, 17), col = c("yellow2", "tomato4"), cex = 0.9)

# Save and close the plot
dev.off()










