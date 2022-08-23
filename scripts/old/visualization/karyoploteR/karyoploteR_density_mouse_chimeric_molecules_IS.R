#---- Installation ----

# Install BiocManager and karyoploteR
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("karyoploteR")

# https://stackoverflow.com/questions/41839214/installation-path-not-writable-r-unable-to-update-packages
#BiocManager::install("karyoploteR")

#install.packages("RColorBrewer")


#---- Library ----

library(karyoploteR)
library(readr)
library(RColorBrewer)



#---- Set data & parameters <--- TO CHANGE! ----

# Set working directory
workdir <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/mouse/IS_chimeric_molecules/"
mouse_genome <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/mouse/mouse_genome.txt"
sample <- list.files(path="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/mouse/IS_chimeric_molecules/", pattern="*_mouse_chimeric_IS.bed")
color <- brewer.pal(n = 9, name = "Set1")
# Display all the colors
display.brewer.all()


#---- Create a function to plot an ideogram ----

# Create a function to plot the ideogram using a custom genome
karyoploteR_plot <- function(workdir, genome) {
  
  # For AAV genome, create a custom genome ideogram
  custom.genome <- toGRanges(genome)
  
  # Plot the ideogram using the custom genome
  kp <- plotKaryotype(plot.type=1, genome = custom.genome)
  
  return(kp)
  
}


#---- Create a function to add IS on mouse genome ----

# Create a function to calcul the density of IS on mouse genome per sample
karyoploteR_IS_plot <- function(workdir, color, path_to_data, pos0, pos1) {
  
  # Import IS data
  data <- read.delim(paste0(workdir, path_to_data), header=TRUE)
  
  # Create a GRanges object to store IS data
  IS.data <- toGRanges(data.frame(chr=data[, "chr"], start=data[, "start"], end=data[, "end"]))
  
  # Add IS data 
  kp <- kpPlotDensity(kp, IS.data, window.size = 10e6, col=color)
  
  # Add axis (with ymax value)
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density)
  
  return(kp)
  
}



#---- [Final] Plot the ideogram with additional data ----

# Source: https://bernatgel.github.io/karyoploter_tutorial//Examples/GeneDensity/GeneDensity.html

# Plotting parameters to reduce blank space
pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20

# Plot the ideogram with the mouse reference genome
kp <- plotKaryotype(genome="mm10", plot.type=4, ideogram.plotter=NULL, labels.plotter=NULL)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=45)
# Execute the function to plot an mouse custom genome ideogram
#kp <- karyoploteR_plot(workdir, mouse_genome)

# Plot additional data: IS data for each sample
## Add regions for each sample
karyoploteR_IS_plot(workdir, color[1], sample[1], 0, 0.2) # 149
karyoploteR_IS_plot(workdir, color[4], sample[4], 0.2, 0.4) # 152
karyoploteR_IS_plot(workdir, color[2], sample[2], 0.4, 0.6) # 150
karyoploteR_IS_plot(workdir, color[3], sample[3], 0.6, 0.8) # 151
karyoploteR_IS_plot(workdir, color[5], sample[5], 0.8, 1.0) # 153
karyoploteR_IS_plot(workdir, color[9], sample[6], 1.0, 1.2) # 154

