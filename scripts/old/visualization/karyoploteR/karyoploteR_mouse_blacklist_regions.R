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
workdir <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/karyoploteR/input/mouse/"
mouse_genome <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/karyoploteR/input/mouse/mouse_genome.txt"
blacklist_regions <- list.files(path="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/karyoploteR/input/mouse/", pattern="mouse_blacklist_regions_*")
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



#---- Create a function to add blacklist regions ----

# Create a function to add coverage data on the ideogram plot
karyoploteR_blacklist_regions_plot <- function(workdir, path_to_similar_regions, path_to_repeat_regions, ymax) {
  
  # Import blacklist regions data
  similar_regions <- read.delim(paste0(workdir, path_to_similar_regions))
  repeat_regions <- read.delim(paste0(workdir, path_to_repeat_regions))
  
  # Create a GRanges object to store blacklist regions informations
  similar.data <- toGRanges(data.frame(chr=similar_regions$chr, start=similar_regions$start, end=similar_regions$end))
  repeat.data <- toGRanges(data.frame(chr=repeat_regions$chr, start=repeat_regions$start, end=repeat_regions$end))
  
  # Add blacklist regions
  kp.repeat <- kpRect(kp, data = repeat.data, y0=0, y1=ymax, col="#FFDDDD", border=NA) # red
  kp.similar <- kpRect(kp, data = similar.data, y0=0, y1=ymax, col="#DDFFDD", border=NA) # green
  
  return(kp.similar)
  
}



#---- [Final] Plot the ideogram with additional data ----

# Plot the ideogram with the mouse reference genome
kp <- plotKaryotype(genome="mm10", layer.margin = 0.1)
# Execute the function to plot an mouse custom genome ideogram
#kp <- karyoploteR_plot(workdir, mouse_genome)

# Add blacklist regions
## Add ymax
## TO DO! Change workdir directory!
kp <- karyoploteR_blacklist_regions_plot(workdir, blacklist_regions[2], blacklist_regions[1], 1)




