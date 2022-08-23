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
blacklist_regions <- list.files(path="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/karyoploteR/input/mouse/", pattern="mouse_blacklist_regions_*")
color <- brewer.pal(n = 9, name = "Set1")
# Display all the colors
display.brewer.all()


#---- Create a function to add blacklist regions ----

# Create a function to add coverage data on the ideogram plot
karyoploteR_blacklist_regions_plot <- function(workdir, path_to_similar_regions, path_to_repeat_regions) {
  
  # Import blacklist regions data
  repeat_regions <- read.delim(paste0(workdir, path_to_repeat_regions))
  similar_regions <- read.delim(paste0(workdir, path_to_similar_regions))
  
  # Create a GRanges object to store blacklist regions informations
  repeat.data <- toGRanges(data.frame(chr=repeat_regions$chr, start=repeat_regions$start, end=repeat_regions$end))
  similar.data <- toGRanges(data.frame(chr=similar_regions$chr, start=similar_regions$start, end=similar_regions$end))
  
  # Add blacklist regions
  kp.repeat <- kpRect(kp, data = repeat.data, y0=0, y1=1, col="#FFDDDD", border=NA) # red
  kp.similar <- kpRect(kp, data = similar.data, y0=0, y1=1, col="#FFDDDD", border=NA) # green #DDFFDD
  
  return(kp.similar)
  
}


#---- Alternative: use BAM input file ----

# Set the filename of the BAM alignment from mouse mapping result (after filtering)
bam.149 <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/1.mapping_cleaning/output/mouse/149/149_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam"
bam.150 <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/1.mapping_cleaning/output/mouse/150/150_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam"
bam.151 <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/1.mapping_cleaning/output/mouse/151/151_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam"
bam.152 <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/1.mapping_cleaning/output/mouse/152/152_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam"
bam.153 <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/1.mapping_cleaning/output/mouse/153/153_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam"
bam.154 <- "/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/1.mapping_cleaning/output/mouse/154/154_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam"

# Density (here: 10 kb window)
window.size <- 1e4

# Plot the ideogram with the mouse reference genome
kp <- plotKaryotype(genome="mm10")

# Add blacklist regions
## ERROR! Only repeat elements are shown and not similar regions
kp <- karyoploteR_blacklist_regions_plot(workdir, blacklist_regions[2], blacklist_regions[1])

# Add coverage data
#kp <- kpPlotBAMDensity(kp, data = bam.149, window.size = window.size, r0=0, r1=0.2, ymax=50000, col=color[1])
#kp <- kpPlotBAMDensity(kp, data = bam.152, window.size = window.size, r0=0, r1=0.2, ymax=50000, col=color[4])
#kp <- kpPlotBAMDensity(kp, data = bam.150, window.size = window.size, r0=0, r1=0.2, ymax=50000, col=color[2])
#kp <- kpPlotBAMDensity(kp, data = bam.151, window.size = window.size, r0=0, r1=0.2, ymax=50000, col=color[3])
#kp <- kpPlotBAMDensity(kp, data = bam.153, window.size = window.size, r0=0, r1=0.2, ymax=50000, col=color[5])
#kp <- kpPlotBAMDensity(kp, data = bam.154, window.size = window.size, r0=0, r1=0.2, ymax=50000, col=color[9])

