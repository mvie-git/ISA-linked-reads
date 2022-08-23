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
workdir <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/AAV_ITRs/3_ITR/"
workdir_blacklist_regions <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/"
AAV_genome <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/genome/3_ITR_vector_genome.txt"
AAV_cytobands <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/genome/3_ITR_vector_cytobands.txt"
sample <- list.files(path="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/AAV_ITRs/3_ITR/", pattern="*_vector_coverage_position.txt")
blacklist_regions <- list.files(path="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/", pattern="vector_blacklist_regions_*")
color <- brewer.pal(n = 9, name = "Set1")
# Display all the colors
display.brewer.all()


#---- Create a function to plot an ideogram ----

# Create a function to plot the ideogram using a custom genome
karyoploteR_plot <- function(workdir, genome, cytobands) {
  
  # For AAV genome, create a custom genome ideogram
  custom.genome <- toGRanges(genome)
  
  # Provide cytobands information to represent the different regions of the vector
  ## Options:
  ### name: ID of each cytoband information
  ### gieStain: reference names of each region of interest
  custom.cytobands <- toGRanges(cytobands)

  # Plot the ideogram using the custom genome
  kp <- plotKaryotype(plot.type=2, genome = custom.genome, cytobands = custom.cytobands)
  
  # Add the cytoband names to the ideogram
  kp <- kpAddCytobandLabels(kp, cex=0.5)
  
  return(kp)
  
}


#---- Create a function to add coverage data ----

# Create a function to add coverage data on the ideogram plot
karyoploteR_coverage_plot <- function(workdir, color, path_to_data, pos0, pos1, sample_name) {
  
  # Import coverage data
  data <- read.delim(paste0(workdir, path_to_data), header=FALSE)
  
  # Create a GRanges object to store coverage data
  coverage.data <- toGRanges(data.frame(chr="1", start=data[, "V2"], end=data[, "V2"]+1, y=data[, "V3"]))
  
  # Add coverage data
  kp.coverage <- kpBars(kp, coverage.data,  y1=coverage.data$y, ymax=max(coverage.data$y), col=color, border=NA, r0=pos0, r1=pos1)
  
  # Add labels (sample name)
  #kp.IS.sample <- kpAddLabels(kp.coverage, labels=paste0("sample ", sample_name), r0=pos0, r1=pos1, data.panel=1, cex=0.9, col=color, side="right")
  
  # Add axis
  kp.coverage.axis <- kpAxis(kp.coverage, ymin=0, ymax=max(coverage.data$y), r0=pos0, r1=pos1, cex=0.7)
  
  return(kp.coverage.axis)
  
}


#---- Create a function to add blacklist regions ----

# Create a function to add coverage data on the ideogram plot
karyoploteR_blacklist_regions_plot <- function(workdir, path_to_repeat_regions, path_to_similar_regions, ymax) {
  
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

### TO DO: reinitialize each time the ideogram plot before each sample!!!
# Execute the function to plot an AAV custom genome ideogram
kp <- karyoploteR_plot(workdir, AAV_genome, AAV_cytobands)

# Add blacklist regions
## Add ymax
kp <- karyoploteR_blacklist_regions_plot(workdir_blacklist_regions, blacklist_regions[1], blacklist_regions[2], 1.25)

# Plot additional data: coverage and blacklist regions for each sample
## Add positions for each sample
karyoploteR_coverage_plot(workdir, color[1], sample[1], 0, 0.20, "149") # 149
karyoploteR_coverage_plot(workdir, color[4], sample[4], 0.25, 0.45, "152") # 152
karyoploteR_coverage_plot(workdir, color[2], sample[2], 0.50, 0.70, "150") # 150
karyoploteR_coverage_plot(workdir, color[3], sample[3], 0.75, 0.95, "151") # 151
karyoploteR_coverage_plot(workdir, color[5], sample[5], 1.0, 1.20, "153") # 153
karyoploteR_coverage_plot(workdir, color[9], sample[6], 1.25, 1.45, "154") # 154


