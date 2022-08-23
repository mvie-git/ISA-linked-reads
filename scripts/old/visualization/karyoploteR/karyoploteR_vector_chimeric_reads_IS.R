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
workdir_data <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/IS_chimeric_molecules_reads/"
workdir_blacklist <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/"
AAV_genome <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/genome/AAV_vector_genome.txt"
AAV_cytobands <- "/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/genome/AAV_vector_cytobands.txt"
sample <- list.files(path="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/visualization/karyoploteR/input/vector/IS_chimeric_molecules_reads/", pattern="*_AAV_chimeric_IS.bed")
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


#---- Create a function to add regions of the AAV vector involved in IS ----

# Create a function to add regions of the AAV vector involved in IS on the ideogram plot
karyoploteR_IS_plot <- function(workdir, color, path_to_data, pos0, pos1, sample_name) {
  
  # Import IS data
  data <- read.delim(paste0(workdir, path_to_data), header=TRUE)
 
   # Replace "chr" term with none
  data$chr <- gsub("chr", "", data$chr)
  
  # Calcul the end position
  data$end <- data$start + data$length
  
  # Remove the length column
  data$length <- NULL
  
  # Create a GRanges object to store IS data
  IS.data <- toGRanges(data.frame(chr="1", start=data[, "start"], end=data[, "end"]))
  
  # Add IS data
  kp.IS <- kpPlotRegions(kp, data=IS.data, col=color, layer.margin = 0.1, border=NA, r0=pos0, r1=pos1)
  
  # Add labels (sample name)
  kp.IS.sample <- kpAddLabels(kp.IS, labels=paste0("sample ", sample_name), r0=pos0, r1=pos1, data.panel=1, cex=0.9, col=color, side="left")
  
  # Add labels (number of total reads per sample)
  kp.IS.sample.nb_reads <- kpAddLabels(kp.IS.sample, labels=nrow(data), r0=pos0, r1=pos1, data.panel=1, cex=0.9, col=color, side="right")
  
  return(kp.IS.sample.nb_reads)
  
}



#---- Create a function to add blacklist regions ----

# Create a function to add coverage data on the ideogram plot
karyoploteR_blacklist_regions_plot <- function(workdir,path_to_repeat_regions, path_to_similar_regions, ymax) {
  
  # Import blacklist regions data
  similar_regions <- read.delim(paste0(workdir, path_to_similar_regions))
  repeat_regions <- read.delim(paste0(workdir, path_to_repeat_regions))
  
  # Create a GRanges object to store blacklist regions informations
  similar.data <- toGRanges(data.frame(chr=similar_regions$chr, start=similar_regions$start, end=similar_regions$end))
  repeat.data <- toGRanges(data.frame(chr=repeat_regions$chr, start=repeat_regions$start, end=repeat_regions$end))
  
  # Add blacklist regions
  kp.similar <- kpRect(kp, data = similar.data, y0=0, y1=ymax, col="#DDFFDD", border=NA) # green
  kp.repeat <- kpRect(kp, data = repeat.data, y0=0, y1=ymax, col="#FFDDDD", border=NA) # red
  
  return(kp.similar)
  
}


#---- [Final] Plot the ideogram with additional data ----

### TO DO: reinitialize each time the ideogram plot before each sample!!!
# Execute the function to plot an AAV custom genome ideogram
kp <- karyoploteR_plot(workdir, AAV_genome, AAV_cytobands)

# Add blacklist regions
## Add ymax
## TO DO! Change workdir directory!
kp <- karyoploteR_blacklist_regions_plot(workdir_blacklist, blacklist_regions[1], blacklist_regions[2], 1.5)

# Plot additional data: IS data for each sample
## Add regions for each sample
karyoploteR_IS_plot(workdir_data, color[1], sample[1], 0, 0.2, "149") # 149
karyoploteR_IS_plot(workdir_data, color[4], sample[4], 0.25, 0.45, "152") # 152
karyoploteR_IS_plot(workdir_data, color[2], sample[2], 0.50, 0.70, "150") # 150
karyoploteR_IS_plot(workdir_data, color[3], sample[3], 0.75, 0.95, "151") # 151
karyoploteR_IS_plot(workdir_data, color[5], sample[5], 1.0, 1.20, "153") # 153
karyoploteR_IS_plot(workdir_data, color[9], sample[6], 1.25, 1.45, "154") # 154

