# Source
#https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html

#---- Installation ----

# Manually: old version (working with Rstudio 3.6)
library(tidyr)
library(dplyr)
library(chromoMap)


#---- TO CHANGE! Input data file: format ----

# Import the chromosome file <------ CHANGE IT!
chromosome <- read.delim("/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/chromoMap/input/reference/chromosome_vector.txt", header=FALSE)


# Prepare the annotation file
## Import the data extracted from BAM alignment file <------ CHANGE IT!
data <- read.delim("/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/chromoMap/input/vector/153/153_pre_annotation_file.txt", header=FALSE)

# Replace "chr" term with none
data$V3 <- gsub("chr", "", data$V3)

# Keep only the useful columns
annotation = subset(data, select = c(1, 2, 3, 4))

# Calcul the end position
annotation$V5 <- annotation$V4 + annotation$V2

# Remove the length column
annotation$V2 <- NULL

# Clean the first table
rm(data)






#---- Input data file: check ----

# Check the head of both files
## Tab-delimited
head(chromosome, sep = "\t")
head(annotation, sep = "\t")


#---- Input data file: export the input files ----

write.table(x = chromosome, file ="input/chromosome.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.table(x = annotation, file ="input/annotation.txt", sep="\t", col.names=FALSE, row.names=FALSE)


#---- ChromoMap ----

# Passing files

# Point annotation annotates an element on a single locus, ignoring its size (overlapping elements can be seen here)
chromoMap("input/chromosome.txt", "input/annotation.txt")

# Segment annotation consider the size and visualize the annotation as a segment (display gene structure)
chromoMap("input/chromosome.txt", "input/annotation.txt", 
          segment_annotation = T)

# Add vertical grid lines to separate regions of the AAV vector (without text)
chromoMap("input/chromosome.txt", "input/annotation.txt", 
          ref_line = T,
          refl_pos = 20,
          vertical_grid = T,
          grid_array = c(130, 184, 2883, 2958, 3398, 3412, 4485, 4513, 4675, 4729, 4856),
          segment_annotation = T)
