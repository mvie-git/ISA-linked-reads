# Source
#https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html

#---- Installation ----

# Manually: old version (working with Rstudio 3.6)
library(tidyr)
library(dplyr)
library(chromoMap)


#---- Input data file: format ----

# Import the chromosome file <------ CHANGE IT!
chromosome <- read.delim("/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/chromoMap/input/reference/chromosome_vector.txt", header=FALSE)


# Prepare the annotation file
## Import the data extracted from BAM alignment file <------ CHANGE IT!
data <- read.delim("/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA/visualization/chromoMap/input/vector/vector_blacklist_regions_short.bed", header=FALSE)

# Subset data
#data <- data[data$V1=="repeat", ]

# Replace "chr" term with none
data$V2 <- gsub(" chr", "", data$V2)

# Add a column ID
data$V5 <- paste0("region", seq.int(nrow(data)))

# Replace the optional data column with group at the end of the table
annotation <- data %>% 
  select(-V1, everything())

annotation <- annotation %>%
  select(V5, everything())

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

# Annotations are categorized into groups (similar or repeat 'blacklist' regions)
chromoMap("input/chromosome.txt", "input/annotation.txt", 
          segment_annotation = T,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("red","green")),
          vertical_grid = T,
          legend=T, 
          lg_x = 100,
          lg_y = 325)
