#---- Import data ----

# Import library
library(dplyr)

# Data
data <- read.delim("~/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/similarity_vector_mouse_genome/BLAST/blast_output_vector_mouse.txt", header=FALSE)

# Rename the columns
names(data) <- c("title","chromosome","identity","?1","?2","?3","v_start","v_end","m_start","m_end","expect","score")



#---- Extract only vector start and end column ----

data <- data[c("v_start", "v_end")]

# Add a first column called "chr" to specify chr1 on vector
data$chr <- "chr1"

# Move the chr column to first position
data <- data %>%
  select("chr", everything())



#---- Export data to a bed file to filter reads on the alignment BAM file ----

write.table(data, file = "blast_output_vector_filter.bed", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)



