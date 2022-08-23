#!/bin/bash

# Version 1.0
# This script was last modified on 04/11/2021
# Script Goal: Prepare annotation file for chromoMap R package on MOUSE genome
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: 
# Visualize regions of mouse reference genome where these filtered reads mapped OR
# Visualize regions of mouse reference genome which are potential integration sites
# 
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the mouse mapping BAM file of interest in quotes
#
# Output: Pre-annotation file for Rstudio script


#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/visualization/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
sample=154 #------------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/visualization/chromoMap/output

# Create a folder to store all the results for this sample
mkdir -p $workdir/visualization/chromoMap/input/mouse/${sample}

# Give the file name of the mouse mapping BAM file filtered (no false positives)
#mapping_file="${workdir}/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam"
# Give the file name of the mouse mapping BAM file filtered (in order to keep only information on chimeric reads)
mapping_file="$workdir/2.identify_chimeric_reads/output/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam"


# Set the number of threads to use
nb_threads=10


#----------------- END: No changes are needed below this line --------------------#



##################################################################################
######                             Main Script                              ######
##################################################################################

#----------------- Prepare annotation file for chromoMap R package --------------------#

# Select columns of interest (read name, chr name, start, length) from the mouse mapping file after filtering
../../bin/samtools view -@ $nb_threads $mapping_file | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > $workdir/visualization/chromoMap/input/mouse/${sample}/${sample}_pre_annotation_file.txt

# Launch R script "chromoMap_mouse_IS.R"

#------------------------------------------- END --------------------------------------------------#