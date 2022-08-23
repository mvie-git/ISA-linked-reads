#!/bin/bash

# Version 1.0
# This script was last modified on 12/11/2021
# Script Goal: Prepare annotation file for Rstudio script to plot graphs of the coverage per position
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: 
# Visualize regions of AAV vector sequence or mouse genome reference where these filtered reads mapped
# 
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file of interest in quotes
#
# Output: Pre-annotation file for Rstudio script


#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/visualization/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Coverage for vector or mouse
type=vector

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/visualization/coverage_per_position/input/${type}

# Set the number of threads to use
nb_threads=10


#----------------- END: No changes are needed below this line --------------------#



##################################################################################
######                             Main Script                              ######
##################################################################################

#----------------- Prepare annotation file for coverage graph --------------------#

# Create the output file to put the results from this script
# If doesn't exist
if [ ! -e "$workdir/visualization/coverage_per_position/input/${type}/coverage_per_position_all_sample.csv" ] ; then
    echo -e 'sample,chromosome,position,coverage' > "$workdir/visualization/coverage_per_position/input/${type}/coverage_per_position_all_sample.csv"
fi

## Looping through the content of a sample index list file: sample_list.txt
while read sample; do

  echo "<============ Starting sample n°$sample ============>"

    # Calculate coverage depth for all the position of the vector sequence or mouse genome (-a option include zero coverage positions)
    samtools depth -a $workdir/1.mapping_cleaning/output/${type}/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam > coverage_per_position_sample.tmp

    # Add which sample is processing in the first column & convert to csv format
    awk '{print '${sample}', $0}' coverage_per_position_sample.tmp | awk '{print $1 "," $2 "," $3 "," $4}' >> "$workdir/visualization/coverage_per_position/input/${type}/coverage_per_position_all_sample.csv"

done <$workdir/visualization/coverage_per_position/sample_list.txt

# Clean up!
rm *.tmp

    #------------------------------------------- END --------------------------------------------------#
