#!/bin/bash

### ERRATUM!!!

# Version 1.0
# This script was last modified on 04/11/2021
# Script Goal: Prepare blacklist regions file for chromoMap R package
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: visualize blacklist regions of AAV reference sequence for downstream analysis
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the BED file of interest in quotes containing the regions identified in AAV vector which are too similar to mouse genome
# 3. Write the name of the BED file of interest in quotes containing the regions identified in AAV vector to be repeat elements present in the mouse genome
#
# Output: a BED file with a specification of the type of 'blacklist' region (similar or repeat)



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/visualization/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/visualization/chromoMap/output

# Create a folder to store all the results for this script (only if it doesn't exist already!)
mkdir -p $workdir/visualization/chromoMap/input/vector

# Give the file name of the BED containing the similar regions to the mouse genome
similar_regions="${workdir}/1.mapping_cleaning/output/vector_similar_regions_to_remove_short.bed"

# Give the file name of the BED containing the repeat elements identified in the mouse genome
repeat_regions="${workdir}/1.mapping_cleaning/output/vector_repeat_regions_to_remove.bed"


#----------------- END: No changes are needed below this line --------------------#



##################################################################################
######                             Main Script                              ######
##################################################################################

#----------------- Prepare annotation file for chromoMap R package --------------------#

# Add a column to specify the type of blacklist region: similar/repeat and concatenate in one single file!
awk '{print "similar\t",$0}' $similar_regions > $workdir/visualization/chromoMap/input/vector/vector_similar_regions_to_remove_short.tmp
awk '{print "repeat\t",$0}' $repeat_regions > $workdir/visualization/chromoMap/input/vector/vector_repeat_regions_to_remove_short.tmp
cat $workdir/visualization/chromoMap/input/vector/vector_similar_regions_to_remove_short.tmp $workdir/visualization/chromoMap/input/vector/vector_repeat_regions_to_remove_short.tmp > $workdir/visualization/chromoMap/input/vector/vector_blacklist_regions_short.bed

# Clear!
rm $workdir/visualization/chromoMap/input/vector/*.tmp

# Launch R script "chromoMap_vector_blasklist_regions.R"

#------------------------------------------- END --------------------------------------------------#