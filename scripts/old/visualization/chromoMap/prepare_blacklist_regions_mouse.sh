#!/bin/bash

# Version 1.0
# This script was last modified on 04/11/2021
# Script Goal: Prepare blacklist regions file for chromoMap R package
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: visualize blacklist regions of mouse reference sequence for downstream analysis
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the BED file of interest in quotes containing the regions identified in mouse genome which are too similar to AAV mouse sequence
# 3. Write the name of the BED file of interest in quotes containing the regions identified in mouse genome to be repeat elements
#
# Output: a BED file with a specification of the type of 'blacklist' region (similar or repeat)



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/visualization/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/visualization/chromoMap/output

# Create a folder to store all the results for this script (only if it doesn't exist already!)
mkdir -p $workdir/visualization/chromoMap/input/mouse

# Give the file name of the BED containing the similar regions to the AAV mouse sequence
similar_regions="${workdir}/1.mapping_cleaning/output/mouse_similar_regions_to_remove.bed"

# Give the file name of the BED containing the repeat elements identified in the mouse genome
repeat_regions="${workdir}/1.mapping_cleaning/output/mouse_repeat_regions_to_remove.bed"


#----------------- END: No changes are needed below this line --------------------#



##################################################################################
######                             Main Script                              ######
##################################################################################

#----------------- Prepare annotation file for chromoMap R package --------------------#

# Add a column to specify the type of blacklist region: similar/repeat and concatenate in one single file!
awk '{print "similar\t",$0}' $similar_regions > $workdir/visualization/chromoMap/input/mouse/mouse_similar_regions_to_remove.tmp
awk '{print "repeat\t",$0}' $repeat_regions > $workdir/visualization/chromoMap/input/mouse/mouse_repeat_regions_to_remove.tmp
cat $workdir/visualization/chromoMap/input/mouse/mouse_similar_regions_to_remove.tmp $workdir/visualization/chromoMap/input/mouse/mouse_repeat_regions_to_remove.tmp > $workdir/visualization/chromoMap/input/mouse/mouse_blacklist_regions.bed

# Clear!
rm $workdir/visualization/chromoMap/input/mouse/*.tmp

# Launch R script "chromoMap_mouse_blasklist_regions.R"

#------------------------------------------- END --------------------------------------------------#