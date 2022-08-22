#!/bin/bash

# Version 1.0
# This script was last modified on 05/11/2021 (UPDATE! Adapt from IS analysis on chimeric reads)
# Script Goal: Integration site analysis from "chimeric" molecules
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After filtering reads to keep only "chimeric" molecules on the alignment file, identify potential integration sites
# Extract chromosome name, start, end, length, barcode sequence and number of reads from the molecules recovered files on the AAV and mouse genome
# Merge the informations on the two genomes for each "chimeric" molecule
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file (with only "chimeric" molecules) of interest in quotes
# 3. Write the name of the mouse mapping BAM file (with only "chimeric" molecules) of interest in quotes
#
# Output: The script generates two lists
# 1. List of non-redundant integration sites on the mouse genome
# 2. List of regions of the AAV genome involved in these potential integration events
# 3. Merge the two lists to have information on the "chimeric" molecule from both genomes
# 
# Concept:
# Extract information on both molecules files from "chimeric" molecules and identify non-redundant integration sites and regions of the AAV genome involved
# 
# Note:



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/5.integration_site_analysis/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
sample=149 #------------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/5.integration_site_analysis/output/molecules

# Give the file name of the AAV/mouse mapping BAM file (with only "chimeric" molecules)
AAV_molecules_file="${workdir}/4.identify_chimeric_molecules/output/${sample}/${sample}_mapping_AAV.molecules.chimeric.bed"
mouse_molecules_file="${workdir}/4.identify_chimeric_molecules/output/${sample}/${sample}_mapping_mouse.molecules.chimeric.bed"

# Set the number of threads to use
nb_threads=10

# Create a folder to store all the results for this sample
mkdir -p $workdir/5.integration_site_analysis/output/molecules/${sample}

#----------------- END: No changes are needed below this line --------------------#



##################################################################################
######                             Main Script                              ######
##################################################################################

#---------------------- Merge the informations from the AAV and mouse chimeric molecules files ----------------------#

# Remove the header
sed 1d $AAV_molecules_file > $workdir/5.integration_site_analysis/output/molecules/${sample}/file1.tmp
sed 1d $mouse_molecules_file > $workdir/5.integration_site_analysis/output/molecules/${sample}/file2.tmp

# Sort the files before join them (by barcode sequence column)
sort -k 6 $workdir/5.integration_site_analysis/output/molecules/${sample}/file1.tmp > $workdir/5.integration_site_analysis/output/molecules/${sample}/file1_sorted.tmp
sort -k 6 $workdir/5.integration_site_analysis/output/molecules/${sample}/file2.tmp > $workdir/5.integration_site_analysis/output/molecules/${sample}/file2_sorted.tmp

# Change delimiter of a file: tab '\t' to comma ','
sed -i 's/\t/,/g' $workdir/5.integration_site_analysis/output/molecules/${sample}/file1_sorted.tmp
sed -i 's/\t/,/g' $workdir/5.integration_site_analysis/output/molecules/${sample}/file2_sorted.tmp


# Join lines of two files on a common field (barcode sequence) with awk
## It matches exactly the key column from both files
### Source: https://unix.stackexchange.com/questions/113898/how-to-merge-two-files-based-on-the-matching-of-two-columns
awk -F , 'NR==FNR {h[$5] = $0; next} {print $0, h[$5]}' $workdir/5.integration_site_analysis/output/molecules/${sample}/file2_sorted.tmp $workdir/5.integration_site_analysis/output/molecules/${sample}/file1_sorted.tmp > $workdir/5.integration_site_analysis/output/molecules/${sample}/file3_sorted.tmp

# Change delimiter of the file because there is tab and comma: convert all to tab delimiter
sed -i 's/   /\t/g' $workdir/5.integration_site_analysis/output/molecules/${sample}/file3_sorted.tmp
sed -i 's/ /\t/g' $workdir/5.integration_site_analysis/output/molecules/${sample}/file3_sorted.tmp
sed -i 's/,/\t/g' $workdir/5.integration_site_analysis/output/molecules/${sample}/file3_sorted.tmp

# Set header names: AAV_chromosome\tAAV_start\tAAV_end\tAAV_bx_sequence\tAAV_nb_reads\tAAV_length\tmouse_chromosome\tmouse_start\tmouse_end\tmouse_bx_sequence\tmouse_nb_reads\tmouse_length
sed '1i AAV_chromosome\tAAV_start\tAAV_end\tAAV_bx_sequence\tAAV_nb_reads\tAAV_length\tmouse_chromosome\tmouse_start\tmouse_end\tmouse_bx_sequence\tmouse_nb_reads\tmouse_length' $workdir/5.integration_site_analysis/output/molecules/${sample}/file3_sorted.tmp > $workdir/5.integration_site_analysis/output/molecules/${sample}/${sample}_AAV_mouse_chimeric_molecules_IS.txt

# Clean
rm $workdir/5.integration_site_analysis/output/molecules/${sample}/file1.tmp
rm $workdir/5.integration_site_analysis/output/molecules/${sample}/file2.tmp
rm $workdir/5.integration_site_analysis/output/molecules/${sample}/file1_sorted.tmp
rm $workdir/5.integration_site_analysis/output/molecules/${sample}/file2_sorted.tmp
rm $workdir/5.integration_site_analysis/output/molecules/${sample}/file3_sorted.tmp



#------------------------------------------- END --------------------------------------------------# 