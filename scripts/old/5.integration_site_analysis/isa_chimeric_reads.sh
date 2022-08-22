#!/bin/bash

# Version 1.0
# This script was last modified on 27/10/2021 (UPDATE! Calcul of the end position of each 'chimeric' read)
# Script Goal: Integration site analysis from "chimeric" reads
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After filtering reads to keep only "chimeric" reads on the alignment file, identify potential integration sites
# Extract chromosome name, start, length, read sequence from the alignment files on the AAV and mouse genome
# Calculate the end position
# Collapse all the "chimeric" reads containing the same insertion sites
# Merge the informations on the two genomes for each "chimeric" read
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file (with only "chimeric" reads) of interest in quotes
# 3. Write the name of the mouse mapping BAM file (with only "chimeric" reads) of interest in quotes
#
# Output: The script generates two lists
# 1. List of redundant integration sites on the mouse genome
# 2. List of regions of the AAV genome involved in these potential integration events
# 3. Merge the two lists to have information on the "chimeric" read from both genomes
# 
# Concept:
# Extract information on both alignment files from "chimeric" reads and identify non-redundant integration sites and regions of the AAV genome involved
# 
# Note:



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/5.integration_site_analysis/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
#sample=151 #------------------------> !CHANGE IT!#

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

    echo "*** Start $sample ***"

    # Create the output directory (only if it doesn't exist already!)
    mkdir -p $workdir/5.integration_site_analysis/output

    # Give the file name of the AAV/mouse mapping BAM file (with only "chimeric" reads)
    AAV_mapping_file="${workdir}/2.identify_chimeric_reads/output/short/${sample}/${sample}_AAV_mapping_filtered.chimeric.bam"
    mouse_mapping_file="${workdir}/2.identify_chimeric_reads/output/short/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam"

    # Set the number of threads to use
    nb_threads=10

    # Create a folder to store all the results for this sample
    mkdir -p $workdir/5.integration_site_analysis/output/short/${sample}

    #----------------- END: No changes are needed below this line --------------------#



    ##################################################################################
    ######                             Main Script                              ######
    ##################################################################################

    #----------------- Extract the important informations from the AAV chimeric reads alignment file --------------------#

    # Create the output file to put the results from this script
    # If doesn't exist
    if [ ! -e "$workdir/5.integration_site_analysis/output/short/${sample}/${sample}_AAV_chimeric_IS.txt" ] ; then
        echo -e 'name\tchromosome\tstart\tMAPQ\tsequence' > "$workdir/5.integration_site_analysis/output/short/${sample}/${sample}_AAV_chimeric_IS.txt"
    fi

    # Extract chromosome name, start, length, read sequence from the alignment files on the AAV and mouse genome
    ../bin/samtools view -@ $nb_threads $AAV_mapping_file | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $10}' >> $workdir/5.integration_site_analysis/output/short/${sample}/${sample}_AAV_chimeric_IS.txt

      

    #----------------- Extract the important informations from the mouse chimeric reads alignment file --------------------#

    # Create the output file to put the results from this script
    # If doesn't exist
    if [ ! -e "$workdir/5.integration_site_analysis/output/short/${sample}/${sample}_mouse_chimeric_IS.txt" ] ; then
        echo -e 'name\tchromosome\tstart\tMAPQ\tsequence' > "$workdir/5.integration_site_analysis/output/short/${sample}/${sample}_mouse_chimeric_IS.txt"
    fi

    # Extract chromosome name, start, length, read sequence from the alignment files on the AAV and mouse genome
    ../bin/samtools view -@ $nb_threads $mouse_mapping_file | awk '{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $10}' >> $workdir/5.integration_site_analysis/output/short/${sample}/${sample}_mouse_chimeric_IS.txt



    #---------------------- Collapse all the "chimeric" reads containing the same insertion sites ----------------------#

    # Same insertion site: same chromosome and same start position?



    #---------------------- Merge the informations from the AAV and mouse chimeric reads alignment files ----------------------#

    # Remove the header
    sed 1d $workdir/5.integration_site_analysis/output/short/${sample}/${sample}_AAV_chimeric_IS.txt > $workdir/5.integration_site_analysis/output/short/${sample}/file1.tmp
    sed 1d $workdir/5.integration_site_analysis/output/short/${sample}/${sample}_mouse_chimeric_IS.txt > $workdir/5.integration_site_analysis/output/short/${sample}/file2.tmp

    # Sort the files before join them (by read sequence column)
    sort -k 6 $workdir/5.integration_site_analysis/output/short/${sample}/file1.tmp > $workdir/5.integration_site_analysis/output/short/${sample}/file1_sorted.tmp
    sort -k 6 $workdir/5.integration_site_analysis/output/short/${sample}/file2.tmp > $workdir/5.integration_site_analysis/output/short/${sample}/file2_sorted.tmp

    # Change delimiter of a file: tab '\t' to comma ','
    sed -i 's/\t/,/g' $workdir/5.integration_site_analysis/output/short/${sample}/file1_sorted.tmp
    sed -i 's/\t/,/g' $workdir/5.integration_site_analysis/output/short/${sample}/file2_sorted.tmp


    # Join lines of two files on a common field (read sequence) with awk
    ## It matches exactly the key column from both files
    ### Source: https://unix.stackexchange.com/questions/113898/how-to-merge-two-files-based-on-the-matching-of-two-columns
    awk -F , 'NR==FNR {h[$5] = $0; next} {print $0, h[$5]}' $workdir/5.integration_site_analysis/output/short/${sample}/file2_sorted.tmp $workdir/5.integration_site_analysis/output/short/${sample}/file1_sorted.tmp > $workdir/5.integration_site_analysis/output/short/${sample}/file3_sorted.tmp

    # Change delimiter of the file because there is tab and comma: convert all to tab delimiter
    sed -i 's/   /\t/g' $workdir/5.integration_site_analysis/output/short/${sample}/file3_sorted.tmp
    sed -i 's/ /\t/g' $workdir/5.integration_site_analysis/output/short/${sample}/file3_sorted.tmp
    sed -i 's/,/\t/g' $workdir/5.integration_site_analysis/output/short/${sample}/file3_sorted.tmp

    # Calculate the end position of each read in the AAV and mouse genome
    awk 'NR>1{print $1 "\t" $2 "\t" $3 "\t" $3+length($5) "\t" length($5) "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $8+length($10) "\t" length($10) "\t" $9 "\t" $10}' $workdir/5.integration_site_analysis/output/short/${sample}/file3_sorted.tmp > $workdir/5.integration_site_analysis/output/short/${sample}/file4_sorted.tmp

    # Gather unique chimeric alignment
    cat $workdir/5.integration_site_analysis/output/short/${sample}/file4_sorted.tmp | sort | uniq > $workdir/5.integration_site_analysis/output/short/${sample}/file4_sorted_uniq.tmp

    # Set header names: AAV_name\tAAV_chromosome\tAAV_start\tAAV_end\tAAV_MAPQ\tAAV_sequence\tmouse_name\tmouse_chromosome\tmouse_start\tmouse_end\tmouse_MAPQ\tmouse_sequence
    sed '1i AAV_name\tAAV_chromosome\tAAV_start\tAAV_end\tAAV_length\tAAV_MAPQ\tAAV_sequence\tmouse_name\tmouse_chromosome\tmouse_start\tmouse_end\tmouse_length\tmouse_MAPQ\tmouse_sequence' $workdir/5.integration_site_analysis/output/short/${sample}/file4_sorted_uniq.tmp > $workdir/5.integration_site_analysis/output/short/${sample}/${sample}_AAV_mouse_chimeric_IS.txt

    # Clean
    rm $workdir/5.integration_site_analysis/output/short/${sample}/*.tmp


    #------------------------------------------- END --------------------------------------------------# 

done < sample_list.txt