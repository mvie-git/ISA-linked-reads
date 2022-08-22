#!/bin/bash

# Version 1.0
# This script was last modified on 19/11/2021 (TO DO: review this version to change it accordly to molecules!!!)
# Script Goal: Structure analysis of the potential "chimeric" molecules identified previously...
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After identifying unique "chimeric" molecules, blast them against both genomes (AAV/mouse) to check if there are molecules with reads map partially to each genome
# Create a BLAST database locally for each genome (AAV and mouse)
# Merge the two BLAST datablases together to search in multiple databases at the same time
# Modify the file containing the list of potential "chimeric" molecules in a FASTA format
# Align those potential "chimeric" molecules against the hybrid genome (AAV/mouse) reference
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the file containing the list of potential "chimeric" molecules of interest
# 3. Write the name of the AAV vector sequence fasta file
# 4. Write the name of the mouse genome reference fasta file
#
# Output: The script generates an alignment output from the BLAST command line
# 
# Concept:
# BLAST the "chimeric" molecules against the two genomes (AAV/mouse) to identify if there are reads mapped to mouse genome and AAV vector!
# 
# Note:



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/6.integration_structure_analysis/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/6.integration_structure_analysis/output

# Give the name of the AAV vector sequence fasta file
vector_reference="${workdir}/0.preparing_input_files/vector_genome.fasta"

# Give the name of the mouse genome reference fasta file
mouse_reference="${workdir}/0.preparing_input_files/mouse_genome.fa"

# Set the number of threads to use
nb_threads=10


#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

    echo "*** Start $sample ***"

    # Create a folder to store all the results for this sample
    mkdir -p $workdir/6.integration_structure_analysis/output/${sample}

    # Give the name of the file containing the list of potential "chimeric" reads of interest
    potential_chimeric_molecules="${workdir}/4.identify_chimeric_molecules/output/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt" #----------------> !Change IT!#


    #----------------- END: No changes are needed below this line --------------------#



    ##################################################################################
    ######                             Main Script                              ######
    ##################################################################################

    #----------------- Align those potential "chimeric" molecules against the hybrid genome (AAV/mouse) reference --------------------#

    # Extract reads from AAV and mouse BAM files with those barcodes
    

    # Add an unique ID for each potential "chimeric" read
    sed 1d $potential_chimeric_reads | awk 'BEGIN{FS=OFS="\t"} {print "read" ++count, $0}' > ${workdir}/6.integration_structure_analysis/output/${sample}/${sample}_reads_chimeric_ID.txt

    # Modify the file containing the list of potential "chimeric" reads in a FASTA format
    sed 1d $potential_chimeric_reads | awk '{print $7}' | awk '{print ">read" ++count ORS $0}' > ${workdir}/6.integration_structure_analysis/output/${sample}/${sample}_reads_chimeric.fasta

    # Align this list of potential chimeric reads against the hybrid genome
    ## Tool: blastn
    ### Options:
    ### -db: BLAST database name
    ### -query: input file => "chimeric" reads to align against AAV/mouse genome BLAST database
    ### -outfmt: formatting options (alignment view options): 6 (tabular), 7 (tabular with comment lines)
    ../bin/blastn -db blast_db/hybrid_genome -query ${workdir}/6.integration_structure_analysis/output/${sample}/${sample}_reads_chimeric.fasta -outfmt 6 > ${workdir}/6.integration_structure_analysis/output/${sample}/${sample}_blast_output.txt

    ### Options: 
    ### -blastn-short: optimized for sequences less than 30 nucleotides
    ../bin/blastn -db blast_db/hybrid_genome -query ${workdir}/6.integration_structure_analysis/output/${sample}/${sample}_reads_chimeric.fasta -outfmt 6 -task "blastn-short" > ${workdir}/6.integration_structure_analysis/output/${sample}/${sample}_blast_short_output.txt


    #------------------------------------------- END --------------------------------------------------# 

done < sample_list.txt