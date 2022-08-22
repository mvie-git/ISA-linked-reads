#!/bin/bash

# Version 1.0
# This script was last modified on 01/12/2021
# Script Goal: Chimeric molecules
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After filtering the output from AAV and mouse mapping results, identify which barcodes are "chimeric" (containing reads mapped either to rAAV sequence or mouse genome)
# Extract list of unique barcodes present in each genome (AAV/mouse)
# Compare the two lists in order to identify "chimeric" barcodes
# Subset the AAV/mouse mapping BAM file (filtered) in order to keep only information on chimeric molecules
# Subset the AAV/mouse molecule BED file in order to keep only information on chimeric molecules
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file (filtered) of interest in quotes
# 3. Write the name of the mouse mapping BAM file (filtered) of interest in quotes
# 4. Write the name of the AAV molecule BED file of interest in quotes
# 5. Write the name of the mouse molecule BED file of interest in quotes
#
# Output: The script generates two BAM files
# 1. List of unique reads present in both genomes (AAV/mouse) = "chimeric" reads
# 2. Subset of the AAV mapping BAM file (filtered) with only "chimeric" reads
# 3. Subset of the mouse mapping BAM file (filtered) with only "chimeric" reads
# 4. Subset of the AAV molecule BED file with only "chimeric" barcodes
# 5. Subset of the mouse molecule BED file with only "chimeric" barcodes
# 
# Concept:
# Extract unique barcodes from AAV molecule BED file and compare it to the list of unique barcodes from mouse molecule BED file
# Extract unique barcodes from AAV mapping BAM file (filtered) and compare it to the list of unique barcodes from mouse mapping BAM file (filtered)
# 
# Note:
# Files from the mouse mapping are heavy!



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/2.identify_chimeric_reads/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
#sample=149 #------------------------> !CHANGE IT!#

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create the output directory (only if it doesn't exist already!)
	mkdir -p $workdir/4.identify_chimeric_molecules/output/short

	# Give the file name of the AAV/mouse mapping BAM file (filtered)
	AAV_mapping_file="${workdir}/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam"
	mouse_mapping_file="${workdir}/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam"

	# Give the file name of the AAV/mouse molecule BED file
	AAV_molecule_file="${workdir}/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.bed"
	mouse_molecule_file="${workdir}/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.bed.gz"

	# Set the number of threads to use
	nb_threads=10

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/4.identify_chimeric_molecules/output/short/${sample}

	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################

	#----------------- Extract list of unique barcodes present in both genomes (AAV/mouse) = "chimeric" barcodes --------------------#
	
	# Extract list of unique barcodes mapped to AAV after filtering false positives
	sed 1d $AAV_molecule_file | awk '{print $4}' | sort | uniq > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_barcodes_after_filtering.txt
	wc -l $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_barcodes_after_filtering.txt

	# Extract list of unique reads mapped to mouse after filtering false positives
	zcat $mouse_molecule_file | sed 1d | awk '{print $4}' | sort | uniq > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_mouse_barcodes_after_filtering.txt
	wc -l $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_mouse_barcodes_after_filtering.txt

	# Extract barcodes present in only AAV molecule BED file and barcodes common in both genomes (-2: supress column with barcodes only in mouse genome)
	comm -2 $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_barcodes_after_filtering.txt $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_mouse_barcodes_after_filtering.txt > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_after_filtering_comm_output.tmp

	# Extract the first column with barcodes only present in AAV molecule BED file (meaning AAV only barcodes)
	awk -F "\t" '{print $1}' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_after_filtering_comm_output.tmp > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_AAV_only.txt

	# Remove empty lines and count the number of barcodes with reads mapped only AAV genome (meaning AAV only barcodes)
	sed -i '/^$/d' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_AAV_only.txt 
	wc -l $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_AAV_only.txt 

	# Extract the second column with barcodes present in both molecules BED files (AAV and mouse)
	awk -F "\t" '{print $2}' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_after_filtering_comm_output.tmp > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt

	# Remove empty lines and count the number of barcodes with reads mapped to both genomes (AAV and mouse)
	sed -i '/^$/d' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt
	wc -l $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt

	# Clean up!
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_after_filtering_comm_output.tmp
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_barcodes_after_filtering.txt # Takes a lot of disk place
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_mouse_barcodes_after_filtering.txt # Takes a lot of disk place



	#-------- Subset the AAV mapping BAM file (filtered) in order to keep only information on chimeric barcodes --------#

	# INFO! If you doesn't want to decompress the BAM into SAM file, you need to extract the header of the samtools view command (SAM file) to put it back before convert it back to BAM format

	# Save the header lines 
	../bin/samtools view -H $AAV_mapping_file > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_header_AAV_mapping_filtered.sam

	# Extract mapped reads from the BAM filtered file on the AAV genome with the list of "chimeric" barcodes
	## Filter alignments using list of "chimeric" barcodes (-w option for exact string match)
	## Add a prefix to the barcode: BX:Z:
	sed -i -e 's/^/BX:Z:/' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt
	../bin/samtools view -@ $nb_threads $AAV_mapping_file | grep -w -f $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_body_AAV_mapping_filtered.sam
	sed -i -e 's/BX:Z://' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt

	# Combine header and body 
	cat $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_header_AAV_mapping_filtered.sam $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_body_AAV_mapping_filtered.sam > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.sam 

	# Convert output.filtered.sam to BAM format 
	../bin/samtools view -@ $nb_threads -b $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.sam  > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.bam

	# Sort and index the BAM file 
	../bin/samtools sort -@ $nb_threads $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.bam -o $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_AAV_mapping_filtered.chimeric.bam
	../bin/samtools index $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_AAV_mapping_filtered.chimeric.bam


	# Subset the AAV molecule BED file in order to keep only information on chimeric barcodes
	grep -w -f $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt ${workdir}/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.bed > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.chimeric.bed


	# Clean up!
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_header_AAV_mapping_filtered.sam
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_body_AAV_mapping_filtered.sam
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.sam
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.bam



	#-------- Subset the mouse mapping BAM file (filtered) in order to keep only information on chimeric barcodes --------#

	# INFO! If you doesn't want to decompress the BAM into SAM file, you need to extract the header of the samtools view command (SAM file) to put it back before convert it back to BAM format

	# Save the header lines 
	../bin/samtools view -H $mouse_mapping_file > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_header_mouse_mapping_filtered.sam

	# Extract mapped reads from the BAM filtered file on the mouse genome with the list of "chimeric" barcodes
	## Filter alignments using list of "chimeric" barcodes (-w option for exact string match)
	## Add a prefix to the barcode: BX:Z:
	sed -i -e 's/^/BX:Z:/' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt
	../bin/samtools view -@ $nb_threads $mouse_mapping_file | grep -w -f $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_body_mouse_mapping_filtered.sam
	sed -i -e 's/BX:Z://' $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt

	# Combine header and body 
	cat $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_header_mouse_mapping_filtered.sam $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_body_mouse_mapping_filtered.sam > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.sam 

	# Convert output.filtered.sam to BAM format 
	../bin/samtools view -@ $nb_threads -b $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.sam  > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.bam

	# Sort and index the BAM file 
	../bin/samtools sort -@ $nb_threads $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.bam -o $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam
	../bin/samtools index $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam


	# Subset the mouse molecule BED file in order to keep only information on chimeric barcodes
	zgrep -w -f $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_list_AAV_mouse_barcodes_chimeric.txt ${workdir}/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.bed > $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.chimeric.bed
	

	# Clean up!
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_header_mouse_mapping_filtered.sam
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_body_mouse_mapping_filtered.sam
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.sam
	rm $workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.bam


	# BONUS!
	# Filter out molecules (mouse regions) with only one read and a length lower than 1000
	# Calculate the density of reads per genomic region
	# Sort each genomic region by increasing density values
	cat 149_mapping_mouse.molecules.chimeric.bed | awk '{ if( $5>1 && $6>1000 ) print $0 }' | awk '$7=$5/$6 {print $0}'| sort -n -k7


	#------------------------------------------- END --------------------------------------------------# 

done < sample_list.txt
