#!/bin/bash

# Version 1.0
# This script was last modified on 07/12/2021
# Script Goal: Local assembly
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After identifying potential "chimeric" molecules, assemble the corresponding mouse locus
# Extract list of "chimeric" barcodes meaning present in each genome (AAV/mouse)
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the mouse molecule BED file with only chimeric barcodes in quotes
# 
# Concept:
# 
# Note:



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/2.identify_chimeric_reads/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
#sample=149 #------------------------> !CHANGE IT!#

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create the output directory (only if it doesn't exist already!)
	mkdir -p $workdir/7.local_assembly/output/short

	# Give the file name of the mouse molecule BED file with only chimeric barcodes
	mouse_molecule_file_chimeric="$workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.chimeric.bed"

	# Give the file name of the mouse mapping BAM file with only chimeric barcodes
	mouse_mapping_file_chimeric="$workdir/4.identify_chimeric_molecules/output/short/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam"


	# Set the number of threads to use
	nb_threads=10

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/7.local_assembly/output/short/${sample}

	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################

	#----------------- Extract list of "chimeric" barcodes meaning present in each genome (AAV/mouse) --------------------#

	# Extract coordinates for genomic regions of "chimeric" molecules
	awk '{print $1 ":" $2 "-" $3}' $mouse_molecule_file_chimeric > "$workdir/7.local_assembly/output/short/${sample}_mapping_mouse.molecules.chimeric.txt"

	# 

	# Compteur
	cpt=0

	while read molecule; do 

		# Extract reads from a BAM file that falls within a region of interest and extract barcodes associated with those reads
		../bin/samtools view $mouse_mapping_file_chimeric $molecule | grep -o "BX:Z:[ACTG]*" > "${workdir}/7.local_assembly/output/short/${sample}/${sample}_list_barcodes.tmp";
	
		# Subset the BAM file to keep only reads associated with these barcodes
		# INFO! If you doesn't want to decompress the BAM into SAM file, you need to extract the header of the samtools view command (SAM file) to put it back before convert it back to BAM 			format

		# Save the header lines 
		../bin/samtools view -H $mouse_mapping_file_chimeric > "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_header_mapping.sam"

		# Extract mapped reads from the BAM file with the list of chosen barcodes
		## Filter alignments using list of chosen barcodes (-w option for exact string match)
		../bin/samtools view -@ $nb_threads $mouse_mapping_file_chimeric | grep -w -f "${workdir}/7.local_assembly/output/short/${sample}/${sample}_list_barcodes.tmp" > "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_body_mapping.sam"

		# Combine header and body 
		cat "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_header_mapping.sam" "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_body_mapping.sam" > "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.sam" 

		# Convert output.subset.sam to BAM format 
		../bin/samtools view -@ $nb_threads -b "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.sam"  > "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.bam"

		# Sort and index the BAM file 
		../bin/samtools sort -@ $nb_threads "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.bam" -o "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.subset.bam"
		
		# R! Sort the BAM file by barcode?
		../bin/samtools index "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.subset.bam"

		# Clean up!
		rm "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_header_mapping.sam"
		rm "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_body_mapping.sam"
		rm "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.sam"
		rm "${workdir}/7.local_assembly/output/short/${sample}/${sample}_output_all_mapping.bam"


	done < "$workdir/7.local_assembly/output/short/${sample}_mapping_mouse.molecules.chimeric.txt"

	#------------------------------------------- END --------------------------------------------------# 

done < sample_list.txt