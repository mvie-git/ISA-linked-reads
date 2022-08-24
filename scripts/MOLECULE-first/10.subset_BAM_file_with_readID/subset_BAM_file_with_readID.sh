# Subset a BAM file by a list of read IDs of interest

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"
nb_threads=12
sambamba=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/bin/sambamba


## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create an output directory
	mkdir -p ${sample}


	# Input files
	input_BAM=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/output/${sample}/${sample}_mapping_mouseANDvector.sorted.bam
	index_file=${workdir}/9.retrieve_all_readID_from_barcodes/${sample}/${sample}_subset_index_file_all_reads.txt

	# Extract only the readID list of the filtered index file
	awk '{print $1}' $index_file | sort | uniq > ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_readID_list.txt

	# Subset the BAM file with a list of readID of interest
	samtools view -@ $nb_threads $input_BAM | grep -w -f ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_readID_list.txt > ${sample}_subset_BAM.sam

	# Extract the header
	samtools view -H $input_BAM > header.sam

	# Concatenate header+body and convert it to BAM format
	cat header.sam ${sample}_subset_BAM.sam | samtools view -@ $nb_threads -Sb > ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.bam

	# Sort the new BAM file
	samtools sort ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.bam > ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.bam

	# Index the new BAM file
	samtools index ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.bam

	# Cleaning
	rm header.sam ${sample}_subset_BAM.sam ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.bam


	# (Optional) Filter out unmapped reads, duplicate reads, multi-map reads, MAPQ0
		# Filter out alignments which correspond to the FLAG option
	## Tool: samtools view
	## Options:
	### -b: output in the BAM format
	### -F: filter out alignments with all bits in INT present in the FLAG field
	### -@: number of threads
	
	### FLAG: 4 -> remove unmapped reads from the AAV mapping BAM file
	samtools view -@ $nb_threads -b -F 4 ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.bam > ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.mapped.bam

	# Identify and remove duplicate reads (PCR and optical duplicates) from the AAV mapping BAM file
	## Tool: sambamba markdup
	## Options:
	### -r: remove duplicates instead of just marking them
	### -t: number of threads
	### -p: show progressbar in the terminal
	$sambamba markdup -r -t $nb_threads -p ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.mapped.bam ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.bam

	# ### FLAG: 256 -> remove secondary alignments from the AAV mapping BAM file
	# samtools view -@ $nb_threads -b -F 256 ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.bam > ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.bam 

	# ### FLAG: 2048 -> remove supplementary alignments from the AAV mapping BAM file
	# samtools view -@ $nb_threads -b -F 2048 ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.bam > ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.no_supplementary.bam

	# ### Options: -q 1 -> remove reads with a mapping quality equal to zero (lower than 1)
	# samtools view -b -q 1 ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.no_supplementary.bam > ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam


	# # Sort and index the output filtered BAM file
	# samtools sort -@ $nb_threads ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam -o ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.no_supplementary.MAPQsup0.sorted.bam
	# samtools index ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.no_secondary.no_supplementary.MAPQsup0.sorted.bam

	# Sort and index the output filtered BAM file
	samtools sort -@ $nb_threads ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.bam -o ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.sorted.bam
	samtools index ${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.sorted.bam

done < sample_list_10.txt
