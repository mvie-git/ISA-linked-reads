# Subset a BAM file by a list of read IDs of interest

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"
nb_threads=12

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create an output directory
	mkdir -p ${sample}

	# Input files
	input_BAM=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/output/${sample}/${sample}_mapping_mouseANDvector.sorted.bam
	index_file=${workdir}/4.retrieve_all_readID_from_barcodes/${sample}/${sample}_subset_index_file_all_reads.txt

	# Extract only the readID list of the filtered index file
	awk '{print $1}' $index_file | sort | uniq > ${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_readID_list.txt

	# Subset the BAM file with a list of readID of interest
	samtools view -@ $nb_threads $input_BAM | grep -w -f ${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_readID_list.txt > ${sample}_subset_BAM.sam

	# Extract the header
	samtools view -H $input_BAM > header.sam

	# Concatenate header+body and convert it to BAM format
	cat header.sam ${sample}_subset_BAM.sam | samtools view -@ $nb_threads -Sb > ${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.bam

	# Sort the new BAM file
	samtools sort ${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.bam > ${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.bam

	# Index the new BAM file
	samtools index ${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.bam

	# Cleaning
	rm header.sam ${sample}_subset_BAM.sam ${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.bam


done < sample_list_5.txt
