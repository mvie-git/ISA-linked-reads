# Select from a BAM file a list of read IDs of interest

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create output directory
	mkdir -p ${sample}

	# Input files
	#input_BAM=$workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam
	input_BAM=$workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.sorted.bam

	# Extract read ID column only
	samtools view ${input_BAM} | awk '{print $1}' | sort | uniq > ${workdir}/2.select_readID_of_interest/${sample}/${sample}_readID_list.txt

done < sample_list_2.txt
