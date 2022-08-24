# Select from a BAM file a list of read IDs of interest

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create output directory
	mkdir -p ${sample}

	# Input files
	#input_BAM=${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.no_supplementary.MAPQsup0.sorted.bam
	input_BAM=${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.sorted.bam

	# Extract read ID column only
	samtools view ${input_BAM} | awk '{print $1}' | sort | uniq > ${workdir}/7.select_readID_of_interest/${sample}/${sample}_readID_list.txt

done < sample_list_7.txt
