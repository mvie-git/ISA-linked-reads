# Sort and index the new BAM file

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"
nb_threads=12

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Input files
	input_BAM=${workdir}/output/pysam/${sample}_subset_BAM.BX_tag.bam
	output_BAM=${workdir}/output/pysam/${sample}_subset_BAM.BX_tag.sorted.bam

	# Sort and index the new BAM tagged file
	samtools sort -@ $nb_threads $input_BAM -o $output_BAM
	samtools index $output_BAM

	# Cleaning...
#	rm $output_filename

done < sample_list_7.txt
