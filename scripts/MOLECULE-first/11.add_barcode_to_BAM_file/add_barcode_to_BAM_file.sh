# Add the barcode information as a tag in the new output BAM file with a tool named pysam which allow to parse a BAM file

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"
nb_threads=12

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create an output directory
	mkdir -p ${sample}

	# Input files
	input_BAM=${workdir}/10.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.no_duplicates.sorted.bam
	index_file=${workdir}/8.retrieve_barcodes_from_readID/${sample}/${sample}_subset_index_file.txt
	output_filename=${workdir}/11.add_barcode_to_BAM_file/${sample}/${sample}_subset_BAM_filtered.BX_tag.bam

	# Run pysam script
	python3 /media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/pysam/add_barcodes_to_BAM_files_pysam.py $index_file $input_BAM $output_filename

done < sample_list_11.txt
