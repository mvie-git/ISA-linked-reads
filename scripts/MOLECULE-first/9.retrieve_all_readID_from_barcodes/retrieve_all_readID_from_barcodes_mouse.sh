# Search through the whole index to retrieve all the reads assigned to each barcode of interest

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create an output directory
	mkdir -p ${sample}

	# Input files
	index_file=${workdir}/0.prepare_index_file/${sample}/${sample}_index.txt
	subset_index_file=${workdir}/8.retrieve_barcodes_from_readID/${sample}/${sample}_subset_index_file.txt

	# Extract only the barcode list of the filtered index file
	awk '{print $2}' $subset_index_file | sort | uniq > ${workdir}/9.retrieve_all_readID_from_barcodes/${sample}/${sample}_barcodes_list.txt

	# Subset the index file keeping only barcodes assigned to the list of read ID of interest
	LC_ALL=C grep -w -f ${workdir}/9.retrieve_all_readID_from_barcodes/${sample}/${sample}_barcodes_list.txt $index_file > ${workdir}/9.retrieve_all_readID_from_barcodes/${sample}/${sample}_subset_index_file_all_reads.txt

done < sample_list_9.txt
