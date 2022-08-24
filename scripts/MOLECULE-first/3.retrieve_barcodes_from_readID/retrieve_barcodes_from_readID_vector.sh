# Search through the whole index to retrieve the barcodes assigned to each read ID of interest

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create an output directory
	mkdir -p ${sample}

	# Input files
	index_file=${workdir}/0.prepare_index_file/${sample}/${sample}_index.txt
	readID_list=${workdir}/2.select_readID_of_interest/${sample}/${sample}_readID_list.txt

	# Subset the index file keeping only barcodes assigned to the list of read ID of interest
	#LC_ALL=C grep -w -f $readID_list $index_file > ${workdir}/3.retrieve_barcodes_from_readID/${sample}/${sample}_subset_index_file.txt
	grep -w -f $readID_list $index_file > ${workdir}/3.retrieve_barcodes_from_readID/${sample}/${sample}_subset_index_file.txt

done < sample_list_3.txt
