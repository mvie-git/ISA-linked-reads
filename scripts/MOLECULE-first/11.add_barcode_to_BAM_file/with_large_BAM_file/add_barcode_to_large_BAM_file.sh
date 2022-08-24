# Subset a BAM file by a list of read IDs of interest

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM/11.add_barcode_to_BAM_file/TEST"
nb_threads=11

# FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create an output directory
	#mkdir -p ${workdir}/output/${sample}

	# Input files
	input_BAM=${workdir}/input/152.sorted.bam
	index_file=${workdir}/input/subset_152_index${sample}.txt

	# Extract only the readID list of the filtered index file
	echo "Extract read ID to create an index list"
	time awk '{print $1}' $index_file > "${workdir}/input/152_readID_list${sample}.txt"

	# Subset the BAM file with a list of readID of interest
	## Stop grep through the entire BAM file after 2 matches with readID on the index list (R1 and R2)
	echo "Subset BAM file with index file"
	time samtools view -@ $nb_threads $input_BAM | grep -m 90000000 -w -f ${workdir}/input/152_readID_list${sample}.txt > subset_BAM.sam

	# Extract the header
	samtools view -H $input_BAM > header.sam

	# Concatenate header+body and convert it to BAM format
	echo "Concatenate and convert to BAM format"
	time cat header.sam subset_BAM.sam | samtools view -@ $nb_threads -Sb > subset_BAM_${sample}.bam

	# Sort the new BAM file
	echo "Sort new BAM file"
	time samtools sort -@ $nb_threads subset_BAM_${sample}.bam > ${workdir}/input/subset_BAM_${sample}.sorted.bam

	# Index the new BAM file
	echo "Index new BAM file"
	time samtools index ${workdir}/input/subset_BAM_${sample}.sorted.bam

	# Cleaning
	rm subset_BAM.sam header.sam subset_BAM_${sample}.bam


	#----- ADD BARCODE TAG to input BAM file----#

	# Input files
	input_BAM=${workdir}/input/subset_BAM_${sample}.sorted.bam
	index_file=${workdir}/input/subset_152_index${sample}.txt
	output_filename=${workdir}/output/152.BX_tag_${sample}.bam

	# Run pysam script
	echo "Add barcode tag"
	time python3 /media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/pysam/add_barcodes_to_BAM_files_pysam.py $index_file $input_BAM $output_filename


	# Cleaning
	rm "${workdir}/input/152_readID_list${sample}.txt"

done < sample_list.txt
