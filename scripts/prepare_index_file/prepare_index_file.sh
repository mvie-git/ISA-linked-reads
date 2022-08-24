# Extract barcode for each read ID present in each sample (format: 1.Read ID; 2.Barcode)

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"

## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create the output directory
	mkdir -p $workdir/0.prepare_index_file/${sample}

	# Input files
	index_file=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/tell_seq/output/${sample}/demultiplex/Full/demultiplex_AAV_I1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz

	# Conversion FASTQ to FASTA
	zcat $index_file | sed -n '1~4s/^@/>/p;2~4p' > ${workdir}/${sample}_index.fasta

	# Put the read ID and the barcode sequence on the same line | remove '>' character & empty lines
	awk 'BEGIN {RS=">";OFS="\t"} {print $1,$2}' ${workdir}/${sample}_index.fasta | sed 's/>//' | awk 'NF' > ${workdir}/0.prepare_index_file/${sample}/${sample}_index.txt

	# Clean
	rm ${workdir}/${sample}_index.fasta

done < sample_list_0.txt
