# Here, we need to convert BAM file to FASTQ files otherwise go to the last command to launch virus-clip software!

while read sample; do

	echo "*** Start $sample ***"

	# Set up parameters
	BAM_file="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.bam"
	sorted_BAM_file="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.bam"
	FASTQ_R1_file="${sample}_1.fastq"
	FASTQ_R2_file="${sample}_2.fastq"
	input_data_r1="${sample}_1.fastq.gz"
	input_data_r2="${sample}_2.fastq.gz"
	single_paired=2


	# But before the BAM file need to be sorted by name (sort paired read alignment):
	echo "*** Starting the sort of the BAM file by read name ***"
	samtools sort -n $BAM_file -o $sorted_BAM_file
	echo "*** BAM file sorted by name done***"

	# Save fastq reads in separate R1 and R2 files
	echo "*** Starting to convert the BAM file into two FASTQ files ***"
	bedtools bamtofastq -i $sorted_BAM_file \
	-fq $FASTQ_R1_file \
	-fq2 $FASTQ_R2_file
	echo "*** BAM file converted in FASTQ format done***"

	# Compress fastq files
	echo "*** Starting the compression of FASTQ files ***"
	pigz $FASTQ_R1_file
	pigz $FASTQ_R2_file
	echo "*** Compress FASTQ files done ***"

	# Copy FASTQ files into the input folder of virus-clip software
	echo "*** Starting the copy of the compressed FASTQ files to the input folder for Virus-Clip analysis ***"
	mv "${FASTQ_R1_file}.gz" /media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip/input
	mv "${FASTQ_R2_file}.gz" /media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip/input
	echo "*** Copy compressed FASTQ files to the input folder for Virus-Clip analysis done ***"

	# Execution of virus-clip script
	echo "*** Starting the execution of the virus-clip software... ***"
	bash /media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/tools/Virus-Clip/Virus-Clip_201201/virus_clip.sh \
	/media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip/input \
	/media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip \
	$sample \
	$single_paired

	echo "*** Virus-Clip is done. ***"

done < sample_list_virus_clip.txt