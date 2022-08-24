while read sample; do

	echo "*** Start $sample ***"

	# Set up parameters
	output_virus_clip="/media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip/GSDIa_samples_after_AAV_mapping/result/${sample}/remove_FP/virus_clip.filter_FP.out"

	# If the left-element of the soft-clipped read is a mouse, extract the integration loci and calculate the length of the clipped read to guess the end of the mouse integration loci (part 1)
	echo "*** Starting the extraction of the integration loci if the left-element of the soft-clipped read is a mouse ***"
	awk '$1=="mouse" {print $2 "\t" $3 "\t" $4}' $output_virus_clip | awk '{print $0"\t"length($3)}' | awk '$5=$2+$4' | awk '{print $1 "\t" $2 "\t" $5}' > mouse_integration_loci_left.bed
	echo "*** Extract the left-element of soft-clipped reads done ***"

	# If the right-element of the soft-clipped read is a mouse, extract the integration loci and calculate the length of the clipped read to guess the end of the mouse integration loci (part 2)
	echo "*** Starting the extraction of the integration loci if the right-element of the soft-clipped read is a mouse ***"
	awk '$5=="mouse" {print $6 "\t" $7 "\t" $8}' $output_virus_clip | awk '{print $0"\t"length($3)}' | awk '$5=$2+$4' | awk '{print $1 "\t" $2 "\t" $5}' > mouse_integration_loci_right.bed
	echo "*** Extract the right-element of soft-clipped reads done ***"

	# Concatenate the left and right clipped-reads
	echo "*** Starting the concatenation of the left and right clipped-reads ***"
	cat mouse_integration_loci_left.bed mouse_integration_loci_right.bed > /media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip/GSDIa_samples_after_AAV_mapping/result/${sample}/remove_FP/mouse_integration_loci_FP_filtered.bed
	echo "*** Concatenate the soft-clipped done ***"

	# Remove the 'chr' 
	echo "*** Remove the 'chr' ***"
	sed -i 's/chr//g' /media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip/GSDIa_samples_after_AAV_mapping/result/${sample}/remove_FP/mouse_integration_loci_FP_filtered.bed
	echo "*** Removing the 'chr' is done ***"

	# Add a header: Chr, Start, End
	echo "*** Add an header for the new BED file ***"
	sed -i '1 i\Chr\tStart\tEnd' /media/mallaury/MALLAURY/PROJECT/GSDIa_ISA/search_chimeric_reads/virus_clip/GSDIa_samples_after_AAV_mapping/result/${sample}/remove_FP/mouse_integration_loci_FP_filtered.bed
	echo "*** BED file created ***"

	# Clean intermediate files
	echo "*** Cleaning intermediate files ***"
	rm mouse_integration_loci_left.bed
	rm mouse_integration_loci_right.bed


	echo "*** Script is done. ***"

done < sample_list_IS_mouse_genome.txt