workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping"

while read sample; do

	echo "*** Start $sample ***"
	
	input="${workdir}/get_molecules/${sample}/${sample}.molecules.bed"

	# Histogram on number of reads per molecule data
	#sed 1d $input | awk '{print $5}' | awk '{print '${sample}' "\t" $0}' > $workdir/ggplot/r_ggplot_histogram/input/${sample}.csv

	# # Histogram on length of molecules
	sed 1d $input | awk '{print $6}' | awk '{print '${sample}' "\t" $0}' > $workdir/ggplot/r_ggplot_histogram/input/${sample}.csv

done < sample_list.txt


# Gather all data
cat $workdir/ggplot/r_ggplot_histogram/input/* > $workdir/ggplot/r_ggplot_histogram/input/all.csv
