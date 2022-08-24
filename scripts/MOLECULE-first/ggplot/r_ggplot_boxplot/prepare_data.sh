workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping"

while read sample; do

	echo "*** Start $sample ***"
	
	input="${workdir}/get_molecules/${sample}/${sample}.molecules.bed"

	# Boxplot on number of molecules per barcode
	#sed 1d $input | awk '{print $4}' | sort | uniq -c | awk '{print $1}' | awk '{print '${sample}' "\t" $0}' > $workdir/ggplot/r_ggplot_boxplot/input/${sample}.csv

	# Boxplot on the length of chimeric molecules
	sed 1d $input | awk '{print $6}' | awk '{print '${sample}' "\t" $0}' > $workdir/ggplot/r_ggplot_boxplot/input/${sample}.csv


done < sample_list.txt


# Gather all data
cat $workdir/ggplot/r_ggplot_boxplot/input/* > $workdir/ggplot/r_ggplot_boxplot/input/all.csv
