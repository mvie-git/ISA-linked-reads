# Link: https://rdrr.io/bioc/karyoploteR/man/kpPlotBAMDensity.html
library(karyoploteR)


## TEST ----

# Import data
data.bg <- "/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/covtobed/test_chr1.bg"

# Plot the empty karyoplote
kp <- plotKaryotype(genome="mm10", chromosomes = "chr1")

# Add the data from BAM file
kpArea(kp, data=data.bg, y=data.bg$score, ymin=0, ymax=maxdata.bg$score)

