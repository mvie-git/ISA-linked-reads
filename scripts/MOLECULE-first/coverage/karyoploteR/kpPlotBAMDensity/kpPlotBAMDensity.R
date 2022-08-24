# Link: https://rdrr.io/bioc/karyoploteR/man/kpPlotBAMDensity.html
library(karyoploteR)

## EXAMPLE ----

#install.packages("pasillaBamSubset") # Download package 0.34.0 at bioconductor website
library(pasillaBamSubset) #A package with 2 example bam files

un1.bam.file <- untreated1_chr4() # get the name of the first bam
un3.bam.file <- untreated3_chr4() #and the name of the second

window.size <- 1e4 #compute the density with 10kb windows

kp <- plotKaryotype(genome="dm6", chromosomes="chr4") #The pasilla data comes from drosophila
kp <- kpAddBaseNumbers(kp, tick.dist = 1e5)
kp <- kpPlotBAMDensity(kp, data = un1.bam.file, window.size = window.size, r0=0.5, r1=1, ymax=50000, col="darkorange")
kp <- kpPlotBAMDensity(kp, data = un3.bam.file, window.size = window.size, r0=0.5, r1=0, ymax=50000, col="darkorchid") #using r0>r1 we can flip the plot
kpAxis(kp, ymin=0, ymax=50000, r0=0.5, r1=1, labels = c("0", "25K", "50K"))
kpAxis(kp, ymin=0, ymax=50000, r0=0.5, r1=0, labels = c("0", "25K", "50K"))

kpText(kp, chr = "chr4", x=7e5, y=0.85, labels = paste0("Untreated 1 (reads per ", window.size, " bases)"))
kpText(kp, chr = "chr4", x=7e5, y=0.15, labels = paste0("Untreated 3 (reads per ", window.size, " bases)"))


#Or normalizing by the number of mapped reads
kp <- plotKaryotype(genome="dm6", chromosomes="chr4") #The pasilla data comes from drosophila
kp <- kpAddBaseNumbers(kp, tick.dist = 1e5)
kp <- kpPlotBAMDensity(kp, data = un1.bam.file, window.size = window.size, normalize=TRUE, r0=0.5, r1=1, ymax=0.2, col="darkorange")
kp <- kpPlotBAMDensity(kp, data = un3.bam.file, window.size = window.size, normalize=TRUE, r0=0.5, r1=0, ymax=0.2, col="darkorchid") #using r0>r1 we can flip the plot


## TEST ----

# Import data
data.bam <- "/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/output/152_mapping_mouseANDvector.sorted.bam"

# Plot the empty karyoplote
kp <- plotKaryotype(genome="mm10", chromosomes = "chr1")

# Add the data from BAM file
kpPlotBAMDensity(kp, data=data.bam)

