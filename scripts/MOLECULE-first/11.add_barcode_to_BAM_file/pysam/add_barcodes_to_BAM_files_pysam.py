
import pysam
import sys

barcodes = {}

# Retrieve barcodes for each read ID
#with open('output/output_151_index_AAV_plus_reads.txt') as barcodes_file:
with open(sys.argv[1]) as barcodes_file:
  for line in barcodes_file:
    read_id, barcode = line.rstrip().split()
    barcodes[read_id] = barcode

# Set the tag names (here BC for barcode tag)
barcode_tag = 'BX'

# Input BAM file
inbam = pysam.Samfile(sys.argv[2], "rb") # input BAM file
outbam = pysam.Samfile(sys.argv[3], "wb", template=inbam) # output BAM filename

# inbam = pysam.Samfile("output/filtered_AAV_plus_reads.bam", "rb")
# outbam = pysam.Samfile("output/filtered_AAV_plus_reads.BX_tag.bam", "wb", template=inbam)

for read in inbam.fetch(until_eof=True):
    for read in inbam:

#      read.tags += [('BX', barcodes[read.qname])]

      try:
        read.tags += [('BX', barcodes[read.qname])]
      except KeyError:
      # No read ID - either skip read
       continue

      outbam.write(read)

outbam.close()
