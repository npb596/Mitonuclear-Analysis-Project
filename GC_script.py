#!/bin/python

# Import Bio modules necessary to find gene position in a fasta genome file and sys module to take arguments

from Bio import SeqIO
from Bio.SeqUtils import GC
import sys

# Define arguments for (1) reference genome in fasta format, (2) chromosome, (3) start position, (4) end position

ref_path = sys.argv[1]
contig = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

# Define reference genome as well as specific locus position, then find and print GC content of that locus

ref = SeqIO.index(ref_path, "fasta")
ref_seq = ref[contig].seq[start:end]
print(GC(ref_seq))
