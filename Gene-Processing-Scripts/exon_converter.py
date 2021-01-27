#!/usr/bin/env python3

# Import system module so that command-line arguments can be read
# After that define the first argument and second arguments as input and output files respectively
# The input file should be an ensGene file with standard HUGO gene names added at the end

import sys
InFile = open(sys.argv[1], 'r')
OutFile = open(sys.argv[2], 'w')

# Define LineNumber as 0 for good practice according to Haddock and Dunn

LineNumber = 0

# Loop through every line in the input file, strip the line endings off (good practice according to Haddock and Dunn)
# Split lines based on tab-delimitation common to ensGene files
# Split the exon start and ends based on comma separation (again, common to ensGene files) and remove the final value, which is just an empty space

for Line in InFile:
	Line = Line.strip('\n')
	ElementList = Line.split('\t')
	
	Starts=ElementList[9].split(',')
	Starts=Starts[:-1]
	
	Ends=ElementList[10].split(',')
	Ends=Ends[:-1]
	
	LineNumber = LineNumber + 1

# Number of exons is defined in the ensGene file for each coding transcript
# So for each transcript give a bed formatted line with only one exon start and end per line as well as transcript and HUGO gene name at end

	for x in range(0,int(ElementList[8])):
		OutFile.write(ElementList[2]+"\t"+Starts[x]+"\t"+Ends[x]+"\t"+ElementList[1]+"\t"+ElementList[16]+"\n")

# Close the input and output files

InFile.close()
OutFile.close()
