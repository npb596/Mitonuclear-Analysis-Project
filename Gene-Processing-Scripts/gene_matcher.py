#!/usr/bin/env python3

# Import modules sys to take arguments, os to remove temporary files, and random to select a random gene

import sys
import os
import random

# Read first argument with bed file of query genes to be matched. Read second argument with bed file of subject genes to pull when they match query genes.

Query = open(sys.argv[1], 'r')
Subject = open(sys.argv[2], 'r')

# Loop through all lines of query bed file.

for QLine in Query:

# Strip each line of query file and tab-separate columns. Then open temporary file for subject genes that are on the same chromosome as query genes.

	QLine = QLine.strip('\n')
	ElementList = QLine.split('\t')
	ChrFile = open("Chr.bed", 'w')

# Loop through subject file, tab-separate columns, and print gene info to temp file onl if the genes match query chromosomes. Then close out temp file.

	for SLine in Subject: 
		SLine = SLine.strip('\n')
		ElementSub = SLine.split('\t')
		if ElementList[0] == ElementSub[0]:
			ChrFile.write(ElementSub[0]+"\t"+ElementSub[1]+"\t"+ElementSub[2]+"\t"+ElementSub[3]+"\t"+ElementSub[4]+"\t"+ElementSub[5]+"\n")
	Subject.seek(0)
	ChrFile.close()

# Open previous temp file to read and open up new temp file to write out new genes that match query recombination rates.

	ChrFile = open("Chr.bed", 'r')
	RecomFile = open("Recom.bed", 'w')

# Loop through chromosome temp file, convert recombination rates to floats, then create variables for upper and lower recombination bounds based on query genes.

	for ChrLine in ChrFile:
		ChrLine = ChrLine.strip('\n')
		ElementChr = ChrLine.split('\t')
		QueryRecom = float(ElementList[4])
		SubjectRecom = float(ElementChr[4])
		upper_query=QueryRecom+(QueryRecom*0.3)
		lower_query=QueryRecom-(QueryRecom*0.3)
	
# Only print subject gene info if recombination rate falls within the upper and lower bounds defined above.

		if SubjectRecom < upper_query and SubjectRecom > lower_query:
			RecomFile.write(ElementChr[0]+"\t"+ElementChr[1]+"\t"+ElementChr[2]+"\t"+ElementChr[3]+"\t"+ElementChr[4]+"\t"+ElementChr[5]+"\n")
	RecomFile.close()
	
# Largely similar to the above blocks of code but now matching based on GC content and then matching based on gene length.

	RecomFile = open("Recom.bed", 'r')
	GCFile = open("GC.bed", 'w')
	for RecomLine in RecomFile:
                RecomLine = RecomLine.strip('\n')
                ElementRecom = RecomLine.split('\t')
                QueryGC = float(ElementList[5])
                SubjectGC = float(ElementRecom[5])
                upper_query=QueryGC+(QueryGC*0.3)
                lower_query=QueryGC-(QueryGC*0.3)
                if SubjectGC < upper_query and SubjectGC > lower_query:
                        GCFile.write(ElementRecom[0]+"\t"+ElementRecom[1]+"\t"+ElementRecom[2]+"\t"+ElementRecom[3]+"\t"+ElementRecom[4]+"\t"+ElementRecom[5]+"\n")
	GCFile.close()
	GCFile = open("GC.bed", 'r')
	LengthFile = open("Length.bed", 'w')
	for GCLine in GCFile:
		GCLine = GCLine.strip('\n')
		ElementGC = GCLine.split('\t')
		QueryLength = float(ElementList[2])-float(ElementList[1])
		SubjectLength = float(ElementGC[2])-float(ElementGC[1])
		upper_query=QueryLength+(QueryLength*0.3)
		lower_query=QueryLength-(QueryLength*0.3)
		if SubjectLength < upper_query and SubjectLength > lower_query:
			LengthFile.write(ElementGC[0]+"\t"+ElementGC[1]+"\t"+ElementGC[2]+"\t"+ElementGC[3]+"\t"+ElementGC[4]+"\t"+ElementGC[5]+"\t"+ElementList[3]+"\n")
	LengthFile.close()	

# Out of all matching subject genes pick a random  one for each query gene and output to the file given as third argument. If there are no matching genes, print an error messaage and move on.

	try:
		FinalFile = open(sys.argv[3], 'a')
		FinalFile.write(random.choice(list(open("Length.bed"))))
	except IndexError:
		print("No genes matched criteria, cannot select random gene")

# Close all files used as arguments and remove temporary files

Query.close()
Subject.close()
FinalFile.close()

os.remove("Chr.bed")
os.remove("Recom.bed")
os.remove("GC.bed")
os.remove("Length.bed")
