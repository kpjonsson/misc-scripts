#!/opt/common/CentOS_6-dev/bin/current/python
descr = "Takes a MAF file as input, returns count of indels with microhomology in 5' and 3' flanking sequences per sample.\nRequires ntdpal tool from primer3."
import argparse, string, pysam, subprocess
# import argparse, subprocess, os, re, glob

parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-m', '--mafFile', help = 'Input MAF', required = True)
parser.add_argument('-l', '--flankLength', help = 'Length of flanking regions to search', type = int, default = 50, required = False)
parser.add_argument('-o', '--outFile', help = 'Name output', required = True)
args = parser.parse_args()

### Microhomology // longest commom substring problem
### http://stackoverflow.com/questions/2892931/longest-common-substring-from-more-than-two-strings-python
def longestSubstring(data):
	substrs = lambda x: {x[i:i+j] for i in range(len(x)) for j in range(len(x) - i + 1)}
	s = substrs(data[0])
	for val in data[1:]:
		s.intersection_update(substrs(val))
	return max(s, key=len)

### Reference 
refFasta = pysam.FastaFile('/ifs/depot/assemblies/H.sapiens/GRCh37/gr37.fasta')

### Read file and find index for necessary columns
maf = open(args.mafFile, 'r') 
header = maf.readline()
if header.startswith('#'): # Check if comments in header
	header = maf.readline()

columnLabels = header.strip().split('\t')
columnIndex = {}
for i,x in enumerate(columnLabels):
    columnIndex[x] = i

sampleIndex = columnIndex['Tumor_Sample_Barcode']
typeIndex = columnIndex['Variant_Type']
chrIndex = columnIndex['Chromosome']
startIndex = columnIndex['Start_Position']
endIndex = columnIndex['End_Position']
refIndex = columnIndex['Reference_Allele']
altIndex = columnIndex['Tumor_Seq_Allele2']

### Parse lines in MAF, check for microhomolgy if indel
with open(args.outFile, 'w') as output:
	# Write header
	output.write('Sample\tChrom\tStart\tEnd\tVariant_Type\tRef_Allele\tAlt_Allele\tFive_Prime_Flank\tThree_Prime_Flank\tMicrohomology\tHomology_Length\n')
	
	for line in maf:
		sline = map(string.strip, line.split('\t'))

		# Check if indel
		variantType = sline[typeIndex]
		if variantType not in ['INS', 'DEL']:
			continue

		# If indel
		sample = sline[sampleIndex]
		chrom = sline[chrIndex]
		start = int(sline[startIndex])
		end = int(sline[endIndex])
		refAllele = sline[refIndex]
		altAllele = sline[altIndex]

		# Run ntdpal
		if variantType == 'DEL':
			fiveFlank = refFasta.fetch(reference = chrom, start = (start-1)-args.flankLength, end = (start-1))
		else: 
			fiveFlank = refFasta.fetch(reference = chrom, start = start-args.flankLength, end = start)
		threeFlank = refFasta.fetch(reference = chrom, start = end, end = end+args.flankLength)
		homology = longestSubstring([fiveFlank, threeFlank])
		homLength = len(homology)

		# Parse and print output, but only if homology > 1 bp
		#if homLength > 1:
		outline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (sample, chrom, start, end, variantType, refAllele, altAllele, fiveFlank, threeFlank, homology, homLength)
		output.write(outline)
