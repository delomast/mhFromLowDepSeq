'''
Created on May 4, 2022

Note that this script can produce some duplicates in the output
Should be filtered afterwards to make sure there aren't any

Filtering lamprey vcf to only include SNPs covered
by read 2 (far enough away from cutsites)
Example usage:
	 python3 filterLamp.py -b cov.stats -v newVCF.vcf.gz -f /data/usda/genome_indices/lamprey/ETRm_v1.fasta
@author: Thomas Delomas
'''

import pysam
from pysam import VariantFile
import sys
import re

def Main():
	# get command line arguments
	# need: bed file -b, vcf file -v, fasta file -f (.fai assumed also present), 
	# cutsite -c, offset -o (number of bases in cutsite to the left of cut location)
	# read length -l
	flag = 1
	bedFile = None
	vcfFile = None
	fastaFile = None
	cutsite = "CCTGCAGG" # default is sbfI
	cutOff = 2 # default is sbfI
	readLen = 101
	while flag < len(sys.argv):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == "-b":
			flag += 1
			bedFile = sys.argv[flag]
		elif sys.argv[flag] == "-v":
			flag += 1
			vcfFile = sys.argv[flag]
		elif sys.argv[flag] == "-f":
			flag += 1
			fastaFile = sys.argv[flag]
		elif sys.argv[flag] == "-c":
			flag += 1
			cutsite = sys.argv[flag]
		elif sys.argv[flag] == "-o":
			flag += 1
			cutOff = int(sys.argv[flag])
		elif sys.argv[flag] == "-l":
			flag += 1
			readLen = int(sys.argv[flag])
		else:
			print("Error: option", sys.argv[flag], "not recognized.")
			return
		flag += 1
	# end while loop for command line arguments
	
	with open(bedFile, "r") as bed, open("filteredSNPs.txt", "w") as snpsOut:
		cutRE = re.compile(cutsite, flags = re.I)
		# open vcf
		vcf = VariantFile(vcfFile)
		# open fasta
		refGenome = pysam.FastaFile(fastaFile)
		# loop through all intervals in bed file
		bLine = bed.readline()
		while bLine:
			sep = bLine.rstrip().split("\t")
			sep[1] = int(sep[1])
			sep[2] = int(sep[2])
			# make sure long enough to contain a cutsite and a read 2
			if sep[2] - sep[1] > readLen:
				# find SNPs
				snps = [x.pos for x in vcf.fetch(contig = sep[0], start = sep[1], stop = sep[2])] # 0-based half open, same as bed
				if len(snps) > 0:
					# find cutsite, adding a window of 500 in case cutsite was excluded from the bed interval
					refStart = max(sep[1] - 500, 0)
					refEnd = min(sep[2] + 500, refGenome.get_reference_length(sep[0]))
					m = [x for x in cutRE.finditer(refGenome.fetch(reference = sep[0],  start = refStart, end = refEnd))]
					# only use SNPs in intervals of one cutsite
					if len(m) == 1:
						posCutLoc = m[0].start() + refStart + cutOff - 1 # this is last base BEFORE cut on positive strand
						negCutLoc = m[0].start() + refStart + len(cutsite) - cutOff + 1 # this is last base BEFORE cut on negative strand
						# trim SNPs to avoid read 1 (snp positions in snps are 1-based)
						snps = [x for x in snps if (((x - 1) < (negCutLoc - readLen)) or ((x - 1) > (posCutLoc + readLen)))]
						# write out for microhap estimation
						for s in snps:
							snpsOut.write("\t".join([sep[0], str(s)]) + "\n")
			
			bLine = bed.readline()
		
		# close fasta and vcf
		refGenome.close()
		vcf.close()


if __name__ == '__main__':
	Main()
