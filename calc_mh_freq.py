#!/usr/bin/env python3

# load libraries
import sys
import pysam
from math import log, exp
from scipy.special import logsumexp
import numpy as np
from datetime import datetime # for some quick and dirty profiling

# function to estimate expected heterozygosity
# of a microhap from low coverage sequencing data
# within ONE population
# @param reads a list of iterators. Each iterator (pysam) is for reads for one individual
# @param pos a list of positions of SNPs (1-based, on reference) in the window
# @return either the expected heterozygosity as a float or (if not enough information to estimate) the string "NA"
def calcHe(reads, pos):
	eps = 0.01
	# change pos to 0-based
	pos = [x - 1 for x in pos]
	maxI = len(pos) # used several times below
	alleles = [{} for x in range(0,maxI)]
	# pull relevant sites and Qual scores from reads for each individual
	indReads = []
	for ind in reads: # for each individual

		# TODO
		# match up paired reads
		# extract all reads labeled with names
		# then for any duplicates, go to first and fill in any gaps with the second
		
		tempInd = []
		for r in ind: # for each read
			# get sites of SNPs
			# looping through all matches, if too slow, can change to
			# parsing the CIGAR string and jumping
			matches = r.get_aligned_pairs(matches_only = True)
			i = 0
			j = 0
			maxJ = len(matches)
			tempRead = ["N"] * maxI # if read is missing SNPs, it will have "N" in those positions
			while j < maxJ and i < maxI:
				if matches[j][1] < pos[i]:
					j += 1 # go to next base in the read
				elif matches[j][1] == pos[i]:
					# add the base called
					base = r.query_sequence[matches[j][0]]
					tempRead[i] = base
					if base != "N":
						alleles[i][base] = 1
					i += 1
					j += 1
				else:
					i += 1 # go to next snp, stay on same base in read
			tempInd += [tempRead]
		indReads += [tempInd]
	numReads = np.array([j for j in [len(i) for i in indReads] if j > 0])
	numInds = len(numReads) # number of inds with >= 1 read
	numReads = np.sum(numReads) # total number of reads
	skip = False
	# skip if no reads
	if numReads < 1:
		skip = True	
	# detemine alleles to consider
	allSNPAlleles = []
	for d in alleles:
		if len(d) == 0:
			# this position has depth of 0 in this population
			# recommend filtering loci like this out prior to running this algorithm
			skip = True
			break
		allSNPAlleles += [list(d)]
	if skip:
		# no information on entire locus
		return ["NA", numReads, numInds]
	allMH = allSNPAlleles[0]
	for i in range(1, len(allSNPAlleles)):
		allMHnew = []		
		for a in allSNPAlleles[i]:
			allMHnew += [x + a for x in allMH]
		allMH = allMHnew
	# determine genotypes to consider
	genos = []
	for i in range(0, len(allMH)):
		for j in range(i, len(allMH)):
			genos += [[i, j]] # allele1, allele2 INDICIES WITHIN allMH
	
	# calculate log-likelihoods for each genotype for each individual
	indLLH = np.zeros((len(indReads), len(genos))) # indices: individual, genotype; value: log-likelihood
	for i in range(0, len(indReads)): # for each individual
		ind = indReads[i]
		for j in range(0, len(genos)): # for each genotype
			g = genos[j]
			a1 = allMH[g[0]] # character strings
			a2 = allMH[g[1]] # character strings
			tempLLH = 0
			for r in ind: # for each read	
				LH_a1 = 1 # these are likelihoods, NOT log
				LH_a2 = 1
				for k in range(0, len(r)): # for each SNP pos in the read
					# is it missing? then skip
					if r[k] == "N":
						continue
					# does it match allele 1?
					if r[k] == a1[k]:
						LH_a1 *= 1 - eps
					else:
						LH_a1 *= eps
					# does it match allele 2?
					if r[k] == a2[k]:
						LH_a2 *= 1 - eps
					else:
						LH_a2 *= eps
				# switching to log-likelihood here
				tempLLH += log(0.5 * (LH_a1 + LH_a2))
			# save LLH for ind, genos
			indLLH[i,j] = tempLLH

	# MLE (EM algorithm) for allele frequency
	af = np.repeat(1/len(allMH), len(allMH)) # start off w/ equal frequencies
	maxIter = 1000
	numIter = 0
	tolerance = .0001
	L2 = log(2) # b/c python isn't compiled, need to save this value to avaoid repeated computation
	lastLLH = 0
	while numIter < maxIter:
		# Calculate P(Z|reads, af)
		# first calculate log of prior prob of each genotype
		priorZ = np.zeros(len(genos))
		for i in range(0, len(genos)):
			a1 = genos[i][0]
			a2 = genos[i][1]
			if a1 == a2:
				priorZ[i] = 2 * log(af[a1])
			else:
				priorZ[i] = L2 + log(af[a1]) + log(af[a2])
		# now calculate P(Z| reads, af)
		# and calculate likelihood of total model
		totalLLH = 0
		genoProb = np.zeros((indLLH.shape[0], len(genos)))
		for i in range(0, indLLH.shape[0]):
			# calc (relative) probability of each genotype
			# prior * likelihood
			genoProb[i,:] = priorZ + indLLH[i,:]
			temp_lse = logsumexp(genoProb[i,:])
			totalLLH += temp_lse
			# now normalize
			genoProb[i,:] = np.exp(genoProb[i,:] - temp_lse)
		
		# decide whether to break or not
		if (totalLLH - lastLLH) < tolerance and numIter > 0:
			break
		lastLLH = totalLLH
		
		# Maximization of LLH for af (calculate af given P(Z))
		# first calc number of each genotype (fractional)
		genoProb = np.sum(genoProb, 0)
		# now allele frequency
		af = np.zeros(len(allMH))
		for i in range(0, len(genoProb)):
			af[genos[i][0]] += genoProb[i]
			af[genos[i][1]] += genoProb[i]
		# and normalize
		af = af / (2 * indLLH.shape[0])
		
		numIter += 1
		# end while loop for EM
	
	# calculate He from allele frequency and return
	return [1 - np.sum(np.power(af, 2)), numReads, numInds]


def Main():

	print("start ", str(datetime.now()))
	
	# defaults
	popMap = None # -m, maps bam files (individuals) to populations. bamFilePath \t popName
	snpPos = None # -s, snp positions Chr \t Position (1-based)
	wS = 60 # -w, window size
	maxSNPs = 8 # -ms, maximum number of SNPs in a window (if exceeded, window is skipped)
	HeOut = "He_mh.txt" # -o, output file name
	# get command line inputs
	flag = 1
	while flag < len(sys.argv):	#start at 1 b/c 0 is name of script
		if sys.argv[flag] == "-m":
			flag += 1
			popMap = sys.argv[flag]
		elif sys.argv[flag] == "-s":
			flag += 1
			snpPos = sys.argv[flag]
		elif sys.argv[flag] == "-w":
			flag += 1
			wS = int(sys.argv[flag])
		elif sys.argv[flag] == "-ms":
			flag += 1
			maxSNPs = int(sys.argv[flag])
		elif sys.argv[flag] == "-o":
			flag += 1
			HeOut = sys.argv[flag]
		else:
			print("Error: option", sys.argv[flag], "not recognized.")
			return
		flag += 1
	# end while loop for command line arguments
	if popMap is None or snpPos is None:
		print("Error: both -m and -s must be specified")
		return
	
	print("Window size", wS, "bp")

	# load in file and population associations as a dictionary
	# key of pop, value of pysam.AlignmentFiles (one per individual)
	# this will allow pops to be looped through, pulling relevant individuals
	popDict = {}
	with open(popMap, "r") as fileIn:
		line = fileIn.readline()
		while line:
			sep = line.rstrip().split("\t")
			popDict[sep[1]] = popDict.get(sep[1], []) + [pysam.AlignmentFile(sep[0], "rb")]
			line = fileIn.readline()
	# list of pop names so output can have a consistent order
	masterPop = [p for p in popDict]
	
	print("start window eval ", str(datetime.now()))
	
	# run sliding window across target SNP file
	with open(snpPos, "r") as snpLocations, open(HeOut, "w") as HeOutFile:
	
		# write header on output file
		HeOutFile.write("\t".join(["Chr", "Pos"] + masterPop + ["NumReads_" + p for p in masterPop] + ["NumInds_" + p for p in masterPop]) + "\n")
		
		# set up beginning of first window
		cur = snpLocations.readline().rstrip().split("\t")
		curChr = cur[0] # current chromosome window is on
		cur = [int(cur[1])] # current positions in the window
		nextR = snpLocations.readline() # next available position
		
		# loop through all possible windows
		numSnpsPerWindow = {}
		numWindows = 0
		while nextR:
			nextR = nextR.rstrip().split("\t")
			nextR[1] = int(nextR[1])
			# determine if the next variant is within the window
			if nextR[0] == curChr and (nextR[1] - cur[0]) < wS:
				cur += [nextR[1]] # add to window
			else:
				# evaluate current window
				numWindows += 1
				numSnpsPerWindow[len(cur)] = numSnpsPerWindow.get(len(cur), 0) + 1
				
				# skip if too many SNPs
				if len(cur) <= maxSNPs:
					# get reads for each individual
					# remember you don't want to have multiple fetch iterators being used
					# on the same file simultaneously unless you open multiple file handles
					reads = {} # this is dict, key is pop, value is list of iterators (one iterator per ind)
					for pop in masterPop:
						reads[pop] = [f.fetch(contig=curChr, start=(cur[0] - 1), stop=cur[-1]) for f in popDict[pop]]
						
					# now calculate expHet within each pop
					# this could be done in the previous loop, but seperating it out
					# to allow easy transition to calculating each pop in parallel if
					# later desired
					# create list of [str(He), str(He), ...] in order of masterPop and number of reads and inds
					tempRes = [calcHe(reads[pop], cur) for pop in masterPop]
					He_results = [str(i[0]) for i in tempRes]
					numReads = [str(i[1]) for i in tempRes]
					numInds = [str(i[2]) for i in tempRes]
					# and write to output
					HeOutFile.write("\t".join([curChr, ",".join([str(x) for x in cur])] + He_results + numReads + numInds) + "\n")
					
					# TODO
					# optional output of allele frequencies
				
				# advance to next window
				if curChr == nextR[0]:
					# have to include the next snp, otherwise
					# new window is a subset of the previous window
					cur = [x for x in cur if (nextR[1] - x) < wS]
					cur += [nextR[1]]
				else:
					curChr = nextR[0]
					cur = [nextR[1]]
			nextR = snpLocations.readline()
		# print some summary information
		print("number of windows: ", numWindows)
		print("Distribution of SNPs per window")
		print("SNPs ", "NumberOfWindows")
		for k in sorted(list(numSnpsPerWindow)):
			print(k, numSnpsPerWindow[k])
		
		print("end ", str(datetime.now()))

if __name__ == '__main__':
	Main()
