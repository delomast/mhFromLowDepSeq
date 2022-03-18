#!/usr/bin/env python3

# load libraries
import sys
import pysam
from math import log, exp
import numpy as np

# log of the sum of the exponents
# @param a a numpy array, numeric
def lse(a):
	c = a.max()
	return c + np.log(np.sum(np.exp(a - c)))

# convert Q on the phred scale to P - prob the base was sequenced wrong
# @param Q an int, note that this has already been converted from ASCII
#    and had the 33 subtracted from it
def phredQtoP(Q):
	return 10**(-Q/10)
	
# calculate probability the basecall is wrong given the probability
# that the base on the template was wrong and the Q (from sequencer phred score)
# assumes that if an error is made in either spot, the probabilty of being any other base is equal
# @param eps probability that the base on the template in the sequencer was wrong
# @param e2 prob the base was sequenced wrong, the output of phredQtoP(Q)
def probSubErr(eps, e2):
	# template wrong, sequencer correct + 
	#   template wrong, sequencer wrong and didn't randomly call correct +
	#   template right, sequencer wrong
	return (eps * (1 - e2)) + (0.6666666667 * eps * e2) + ((1 - eps) * e2)


# function to run EM algorithm for estimating either 
# expected heterozygosity or allele frequencies
# of a microhap from low coverage sequencing data
# within ONE population
# @param reads a list of iterators. Each iterator (pysam) is for reads for one individual
# @param pos a list of positions of SNPs (0-based, on reference) in the window
# @param alleleFreq boolean to indicate whether to output allele frequencies (True) or just expectedHet (False)
# @param minAF minimum allele frequency to keep an allele in the analysis
# @return either the expected heterozygosity as a float or (if not enough information to estimate) the string "NA"
def calcEM(reads, pos, alleleFreq, minAF, epsTemplate):
	# change pos to 0-based
	pos = [x - 1 for x in pos]
	maxI = len(pos) # used several times below
	alleles = [{} for x in range(0,maxI)]
	# pull relevant sites and Qual scores from reads for each individual
	indReads = []
	indQuals = []
	for ind in reads: # for each individual

		readNames = []
		tempInd = []
		tempQual = []
		for r in ind: # for each read
			# get base calls for SNPs
			j = r.reference_start # position in the reference, 0 based
			k = 0 # position in the query
			c = 0 # position in the CIGAR string
			# determine if it's the pair of an earlier read
			pair = False
			if r.query_name in readNames:
				pair = True
				matePos = readNames.index(r.query_name)
				tempRead = tempInd[matePos]
				tempReadQual = tempQual[matePos]
			else:
				tempRead = ["N"] * maxI # if read is missing SNPs, it will have "N" in those positions
				tempReadQual = [1] * maxI # if read is missing SNPs, it will have "N" in those positions
			for i in range(0, maxI): # for each SNP
				# check if passed SNP
				if j > pos[i]:
					continue # go to the next SNP
				# if paired, only overwrite "N"s
				if pair and tempRead[i] != "N":
					continue
				for cPos in range(c, len(r.cigartuples)):
					if r.cigartuples[cPos][0] in (0,7,8):
						# match
						# if target SNP is in the range
						if pos[i] < (j + r.cigartuples[cPos][1]) and pos[i] >= j: # >= j is defensive, should never be False b/c would have skipped
							# pull the target base from the query
							base = r.query_sequence[k + pos[i] - j]
							# pull the target quality from the query
							qual = phredQtoP(r.query_qualities[k + pos[i] - j])
							if pair and qual >= tempReadQual[i]:
								# if paired reads overlap at this base, then
								# skip if qual is lower than existing
								break
							else:
								tempRead[i] = base
								tempReadQual[i] = qual
							if base != "N":
								alleles[i][base] = 1
							# we don't advance anything b/c the current positions are
							# where we want to start to look for the next SNP
							# just move to the next SNP
							break
						else:
							# consume both
							k += r.cigartuples[cPos][1]
							j += r.cigartuples[cPos][1]
							c += 1
					elif r.cigartuples[cPos][0] in (1,4):
						# insertion or soft clip, consume query not ref
						k += r.cigartuples[cPos][1]
						c += 1
					elif r.cigartuples[cPos][0] in (2,3):
						# deletion or skip, consume ref not query
						j += r.cigartuples[cPos][1]
						c += 1
					elif r.cigartuples[cPos][0] in (5,6):
						# hard clip or padding, do nothing
						pass
					else:
						raise RuntimeError("unrecognized CIGAR operation")
					# check if passed SNP
					if j > pos[i]:
						break # go to the next SNP

			if pair:
				tempInd[matePos] = tempRead
				tempQual[matePos] = tempReadQual
			else:
				tempInd += [tempRead]
				tempQual += [tempReadQual]
				readNames += [r.query_name]
		# if one or more read, save and use in estimates
		if len(tempInd) > 0:
			indReads += [tempInd]
			indQuals += [tempQual]
	numReads = np.array([len(i) for i in indReads])
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
		lineOut = ["NA", str(numReads), str(numInds)]
		if alleleFreq:
			lineOut += ["NA"]
		return lineOut
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
	genos = np.array(genos)
	
	# calculate log-likelihoods for each genotype for each individual
	indLLH = np.zeros((len(indReads), genos.shape[0])) # indices: individual, genotype; value: log-likelihood
	for i in range(0, len(indReads)): # for each individual
		for j in range(0, genos.shape[0]): # for each genotype
			a1 = allMH[genos[j,0]] # character strings
			a2 = allMH[genos[j,1]] # character strings
			tempLLH = 0
			for m in range(0, len(indReads[i])): # for each read
				r = indReads[i][m]
				q = indQuals[i][m]
				LH_a1 = 1 # these are likelihoods, NOT log
				LH_a2 = 1
				for k in range(0, len(r)): # for each SNP pos in the read
					# is it missing? then skip
					if r[k] == "N":
						continue
					# calc prob of error
					eps = probSubErr(epsTemplate, q[k])
					# does it match allele 1?
					if r[k] == a1[k]:
						LH_a1 *= 1 - eps
					else:
						LH_a1 *= eps / 3
					# does it match allele 2?
					if r[k] == a2[k]:
						LH_a2 *= 1 - eps
					else:
						LH_a2 *= eps / 3
				# switching to log-likelihood here
				tempLLH += log(0.5 * (LH_a1 + LH_a2))
			# save LLH for ind, genos
			indLLH[i,j] = tempLLH

	# MLE (EM algorithm) for allele frequency
	af = np.repeat(1/len(allMH), len(allMH)) # start off w/ equal frequencies
	maxIter = 1000
	tolerance = .0001
	L2 = log(2) # saving to avoid repeated computation by interpreter
	lastLLH = 0
	# initialize some variables
	priorZ = np.zeros(genos.shape[0]) # prior prob of all genos
	genoProb = np.zeros((indLLH.shape[0], genos.shape[0])) # llh and then posterior of individuals(rows) x genotypes(columns)
	for numIter in range(0, maxIter):
		# Calculate P(Z|reads, af)
		# first calculate log of prior prob of each genotype
		afLog = np.log(af)
		for i in range(0, genos.shape[0]):
			if genos[i,0] == genos[i,1]: # homozygous
				priorZ[i] = 2 * afLog[genos[i,0]]
			else:
				priorZ[i] = L2 + afLog[genos[i,0]] + afLog[genos[i,1]]
		# now calculate P(Z| reads, af)
		# and calculate likelihood of total model
		totalLLH = 0
		for i in range(0, indLLH.shape[0]):
			# calc (relative) probability of each genotype
			# prior * likelihood
			genoProb[i,:] = priorZ + indLLH[i,:]
			temp_lse = lse(genoProb[i,:])
			totalLLH += temp_lse
			# now normalize
			genoProb[i,:] = np.exp(genoProb[i,:] - temp_lse)
		
		# decide whether to break or not
		if (totalLLH - lastLLH) < tolerance and numIter > 0:
			break
		lastLLH = totalLLH
		
		# Maximization of LLH for af (calculate af given P(Z))
		# first calc number of each genotype (fractional)
		genoProb2 = np.sum(genoProb, 0)
		# now allele frequency
		af = np.zeros(len(allMH))
		for i in range(0, genoProb2.shape[0]):
			af[genos[i,0]] += genoProb2[i]
			af[genos[i,1]] += genoProb2[i]
		# and normalize
		af = af / (2 * indLLH.shape[0])
		# end for loop for EM
	
	# calculate He from allele frequency and return
	He = [str(1 - np.sum(np.power(af, 2))), str(numReads), str(numInds)]
	
	# return allele frequencies if indicated
	if alleleFreq:
		# allele1:freq,allele2:freq,...
		afStr = [allMH[i] + ":" + str(af[i]) for i in range(0, len(af)) if af[i] >= minAF]
		He += [",".join(afStr)]
	
	return He


def Main():
	
	# defaults
	popMap = None # -m, maps bam files (individuals) to populations. bamFilePath \t popName
	snpPos = None # -s, snp positions Chr \t Position (1-based) MUST BE SORTED SMALLEST TO LARGEST
	wS = 60 # -w, window size
	maxSNPs = 8 # -ms, maximum number of SNPs in a window (if exceeded, window is skipped)
	HeOut = "He_mh.txt" # -o, output file name
	af = False # -af, whether to output allele frequencies instead of He
	minAF = 0.0001 # -minAF minimum allele frequency to report allele frequency for
	epsTemplate = 0.01 # -eps probability a base in the template is wrong
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
		elif sys.argv[flag] == "-af":
			flag += 1
			af = True
		elif sys.argv[flag] == "-minAF":
			flag += 1
			minAF = float(sys.argv[flag])
		elif sys.argv[flag] == "-eps":
			flag += 1
			epsTemplate = float(sys.argv[flag])
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
		
	# run sliding window across target SNP file
	with open(snpPos, "r") as snpLocations, open(HeOut, "w") as HeOutFile:
	
		# write header on output file
		lineOut = ["Chr", "Pos"] + masterPop + ["NumReads_" + p for p in masterPop] + ["NumInds_" + p for p in masterPop]
		if af:
			lineOut += ["AlleleFreq_" + p for p in masterPop]
		HeOutFile.write("\t".join(lineOut) + "\n")
		
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
					tempRes = [calcEM(reads[pop], cur, af, minAF, epsTemplate) for pop in masterPop]
					# and write to output
					lineOut = [curChr, ",".join([str(x) for x in cur])] + [i[0] for i in tempRes] + [i[1] for i in tempRes] + [i[2] for i in tempRes]
					if af:
						lineOut += [i[3] for i in tempRes]
					HeOutFile.write("\t".join(lineOut) + "\n")
				
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
			
		# while loop will end with the last window un-evaluated
		# evaluate last window
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
			tempRes = [calcEM(reads[pop], cur, af, minAF, epsTemplate) for pop in masterPop]
			# and write to output
			lineOut = [curChr, ",".join([str(x) for x in cur])] + [i[0] for i in tempRes] + [i[1] for i in tempRes] + [i[2] for i in tempRes]
			if af:
				lineOut += [i[3] for i in tempRes]
			HeOutFile.write("\t".join(lineOut) + "\n")

		# print some summary information
		print("number of windows: ", numWindows)
		print("Distribution of SNPs per window")
		print("SNPs ", "NumberOfWindows")
		for k in sorted(list(numSnpsPerWindow)):
			print(k, numSnpsPerWindow[k])


# This iteratively runs the EM while dropping alleles to account for 
# windows with lots of SNPs (too many to efficiently consider all possible haplotypes
# 
# function to run EM algorithm for estimating either 
# expected heterozygosity or allele frequencies
# of a microhap from low coverage sequencing data
# within ONE population
# @param reads a list of iterators. Each iterator (pysam) is for reads for one individual
# @param pos a list of positions of SNPs (0-based, on reference) in the window
# @param alleleFreq boolean to indicate whether to output allele frequencies (True) or just expectedHet (False)
# @param minAF minimum allele frequency to keep an allele in the analysis
# @param maxNumHaplotypes maximum number of haplotypes to consider in one round. Note that if
#    adding one more position exceeds this value, "NA" will be returned
# @return either the expected heterozygosity as a float or (if not enough information to estimate) the string "NA"
# def iterEM(reads, pos, alleleFreq, minAF):

	
if __name__ == '__main__':
	Main()
