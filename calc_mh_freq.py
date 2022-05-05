#!/usr/bin/env python3

# load libraries
import sys
import pysam
from math import log, exp
import numpy as np
import numba
import random

# log of the sum of the exponents
# for the individual based function, this function is a significant use of computation time
# @param a numpy array
@numba.jit(nopython=True)
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
@numba.jit(nopython=True)
def probSubErr(eps, e2):
	# template wrong, sequencer correct + 
	#   template wrong, sequencer wrong and didn't randomly call correct +
	#   template right, sequencer wrong
	return (eps * (1 - e2)) + (0.6666666667 * eps * e2) + ((1 - eps) * e2)

# pull the bases and quality scores at each pos 
# from a set of alignments for ONE bam file (typically an individual or pool)
# while maintaing phase from the (optionally paired) reads
# @param ind an iterator (pysam) for reads for one bam file (often one individual)
# @param pos a list of 0-based positions (in the reference) to pull base calls for
# @parm readGroups True to separate into individuals by read groups, False to not
# @param return a list of [listOne, listTwo] where listOne is a list of lists 
#    representing basecalls within each read, e.g. [[A, C, N], [A, C, A], [G, C, A]] with each
#    read being in the same order as pos. N can either be N present in the read OR the read did 
#    not cover that position. listTwo is similar but has Q scores (converted to int 
#    and 33 substracted already by pysam)
def pullReadsQualsPotenAlleles(ind, pos, readGroups, headerRG = None):
	maxI = len(pos) # used several times below
	# pull relevant sites and Qual scores from reads for each individual
	readNames = []
	tempInd = []
	tempQual = []
	rg = []
	for r in ind: # for each read
		# get base calls for SNPs
		j = r.reference_start # position in the reference, 0 based
		k = 0 # position in the query
		c = 0 # position in the CIGAR string
		# determine if it's the pair of an earlier read
		pair = False
		hasABaseCall = False # Whether the read is informative for >= 1 target SNP
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
							# skip if probability base is wrong is higher than existing
							break
						else:
							if base != "N":
								hasABaseCall = True
							tempRead[i] = base
							tempReadQual[i] = qual
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
		if not hasABaseCall:
			continue # skip if no base calls for target SNPs
		if pair:
			tempInd[matePos] = tempRead
			tempQual[matePos] = tempReadQual
		else:
			tempInd += [tempRead]
			tempQual += [tempReadQual]
			readNames += [r.query_name]
			# record read group if needed
			if readGroups:
				tempRG = None
				for tag in r.get_tags():
					if tag[0] == "RG":
						tempRG = tag[1]
						break
				if tempRG is None:
					raise RuntimeError("read found without read group")
				rg += [headerRG[tempRG]]
	
	# sort read groups if needed
	if readGroups:
		tempIndDict = {}
		tempQualDict = {}
		for i in range(0, len(tempInd)):
			tempIndDict[rg[i]] = tempIndDict.get(rg[i], []) + [tempInd[i]]
			tempQualDict[rg[i]] = tempQualDict.get(rg[i], []) + [tempQual[i]]
		tempInd = []
		tempQual = []
		for k in tempIndDict.keys():
			tempInd += [tempIndDict[k]]
			tempQual += [tempQualDict[k]]
	
	return [tempInd, tempQual]

# wonder if this would be much faster by coercing it to a numpy matrix or pandas df and
# then getting unique elements?
# making a list of lists in the order of 
# pos with lists being the unique bases called at that position (meant to provide
# a set of all possible alleles detected at each position)
# @param baseCalls a list of the first output of pullReadsQualsPotenAlleles for multiple individuals
#   with each entry being an individual
# @return a list of lists as described above
def SNPpossibleAlleles(baseCalls):
	# pulling number of SNPs from the first read of the first individual
	# assumes all reads/inds were assesed for all SNPs, which should be the case
	# also assumes all individuals have >= 1 read, which should also be the case
	if len(baseCalls) == 0:
		return [[]]
	maxI = len(baseCalls[0][0])
	alleles = [[] for x in range(0, maxI)]
	for ind in baseCalls: # for each ind
		for r in ind: # for each read
			for i in range(0, maxI): # for each SNP
				if r[i] != "N" and r[i] not in alleles[i]:
					alleles[i] += [r[i]]
	return alleles

# @param allMH a list of all alleles
# @return a numpy array (matrix) w/ one row per genotype
#    and two columns (one per allele) with ints as values
#    indicating the index of the allele within allMH (0-based)
def makeAllGenos(allMH):
	genos = []
	for i in range(0, len(allMH)):
		for j in range(i, len(allMH)):
			genos += [[i, j]] # allele1, allele2 INDICIES WITHIN allMH
	return np.array(genos)

# this is the EM algorithm run within `iterEM()` but implemented
# as a function so that it can use numba for speed
@numba.jit(nopython=True)
def indivEM(genos, indLLH, af, maxIter, tolerance, numAlleles):
	L2 = log(2) # saving to avoid repeated calculation by interpreter, not sure if this is needed with JIT though?
	# initialize some variables
	lastLLH = 0
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
		af = np.zeros(numAlleles)
		for i in range(0, genoProb2.shape[0]):
			af[genos[i,0]] += genoProb2[i]
			af[genos[i,1]] += genoProb2[i]
		# and normalize
		af = af / (2 * indLLH.shape[0])
		# end for loop for EM
	return af

# this is the EM algorithm run within `poolIterEM()` but implemented
# as a function so that it can use numba for speed
@numba.jit(nopython=True)
def poolEM(readLH, maxIter, tolerance, af):
	# initialize some variables
	# note we're working in likelihood except for totalLLH
	# since we're working with individual reads, don't have to worry about underflow
	# unless there are 100+ SNPs or crazy low error rates (like 1e-20)
	lastLLH = 0
	alleleProb = np.zeros((readLH.shape[0], readLH.shape[1])) # lh and then posterior of reads(rows) x alleles(columns)
	for numIter in range(0, maxIter):
		# now calculate P(Z| read, af)
		# and calculate likelihood of total model
		totalLLH = 0
		for i in range(0, readLH.shape[0]):
			# calc (relative) probability of each allele
			# prior * likelihood
			alleleProb[i,:] = af * readLH[i,:]
			temp_sum = np.sum(alleleProb[i,:])
			totalLLH += log(temp_sum)
			# now normalize
			alleleProb[i,:] = alleleProb[i,:] / temp_sum
		
		# decide whether to break or not
		if (totalLLH - lastLLH) < tolerance and numIter > 0:
			break
		lastLLH = totalLLH
		
		# Maximization of LLH for af (calculate af given P(Z))
		# calc number of each allele (fractional) and normalize
		af = np.sum(alleleProb,0) / readLH.shape[0]
		# end for loop for EM
	return af

# calculate LLH for individual genotypes
# as a function, again, b/c python is _slow_
# and I need to speed things up
@numba.jit(nopython=True)
def indGeno(indReads, indQuals, genos, allMH, epsTemplate):
	indLLH = np.zeros(genos.shape[0])
	for j in range(0, genos.shape[0]): # for each genotype
		a1 = allMH[genos[j,0]] # character strings
		a2 = allMH[genos[j,1]] # character strings
		tempLLH = 0
		for m in range(0, indReads.shape[0]): # for each read
			r = indReads[m,:]
			q = indQuals[m,:]
			LH_a1 = 1 # these are likelihoods, NOT log
			LH_a2 = 1
			for k in range(0, len(a1)): # for each SNP pos in the allele
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
		indLLH[j] = tempLLH
	return indLLH

# calculate LH for reads | individual alleles
# as a function, again, b/c python is _slow_
# and I need to speed things up
@numba.jit(nopython=True)
def calcReadLH(indReads, indQuals, allMH, epsTemplate):
	readLH = np.zeros((indReads.shape[0], allMH.shape[0])) # indices: read, allele; value: likelihood (NOT log)
	for i in range(0, indReads.shape[0]): # for each read
		r = indReads[i,:]
		q = indQuals[i,:]
		for j in range(0, allMH.shape[0]): # for each allele
			a1 = allMH[j] # character string
			LH_a1 = 1 # likelihood, NOT log
			for k in range(0, len(a1)): # for each SNP pos in the allele
				# is it missing? then skip
				if r[k] == "N":
					continue
				# calc prob of error
				eps = probSubErr(epsTemplate, q[k])
				# does it match the allele
				if r[k] == a1[k]:
					LH_a1 *= 1 - eps
				else:
					LH_a1 *= eps / 3
			# switching to log-likelihood here
			# save LLH for read, allele
			readLH[i,j] = LH_a1
	return readLH

# This iteratively runs the EM while dropping alleles to account for
# windows with lots of SNPs (too many to efficiently consider all possible haplotypes
# 
# function to run EM algorithm for estimating either 
# expected heterozygosity or allele frequencies
# of a microhap from low coverage sequencing data
# within ONE population
# @param reads a list of iterators. Each iterator (pysam) is for reads for one individual
# @param pos a list of positions of SNPs (1-based, on reference) in the window
# @param alleleFreq boolean to indicate whether to output allele frequencies (True) or just expectedHet (False)
# @param minAF minimum allele frequency to keep an allele in the analysis
# @param epsTemplate the probability that a base in the template is physically wrong (prior to sequencing)
# @param maxNumHaplotypes maximum number of haplotypes to consider in one round. Note that if
#    adding one more position exceeds this value, "NA" will be returned
# @param maxSNPsAdd maximum number of SNPs to add in one iteration of the algorithm
# @return either the expected heterozygosity as a float or (if not enough information to estimate) the string "NA"
def iterEM(reads, pos, alleleFreq, minAF, epsTemplate, maxNumHaplotypes, maxSNPsAdd, maxNumReads, rSeed, readGroups, headerRG = None):
	# set random seed
	random.seed(rSeed)
	# change pos to 0-based
	pos = [x - 1 for x in pos]
	# pull basecalls and quality scores at each position for each individual
	indReads = []
	indQuals = []
	for ind in reads:
		temp = pullReadsQualsPotenAlleles(ind, pos, readGroups, headerRG)
		if readGroups:
			# either the whole thing is empty, [[],[]] or there is all entries are non-zero length
			# b/c only read groups that are present are considered
			if len(temp[0]) > 0:
				for i in range(0, len(temp[0])): # fro each individual
					# subsample if needed
					if(len(temp[0][i]) > maxNumReads):
						toKeep = random.sample(range(0, len(temp[0][i])), maxNumReads)
						temp[0][i] = [temp[0][i][x] for x in toKeep]
						temp[1][i] = [temp[1][i][x] for x in toKeep]
					indReads += [temp[0][i]]
					indQuals += [temp[1][i]]
		else:
			# only save if any reads were present
			if len(temp[0]) > 0:
				# subsample if needed
				if(len(temp[0]) > maxNumReads):
					toKeep = random.sample(range(0, len(temp[0])), maxNumReads)
					temp[0] = [temp[0][x] for x in toKeep]
					temp[1] = [temp[1][x] for x in toKeep]
				indReads += [temp[0]]
				indQuals += [temp[1]]
	# now determine list of potential alleles at each pos
	allSNPAlleles = SNPpossibleAlleles(indReads)
	
	# calculate some summmary numbers
	numReads = np.array([len(i) for i in indReads])
	numInds = len(numReads) # number of inds with >= 1 read
	numReads = np.sum(numReads) # total number of reads
	
	# make sure there is enough data to continue with the estimates
	skip = False
	if numReads < 1:
		skip = True	# skip if no reads
	
	for d in allSNPAlleles:
		if len(d) == 0: # make sure there is one or more read covering each position
			# recommend filtering loci like this out prior to running this algorithm
			# not going to estimate anything for this window since no data for 1+ position
			skip = True
			break
	if skip:
		lineOut = ["NA", str(numReads), str(numInds)]
		if alleleFreq:
			lineOut += ["NA"]
		return lineOut

	nextPosToAdd = 0# next SNP to add to estimation routine
	allMH = [""]
	af = None
	maxIter = 1000 # EM parameter
	tolerance = .0001 # EM parameter
	indReads = [np.array(x) for x in indReads]
	indQuals = [np.array(x) for x in indQuals]
	while nextPosToAdd < len(allSNPAlleles):
		# detemine microhap alleles to consider
		for i in range(0, maxSNPsAdd):
			# add the next SNP (at least one)
			allMHnew = []
			for a in allSNPAlleles[nextPosToAdd]:
				allMHnew += [x + a for x in allMH]
			allMH = allMHnew
			nextPosToAdd += 1
			# if at the end OR adding the next SNP would put 
			# the number of haplotypes over the limit, move to estimation
			if nextPosToAdd >= len(allSNPAlleles) or (len(allMH) * len(allSNPAlleles[nextPosToAdd])) > maxNumHaplotypes:
				break
		
		# if the addition of one SNP violates the limit on number of 
		# haplotypes, then can't estimate this window
		if len(allMH) > maxNumHaplotypes:
			lineOut = ["NA", str(numReads), str(numInds)]
			if alleleFreq:
				lineOut += ["NA"]
			return lineOut
				
		# now we have the list of haplotypes (alleles), so we estimate
		
		# calculate log-likelihoods for each genotype for each individual
		genos = makeAllGenos(allMH) # determine genotypes to consider
		indLLH = np.zeros((len(indReads), genos.shape[0])) # indices: individual, genotype; value: log-likelihood
		allMH_np = np.array(allMH) # to pass to numba function
		for i in range(0, len(indReads)): # for each individual
			indLLH[i,:] = indGeno(indReads[i], indQuals[i], genos, allMH_np, epsTemplate)
		
		# MLE (EM algorithm) for allele frequency
		if af is None:
			af = np.repeat(1/len(allMH), len(allMH)) # start off w/ equal frequencies
		else:
			afOld = af / (len(allMH) / len(af)) # dividing by number of alleles each allele was split into
			af = []
			oldAlleleLength = len(allMHold[0])
			for a in allMH:
				for i in range(0, len(allMHold)):
					if allMHold[i] == a[0:oldAlleleLength]:
						af += [afOld[i]]
						break
			af = np.array(af)
		
		# run EM to estimate allele frequencies
		# written as a function to speed up with numba
		af = indivEM(genos, indLLH, af, maxIter, tolerance, len(allMH))
		
		# now we either extract alleles that are frequent enough to keep and 
		# start the next iteration,
		# or return the results - note no trimming based on minAF if not
		# running another iteration
		if nextPosToAdd < len(allSNPAlleles):
			# there will be another round
			allMHnew = []
			newAF = []
			for i in range(0, len(af)):
				if af[i] >= minAF:
					newAF += [af[i]]
					allMHnew += [allMH[i]]
			af = np.array(newAF)
			af = af / np.sum(af) # normalize for starting values of next round
			allMHold = allMHnew # this is used to identify af starting values in the new round
			allMH = allMHnew # this is used to build new alleles
	# end of iterative while loop
	
	# calculate He from allele frequency and return
	He = [str(1 - np.sum(np.power(af, 2))), str(numReads), str(numInds)]
	
	# return allele frequencies if indicated
	if alleleFreq:
		# allele1:freq,allele2:freq,...
		afStr = [allMH[i] + ":" + str(af[i]) for i in range(0, len(af))]
		He += [",".join(afStr)]
	
	return He

# This is for data that is all one pool of a large number of individuals with
# no barcodes (i.e., no way to assign reads to individuals)
# This iteratively runs the EM while dropping alleles to account for
# windows with lots of SNPs (too many to efficiently consider all possible haplotypes
# 
# function to run EM algorithm for estimating either 
# expected heterozygosity or allele frequencies
# of a microhap from pooled sequencing data
# within ONE population
# @param reads a list of iterators. Each iterator (pysam) is for reads for one bam file, but they
#    are combined and treated as one input pool
# @param pos a list of positions of SNPs (1-based, on reference) in the window
# @param alleleFreq boolean to indicate whether to output allele frequencies (True) or just expectedHet (False)
# @param minAF minimum allele frequency to keep an allele in the analysis
# @param epsTemplate the probability that a base in the template is physically wrong (prior to sequencing)
# @param maxNumHaplotypes maximum number of haplotypes to consider in one round. Note that if
#    adding one more position exceeds this value, "NA" will be returned
# @param maxSNPsAdd maximum number of SNPs to add in one iteration of the algorithm
# @return either the expected heterozygosity as a float or (if not enough information to estimate) the string "NA"
def poolIterEM(reads, pos, alleleFreq, minAF, epsTemplate, maxNumHaplotypes, maxSNPsAdd, maxNumReads, rSeed):
	# set random seed
	random.seed(rSeed)
	# change pos to 0-based
	pos = [x - 1 for x in pos]
	# pull basecalls and quality scores at each position for each individual
	indReads = []
	indQuals = []
	for ind in reads:
		temp = pullReadsQualsPotenAlleles(ind, pos, False)
		# only save if any reads were present
		if len(temp[0]) > 0:
			# note the difference here vs barcoded function
			# Here we are saving these as one long list of reads, not a list of individuals
			indReads += temp[0]
			indQuals += temp[1]
	# subsample if needed
	if(len(indReads) > maxNumReads):
		toKeep = random.sample(range(0, len(indReads)), maxNumReads)
		indReads = [indReads[x] for x in toKeep]
		indQuals = [indQuals[x] for x in toKeep]
	
	# calculate some summmary numbers
	numReads = len(indReads)
	
	# make sure there is enough data to continue with the estimates
	skip = False
	if numReads < 1:
		skip = True	# skip if no reads
	else:
		# now determine list of potential alleles at each pos
		# wrapping in an extra list to match input function expects
		# will throw error if run when numReads = 0
		allSNPAlleles = SNPpossibleAlleles([indReads])
		for d in allSNPAlleles:
			if len(d) == 0: # make sure there is one or more read covering each position
				# recommend filtering loci like this out prior to running this algorithm
				# not going to estimate anything for this window since no data for 1+ position
				skip = True
				break
	if skip:
		lineOut = ["NA", str(numReads)]
		if alleleFreq:
			lineOut += ["NA"]
		return lineOut

	nextPosToAdd = 0 # next SNP to add to estimation routine
	allMH = [""]
	af = None
	maxIter = 1000 # EM parameter
	tolerance = .0001 # EM parameter
	indReads = np.array(indReads)
	indQuals = np.array(indQuals)
	while nextPosToAdd < len(allSNPAlleles):
		# detemine microhap alleles to consider
		for i in range(0, maxSNPsAdd):
			# add the next SNP (at least one)
			allMHnew = []
			for a in allSNPAlleles[nextPosToAdd]:
				allMHnew += [x + a for x in allMH]
			allMH = allMHnew
			nextPosToAdd += 1
			# if at the end OR adding the next SNP would put 
			# the number of haplotypes over the limit, move to estimation
			if nextPosToAdd >= len(allSNPAlleles) or (len(allMH) * len(allSNPAlleles[nextPosToAdd])) > maxNumHaplotypes:
				break
		
		# if the addition of one SNP violates the limit on number of 
		# haplotypes, then can't estimate this window
		if len(allMH) > maxNumHaplotypes:
			lineOut = ["NA", str(numReads)]
			if alleleFreq:
				lineOut += ["NA"]
			return lineOut
				
		# now we have the list of haplotypes (alleles), so we estimate
		
		# calculate likelihoods for each allele for each read
		readLH = calcReadLH(indReads, indQuals, np.array(allMH), epsTemplate)
		
		# MLE (EM algorithm) for allele frequency
		
		# pick starting values
		if af is None:
			af = np.repeat(1/len(allMH), len(allMH)) # start off w/ equal frequencies
		else:
			afOld = af / (len(allMH) / len(af)) # dividing by number of alleles each allele was split into
			af = []
			oldAlleleLength = len(allMHold[0])
			for a in allMH:
				for i in range(0, len(allMHold)):
					if allMHold[i] == a[0:oldAlleleLength]:
						af += [afOld[i]]
						break
			af = np.array(af)
		
		# run EM (written as a function to use numba)
		af = poolEM(readLH, maxIter, tolerance, af)
	
		# now we either extract alleles that are frequent enough to keep and 
		# start the next iteration,
		# or return the results - note no trimming based on minAF if not
		# running another iteration
		if nextPosToAdd < len(allSNPAlleles):
			# there will be another round
			allMHnew = []
			newAF = []
			for i in range(0, len(af)):
				if af[i] >= minAF:
					newAF += [af[i]]
					allMHnew += [allMH[i]]
			af = np.array(newAF)
			af = af / np.sum(af) # normalize for starting values of next round
			allMHold = allMHnew # this is used to identify af starting values in the new round
			allMH = allMHnew # this is used to build new alleles
	# end of iterative while loop
	
	# calculate He from allele frequency and return
	He = [str(1 - np.sum(np.power(af, 2))), str(numReads)]
	
	# return allele frequencies if indicated
	if alleleFreq:
		# allele1:freq,allele2:freq,...
		afStr = [allMH[i] + ":" + str(af[i]) for i in range(0, len(af))]
		He += [",".join(afStr)]
	
	return He


def Main():
	# defaults
	popMap = None # -m, maps bam files (individuals) to populations. bamFilePath \t popName
	snpPos = None # -s, snp positions Chr \t Position (1-based) MUST BE SORTED SMALLEST TO LARGEST
	wS = 125 # -w, window size
	maxSNPs = 25 # -ms, maximum number of SNPs in a window (if exceeded, window is skipped)
	HeOut = "He_mh.txt" # -o, output file name
	af = False # -af, whether to output allele frequencies He
	minAF = 0.001 # -minAF minimum allele frequency to keep in the model when iterating over a large window
	epsTemplate = 0.01 # -eps probability a base in the template is wrong
	maxNumHaplotypes = None # -maxH maximum number of haplotypes to try to estimate for
	maxSNPsAdd = 1 # -maxS maximum number of SNPs to add in one iteration
	pooledData = False # -pool whether or not to treat data as individuals or as a pool
	maxNumReads = None # -maxR maximum number of reads for a window, either in an individual, or overall if -pool (see defaults below)
	rSeed = 7 # -r random seed for subsampling reads
	readGroups = False # -rg to separate reads in the bam files by read group SM into individuals, ignored if -pool is specified. RG ID and SM must be unique within pops
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
		elif sys.argv[flag] == "-r":
			flag += 1
			rSeed = int(sys.argv[flag])
		elif sys.argv[flag] == "-maxR":
			flag += 1
			maxNumReads = int(sys.argv[flag])
		elif sys.argv[flag] == "-maxH":
			flag += 1
			maxNumHaplotypes = int(sys.argv[flag])
		elif sys.argv[flag] == "-maxS":
			flag += 1
			maxSNPsAdd = int(sys.argv[flag])
		elif sys.argv[flag] == "-af":
			af = True
		elif sys.argv[flag] == "-pool":
			pooledData = True
		elif sys.argv[flag] == "-minAF":
			flag += 1
			minAF = float(sys.argv[flag])
		elif sys.argv[flag] == "-eps":
			flag += 1
			epsTemplate = float(sys.argv[flag])
		elif sys.argv[flag] == "-o":
			flag += 1
			HeOut = sys.argv[flag]
		elif sys.argv[flag] == "-rg":
			readGroups = True
		else:
			print("Error: option", sys.argv[flag], "not recognized.")
			return
		flag += 1
	# end while loop for command line arguments
	
	# different default for pooled vs individual data
	if maxNumHaplotypes is None:
		if pooledData:
			maxNumHaplotypes = 256
		else:
			maxNumHaplotypes = 128
	if maxNumReads is None:
		if pooledData:
			maxNumReads = 100
		else:
			maxNumReads = 40
	
	if popMap is None or snpPos is None:
		print("Error: both -m and -s must be specified")
		return
	if maxSNPsAdd < 1:
		print("Error: -maxS must be 1 or greater")
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
	
	# extract read groups from header(s) if needed
	if not pooledData:
		headerRG = {} # key is popName, value is dictionary {key is RG ID, value is RG SM}
		for p in masterPop:
			if readGroups:
				headerRG[p] = {}
				for tempAlignmentFile in popDict[p]:
					tempHeader = tempAlignmentFile.header.to_dict()["RG"]
					for line in tempHeader:
						headerRG[p][line["ID"]] = line["SM"]
			else:
				headerRG[p] = None
	
	# run sliding window across target SNP file
	with open(snpPos, "r") as snpLocations, open(HeOut, "w") as HeOutFile:
	
		# write header on output file
		lineOut = ["Chr", "Pos"] + masterPop + ["NumReads_" + p for p in masterPop]
		if not pooledData:
			lineOut += ["NumInds_" + p for p in masterPop]
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
					if pooledData:
						tempRes = [poolIterEM(reads[pop], cur, af, minAF, epsTemplate, maxNumHaplotypes, maxSNPsAdd, maxNumReads, rSeed) for pop in masterPop]
						lineOut = [curChr, ",".join([str(x) for x in cur])] + [i[0] for i in tempRes] + [i[1] for i in tempRes]
					else:
						tempRes = [iterEM(reads[pop], cur, af, minAF, epsTemplate, maxNumHaplotypes, maxSNPsAdd, maxNumReads, rSeed, readGroups, headerRG[pop]) for pop in masterPop]
						lineOut = [curChr, ",".join([str(x) for x in cur])] + [i[0] for i in tempRes] + [i[1] for i in tempRes] + [i[2] for i in tempRes]
					# and write to output
					if af:
						lineOut += [i[3 - pooledData] for i in tempRes]
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
			if pooledData:
				tempRes = [poolIterEM(reads[pop], cur, af, minAF, epsTemplate, maxNumHaplotypes, maxSNPsAdd, maxNumReads, rSeed) for pop in masterPop]
				lineOut = [curChr, ",".join([str(x) for x in cur])] + [i[0] for i in tempRes] + [i[1] for i in tempRes]
			else:
				tempRes = [iterEM(reads[pop], cur, af, minAF, epsTemplate, maxNumHaplotypes, maxSNPsAdd, maxNumReads, rSeed, readGroups, headerRG[pop]) for pop in masterPop]
				lineOut = [curChr, ",".join([str(x) for x in cur])] + [i[0] for i in tempRes] + [i[1] for i in tempRes] + [i[2] for i in tempRes]
			# and write to output
			if af:
				lineOut += [i[3 - pooledData] for i in tempRes]
			HeOutFile.write("\t".join(lineOut) + "\n")

		# print some summary information
		print("number of windows: ", numWindows)
		print("Distribution of SNPs per window")
		print("SNPs ", "NumberOfWindows")
		for k in sorted(list(numSnpsPerWindow)):
			print(k, numSnpsPerWindow[k])
	
if __name__ == '__main__':
	Main()
