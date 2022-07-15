# mhFromLowDepSeq
 Identify microhaplotypes from low depth sequencing data
 
 Situation: You have either pooled sequencing data (sequencing data for a bunch of individuals 
 with no barcodes or other means of identifying which reads came from which individual) or you 
 have low coverage sequencing data for a bunch of individuals (not enough to reliably call genotypes).
 
 Desired outcome: Allele frequencices and/or expected heterozygosity for microhaps in the population 
 that was sequenced.
 
 Method: Mixture models fit with EM algorithms to infer allele frequencies (see full description in ...).
 
 Limitations: Can only handle substitution SNPs
 
 # A few notes

 To use the method, the only file you need to download from this repository is `calc_mh_freq.py` 
 (and I might also recommend the `example` directory). The 
 other files are only included because they are relevant to the analyses discussed in ...

 
 Note on speed: I recommend running it separately (and concurrently)
 on each chromosome and/or running it separately for each population. If running populations separately,
 you should be able to join the results together directly (or using the Chr and Pos columns as keys if 
 you want to double check). As long as each population is run with the same `-s` input file, `-w` window 
 size, and `-ms` maximum number of SNPs in a window, these columns should be identical (the same windows 
 evaluated in the same order). In some limited speed testing, it seems that using smaller values of `-maxS`, 
 such as 1, increases speed. It also seems that combining all reads into a single BAM file (and identiifying 
 individuals with read groups, if you have data for individuals) is significantly quicker than having one BAM file 
 for each individual (sorting reads by read group internally is quicker than multiple BAM reading operations). 
 Additionally, using the pool method (`-pool`) seems generally faster (at the expense 
 of model complexity). Also note that since it uses numba, there is some fixed compilation time at the start of 
 the script. For analyses of realistic size, this extra time is quickly made up.
 
 Note on memory: memory usage should be pretty low (a few GB or less), as long as you are using 
 smaller values of `-maxS`, do not have hundreds of individuals (for the individual method), and use reasonable values of `-maxH` and 
 `maxR`. If memory usage is high, one common cause is that you are inputing a large number of BAM files and they have 
 large headers (lots of @SQ from thousands of small contigs and/or an unreasonable number of @PG lines). This happens because the header 
 is stored in memory while the file is opened for reading. You can remedy this by combining 
 the bam files into one per population (mark individuals with @RG). Theoretically, you could also make the header 
 smaller by removing irrelevant information, but only do this if you understand what you're doing and any downstream 
 repurcussions that you may cause, particularly if you use the file(s) for something else.

# Manual

## Quick examples: 

For low coverage data from multiple individuals
```
python3 calc_mh_freq.py -m ./example/filePop.txt -s ./example/snpLoc.txt -o output.txt -af -w 65
```

For pooled data (can't assign reads to individuals)
```
python3 calc_mh_freq.py -m ./example/filePop.txt -s ./example/snpLoc.txt -o outputPool.txt -af -w 65 -pool
```

## Command line options
 
- `-m FILE` a tab delimited text file with _no_ header that maps bam files (pools or individuals, first column) to 
populations (second column). bamFilePath \t popName
- `-s FILE` a tab delimited text file with _no_ header that indicates the positions of substitution SNPs (no indels). First 
column is the chromosome name and second column is the 1-based position. Variants _must_ be sorted by chromosome
and by position, with position sorted from smallest to largest. Chr \t Position
- `-w INT` window size for defining microhaplotypes. default 125
- `-ms INT` maximum number of SNPs in a window to be valid (if exceeded, window is skipped and nothing is written to 
output for this window). default 25
- `-o STRING` output file name. default He_mh.txt
- `-af` include this argument if you want the output file to contain allele frequencies in addition to He. default is to not output allele frequencies
- `-minAF FLOAT` minimum allele frequency to keep an allele in the model when pruning. Note that pruning is _not_ performed after the 
last cycle, so the final result may have alleles with frequencies below this minimum. default 0.001
- `-eps FLOAT` probability a base in the template is wrong, should be > 0 and < 1. default 0.01
- `-maxH INT` maximum number of alleles to try to estimate results for. If this is exceeded, the window is skipped and NA is output. default is 128 for 
low coverage data from individuals and 256 for pooled data (option `-pool` is used)
- `-maxS INT` maximum number of SNPs to add in one iteration. Must be > 0. default 1
- `-pool` include this argument to treat data (even if multiple bam files per population are provided) as a pool. default is to treat data as coming from individuals
with no knowledge of which reads came from which individual
- `-maxR INT` The maximum number of reads within a window to consider for a given individual or, if `-pool` is used, for the entire pool. If more reads
are present, a random subsample is taken. default for individual based analyses: 40 reads per individual; for `-pool` analyses: 100 reads total. 
- `-r INT` The seed to use for random subsampling of reads when a window exceeds `-maxR`. default is 7
- `-rg` include this argument to split reads within _ONE_ BAM file into separate samples based on the read group SM tag. Note that 
  __all reads for a given sample must be in one BAM file__. If you use this option, it is __strongly__ recommended that you only input one BAM 
  file for each population.
