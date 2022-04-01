# mhFromLowDepSeq
 Identify microhaplotypes from low depth sequencing data
 
 Under development. Not ready for release.
 
 Situation: You have either pooled sequencing data (sequencing data for a bunch of individuals 
 with no barcodes or other means of identifying which reads came from which individual) or you 
 have low coverage sequencing data for a bunch of individuals (not enough to reliably call genotypes).
 
 Desired outcome: Allele frequencices and/or expected heterozygosity for microhaps in the population 
 that was sequenced.
 
 Method: Mixture models fit with EM algorithms to infer allele frequencies (will add more here later).
 
 Note on speed: This is not the fastest function. I recommend running it separately (and concurrently)
 on each chromosome and/or running it separately for each population. If running populations separately,
 you should be able to join the results together directly (or using the Chr and Pos columns as keys if 
 you want to double check). As long as each population is run with the same `-s` input file, `-w` window 
 size, and `-ms` maximum number of SNPs in a window, these columns should be identical (the same windows 
 evaluated in the same order). In some limited speed testing, it seems that using smaller values of `-maxS`, 
 such as 1, increases speed. Additionally, using the pool method (`-pool`) seems generally faster (at the expense 
 of model complexity).

# Manual

## Quick examples: 

For low coverage data from multiple individuals
```
python3 calc_mh_freq.py -m indsToPopulations.txt -s known_snps.txt -o outputFile.txt -af -w 65
```

For pooled data (can't assign reads to individuals)
```
python3 calc_mh_freq.py -m bamfileToPopulations.txt -s known_snps.txt -o outputFile.txt -af -w 65 -pool
```

## Command line options
 
- `-m FILE` a tab delimited text file with _no_ header that maps bam files (pools or individuals, first column) to 
populations (second column). bamFilePath \t popName
- `-s FILE` a tab delimited text file with _no_ header that indicates the positions of SNPs. First 
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
- `-maxH INT` maximum number of haplotypes to try to estimate results for. If this is exceeded, the window is skipped and NA is output. default is 128 for 
low coverage data from individuals and 256 for pooled data (option -pool is used)
- `-maxS INT` maximum number of SNPs to add in one iteration. Must be > 0. default 4
- `-pool` include this argument to treat data (even if multiple bam files per population are provided) as a pool. default is to treat data as coming from individuals
with no knowledge of which reads came from which individual

