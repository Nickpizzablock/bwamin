# bwamin
Ramen with a 'bw'. This is a student project emulating BWA MEM. 

## Example Commands
`bwa index hg19.fa`
`bwa mem [options] <idxbase> <file1.fq> <file2.fq> > output.sam `

## Files needed:
- `.fa'
- `.fq'
- Find the sequences using ???

# Timeline:
[x] Choose if only bwa mem or all bwa
    - I think just do bwa mem for now
[] Make the algorithm
[] Make the command controls
[] Find out how index works
[] Find out how mem works
[] Extra: Find how gpu works

# Commands (Commands used in class):
[] index         index sequences in the FASTA format
[] mem           BWA-MEM algorithm

    (Commands We didn't use)
    [] fastmap       identify super-maximal exact matches
    [] pemerge       merge overlapping paired ends (EXPERIMENTAL)
    [] aln           gapped/ungapped alignment
    [] samse         generate alignment (single ended)
    [] sampe         generate alignment (paired ended)
    [] bwasw         BWA-SW for long queries (DEPRECATED)

    [] shm           manage indices in shared memory
    [] fa2pac        convert FASTA to PAC format
    [] pac2bwt       generate BWT from PAC
    [] pac2bwtgen    alternative algorithm for generating BWT
    [] bwtupdate     update .bwt to the new format
    [] bwt2sa        generate SA from BWT and Occ

# BWA Commands Help Text
## BWA Index:
```
-a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]
-p STR    prefix of the index [same as fasta name]
-b INT    block size for the bwtsw algorithm (effective with -a bwtsw)[10000000]
-6        index files named as <in.fasta>.64.* instead of <in.fasta>.*
```

## BWA MEM:
Algorithm options:
```
-t INT        number of threads [1]
-k INT        minimum seed length [19]
-w INT        band width for banded alignment [100]
-d INT        off-diagonal X-dropoff [100]
-r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [15]
-y INT        seed occurrence for the 3rd round seeding [20]
-c INT        skip seeds with more than INT occurrences [500]
-D FLOAT      drop chains shorter than FLOAT fraction of the longestoverlapping chain [0.50]
-W INT        discard a chain if seeded bases shorter than INT [0]
-m INT        perform at most INT rounds of mate rescues for each read [50]
-S            skip mate rescue
-P            skip pairing; mate rescue performed unless -S also in use
```

Scoring options:
```
-A INT        score for a sequence match, which scales options -TdBOELU unlessoverridden [1]
-B INT        penalty for a mismatch [4]
-O INT[,INT]  gap open penalties for deletions and insertions [6,6]
-E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
-L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
-U INT        penalty for an unpaired read pair [17]
-x STR        read type. Setting -x changes multiple parameters unlessoverridden [null]
              pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
              ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
              intractg: -B9 -O16 -L5  (intra-species contigs to ref)
```
Input/output options:
```
-p            smart pairing (ignoring in2.fq)
-R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
-H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE[null]
-o FILE       sam file to output results to [stdout]
-j            treat ALT contigs as part of the primary assembly (i.e. ignore<idxbase>.alt file)
-5            for split alignment, take the alignment with the smallest query(not genomic) coordinate as primary
-q            don't modify mapQ of supplementary alignments
-K INT        process INT input bases in each batch regardless of nThreads (forreproducibility) []
-v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]
-T INT        minimum score to output [30]
-h INT[,INT]  if there are <INT hits with score >80.00% of the max score,output all in XA [5,200]
              A second value may be given for alternate sequences.
-z FLOAT      The fraction of the max score to use with -h [0.800000].
              specify the mean, standard deviation (10% of the mean if absent), max
-a            output all alignments for SE or unpaired PE
-C            append FASTA/FASTQ comment to SAM output
-V            output the reference FASTA header in the XR tag
-Y            use soft clipping for supplementary alignments
-M            mark shorter split hits as secondary
-I FLOAT[,FLOAT[,INT[,INT]]]
              specify the mean, standard deviation (10% of the mean if absent), max
              (4 sigma from the mean if absent) and min of the insert size distribution.
              FR orientation only. [inferred]
-u            output XB instead of XA; XB is XA with the alignment score andmapping quality added.
```

Implementation Notes:
bwa index:
- Output files are the following
    - /home/ndtorres/public/genomes/hg19.fa
    - /home/ndtorres/public/genomes/hg19.fa.fai
    - /home/ndtorres/public/genomes/hg19.fa.amb
    - /home/ndtorres/public/genomes/hg19.fa.pac
    - /home/ndtorres/public/genomes/hg19.fa.ann
    - /home/ndtorres/public/genomes/hg19.fa.sa
    - /home/ndtorres/public/genomes/hg19.fa.bwt
