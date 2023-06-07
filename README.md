# bwamin (CSE185 Project Deno)
Ramen with a 'bw'. This is an demonstation project which implements an aligner that takes in `.fa `and `.fq` files and aligns them to output a `.sam` file. Two algorithms (Burrows-Wheeler Transform and Smith-Waterman Algorithm) are presented with their pros and cons. This program will compete against `bwa index` and `bwa min` seen in class. The bwa Github Repo can be accessed here: https://github.com/lh3/bwa

## Installation Procedure
Note: Make sure to use `pip3` and `python3` during this process if applicable
(usually on Macbooks)

### Getting bwamin.py to work
Note2: Xcode is needed on Macs. This process takes very long.
Search Xcode on the App Store and install there.

- Install Python

Python can be installed here: https://www.python.org/downloads/

- Install the 'pyfaidx' library using:
```
# For Mac
pip3 install pyfaidx

# For Other OS
pip install pyfaidx
```
### Benchmarking and Analysis
- Install Anaconda for benchmarking

Anaconda can be installed here: https://www.anaconda.com/download 

- Install BWA and pyfaidx for benchmarking using conda (in mac or linux environment)
```
conda install -c bioconda bwa
conda install -c bioconda pyfaidx
```
- Use IGV

It can be downloaded here:

or just use the webapp. Anyways, the screenshots will be in the notebook.

- Also get pandas and seaborn to view Jupyter Notebooks
```
pip install pandas
pip install seaborn
```

## Basic Usage
An example is the following:
```
# For Mac
python3 bwamin.py [options] [--fa ref.fa] [--fq read.fq]

# For Other OS
python bwamin.py [options] [--fa ref.fa] [--fq read.fq]
```

To index an `.fasta` file:
```
python3 bwamin.py --index --bwt --fa testfastring.fa
```

Then, to get the `.sam` file:
```
python3 bwamin.py --mem --bwt --fa testfastring.fa testfqstring.fq
or
python3 bwamin.py --mem --sw --fa testfastring.fa testfqstring.fq
```
where `testfastring.fa` and `testfqstring.fq` can be replaced with your own `.fa` and `.fq` files.

## Limitations
### For Smith-Waterman options (--sw)
`--sw` does need to index the file, just go straight to `--mem`.

This gives the best case found for each read, however this takes a lot of memory.
This should not be used with large file sizes.

### For Burrows-Wheeler options (--bwt)
Currently, this only checks for exact matches from .fq onto .fa if there is at least one exact match.

This could be used for larger file sizes but limited types of reads.

### For all options
Currently, does not support paired end reads on all options. MAPQ is not implemented well.


## Bwamin Options
### Valid Option Combinations
```
--index --bwt --fa
--mem --bwt --fa --fq
--mem --sw [-A -B -O -E] --fa --fq
```
### Option Details

- `--fasta FILE`, `--fa FILE`: a `.fasta` file, reference genome required for `--index` and `--mem`
- `--fastq FILE`, `--fq FILE`: a `.fasta` file, reads to search in reference genome
- `--index`: stores True, enables index mode for `--bwt`
- `--mem`: stores True, enables mem mode for `--bwt` or `--sw`
- `--bwt`: stores True, enables Burrows-Wheeler Transform 
- `--sw`: stores True, enables Smith-Waterman
- `-A INT`: default 1, matching score
- `-B INT`: default 4, mismatch penalty
- `-O INT`: default 6, gap open penalty
- `-E INT`: default 1, gap extension penalty
- `-out FILE`: default zenith, specify sam file out

## File Format
File formats are the same from BWA but in `zenith.sam`. Do not use the `>` symbol for getting the output.

Bwamin outputs a `.sam` file, different from `bwa mem`. The documentation of the file structure can be found here https://samtools.github.io/hts-specs/SAMv1.pdf

## Contributors
This repo was made by me, Neo Torres (aka Nickpizzablock). I drew inspiration from `bwa`, `Smith-Waterman`, and `Needleman-Wunsch`.

If I made any silly mistakes, @ me plz.


![cat eating bwamin](https://media4.giphy.com/media/Fj0MaDHcLycOk/giphy.gif)

cat eating ~ramen~ bwamin
