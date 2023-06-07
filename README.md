# bwamin (CSE185 Project Deno)
Ramin with a 'bw'. This is an demonstation project which implements an aligner that takes in `.fa `and `.fq` files and aligns them to output a `.sam` file. This program will compete against `bwa index` and `bwa min` seen in class. The bwa Github Repo can be accessed here: https://github.com/lh3/bwa

# Installation Prosedure
Note: Make sure to use `pip3` and `python3` during this process if applicable
(usually on Macbooks)

Note2: Xcode is needed on Macs. This process takes very long.
Search Xcode on the App Store and install there.

- Install Python
Python can be installed here: https://www.python.org/downloads/
- Install Anaconda for benchmarking
Anaconda can be installed here: https://www.anaconda.com/download 
- Install BWA for benchmarking using conda
```
conda install -c bioconda bwa
```
- Install the 'pyfaidx' library using:
```
# For Mac
pip3 install pyfaidx

# For Other OS
pip install pyfaidx
```
# Basic Usage
An example is the following:
```
# For Mac
python3 bwamin.py [options] [.fa] [.fq]

# For Other OS
python bwamin.py [options] [.fa] [.fq]
```

A test example is:
```
python3 bwamin.py --index --bwt testfastring.fa testfqstring.fq
```
where `testfastring.fa` and `testfqstring.fq` can be replaced with your own `.fa` and `.fq` files.

## Note
Under development. Currently, this only checks for exact matches from .fq onto .fa if there is at least one exact match. Currently, does not support paired end reads. Does not output a `.sam ` file.

# Bwamin Options
## Note
Must choose between `--index` or `--mem` but they don't change the output for now.

# File Format
File formats are DIFFERENT from BWA.


Bwamin outputs a `.sam` file, different from `bwa mem`. The documentation of the file structure can be found [here.](https://samtools.github.io/hts-specs/SAMv1.pdf)

# Contributors
This repo was made by me, Neo Torres (aka Nickpizzablock). I drew inspiration from `bwa`, `Smith-Waterman`, and `Needleman-Wunsch`.

If I made any silly mistakes, @ me plz.
