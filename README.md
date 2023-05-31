# bwamin (SCE185 Project Deno)
Ramen with a 'bw'. This is an demonstation project which implements an aligner that takes in `.fa `and `.fq` files and aligns them to output a `.sam` file. This will go against `bwa index` and `bwa min` seen in class. The bwa Github Repo can be accessed [here.](https://github.com/lh3/bwa) 

# Installation Prosedure
Please install the 'pyfaidx' library using:

```
# For Mac
pip3 install pyfaidx

# For Other OS
pip install pyfaidx
```

Currently, this is a command line tool. No installation for now.

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
python3 bwamin.py testfastring.fa testfqstring.fq
```
where `testfastring.fa` and `testfqstring.fq` can be replaced with your own `.fa` and `.fq` files.

## Note
Under development. Currently, this only checks for exact matches from .fq onto .fa if there is at least one exact match. Currently, does not support paired end reads. Does not output a `.sam ` file.

# Bwamin Options
## Note
Options are in progress.

# File Format
Currently, no file is output.

Bwamin outputs a `.sam` file, the same as `bwa mem`. The documentation of the file structure can be found [here.](https://samtools.github.io/hts-specs/SAMv1.pdf)

# Contributors
This repo was made by me, Neo Torres (aka Nickpizzablock). I drew inspiration from `bwa`, `Smith-Waterman`, and `Needleman-Wunsch`.

If I made any silly mistakes, @ me and submit a pull reqest plz.
