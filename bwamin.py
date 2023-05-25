"""
Notes:
Do mem portion. Prob not do index.

I need to have a package done by friday night


"""


import os, argparse

# Parser
parser = argparse.ArgumentParser(description='minimum bwa')

# Tool options
parser.add_argument('--index', action="store_true", help='index a fasta file')
parser.add_argument('--mem', action="store_true", help='index a fasta file')

# Mem options
parser.add_argument('-t', type=int, default=1, help="number of threads [1]")
args = parser.parse_args()

# Check if only 1 option is selected (-h still works here)
optionsList = [args.index, args.mem]
if optionsList.count(True) != 1:
    print('ERROR: Too many / few arguments')
    raise False


# print(optionsList.count(True))
# parser.add_argument("index", type=int, required=True, help="match val")
# parser.add_argument("mem", type=int, required=True, help="mismatch val")
print('hellowaorld')
