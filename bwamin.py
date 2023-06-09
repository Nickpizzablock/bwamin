import os, argparse, sys
from pyfaidx import Fasta
import align, sambuild, bwt

# Parser
parser = argparse.ArgumentParser(description='minimum bwa')

# Tool options
parser.add_argument('--index', action="store_true", help='index a fasta file')
parser.add_argument('--mem', action="store_true", help='index a fasta file')
parser.add_argument('--bwt', action="store_true", help='enables Burrow-Wheeler Transform search')
parser.add_argument('--sw', action="store_true", help='enables Smith-Waterman search')

# File options
parser.add_argument('--fasta', '--fa', type=str, help="fasta input")
parser.add_argument('--fastq', '--fq', type=str, help="fastq input")

# SW Mem options
parser.add_argument('-A', type=int,	default=1, help="Matching score. [1]")
parser.add_argument('-B', type=int,	default=4, help="Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4]")
parser.add_argument('-O', type=int,	default=6, help="Gap open penalty. [6]")
parser.add_argument('-E', type=int,	default=1, help="Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1]")

# Other
parser.add_argument('-out', type=str, default='zenith', help="Specify output file")
args = parser.parse_args()

# Nothing is selected
optionsList = [args.index, args.mem, args.bwt, args.sw]
if optionsList.count(True) == 0:
    print('Please add arguments. Use \'python bwamin.py -h\' for help.', file=sys.stderr)
    sys.exit(1)

# Check if only 1 option is selected 
optionsList = [args.index, args.mem]
if optionsList.count(True) != 1:
    print('ERROR: Pick one: index or mem', file=sys.stderr)
    sys.exit(1)
    

# Check if only 1 option is selected 
optionsList = [args.bwt, args.sw]
if optionsList.count(True) != 1:
    print('ERROR: Pick one: bwt or sw', file=sys.stderr)
    sys.exit(1)

# sw does not have to be indexed
optionsList = [args.index, args.sw]
if optionsList.count(True) == 2:
    print('ERROR: cannot index with sw, use --mem instead', file=sys.stderr)
    sys.exit(1)

# File path checking
if not os.path.exists(args.fasta):
    print("ERROR: Fasta file invalid", file=sys.stderr)
    sys.exit(1)
else:
    faFile = args.fasta

if args.mem and not os.path.exists(args.fastq):
    print("ERROR: Fastq file invalid", file=sys.stderr)
    sys.exit(1)
else:
    fqFile = args.fastq


# BWT option
if args.index and args.bwt:
    try:
        bwtindex = open(faFile + '.myindex','w')
    except IOError:
        print("ERROR: could not find " + faFile + '.myindex', file=sys.stderr)
        print("Please run --index first")
        sys.exit(1)
    faOut = Fasta(faFile)
    for j in faOut.keys():
        bwtindex.write('>' + j + '\n')
        bwtindex.write(bwt.bwt(str(faOut[j])) + '\n')
    bwtindex.close()
    exit()
# SW option
else:
    # Reading Fasta and fastq
    faOut = Fasta(faFile)

# Read the fq file into fqOut
fqOut = align.sortFqFile(fqFile)

# Stating Alignment settings
match = args.A
mismatch = args.B
indel = args.O
gapPenalty = args.E
print('Alignment Settings')
print('match weight: ' + str(match))
print('mismatch weight: ' + str(mismatch))
print('indel weight: ' + str(indel))
print('gapPenalty: ' + str(gapPenalty))


# Try output file
if '.sam' not in args.out:
    args.out = args.out + '.sam'
try:
    zenith = open(args.out,'w')
except IOError:
    print('Output file does\'t want to be opened. Check chmod', file=sys.stderr)

# Making the header on zenith
# https://samtools.github.io/hts-specs/SAMv1.pdf
for i in faOut.keys():
    zenith.write(sambuild.sq(i, len(faOut[i]), ''))
zenith.write(sambuild.hd('1.5', 'unsorted', 'query'))
command = 'python'
for i in sys.argv:
    command = command + ' ' + i
zenith.write(sambuild.pg('bwamin', 'bwamin', '0.001-alpha', command))

# For each read, look at each chromsome and find the best score
# For bwt
if args.bwt:
    # Read corresponding index for bwt
    faDict = align.alignGenome(faFile + '.myindex')
    for i in fqOut:
        foundMatch = False
        for j in faDict:
            leftpos = []
            leftpos = bwt.find(faDict[j], fqOut[i][0])
           
            # Perfect match
            quality = 60
            if leftpos != None and len(leftpos) == 1:
                flag = 0
                pos = leftpos[0]
                
                # Exact matching only
                # flag = 0 because not checking reverse
                # Pos by 0-index
                foundMatch = True
                zenith.write(sambuild.readToString(i.split(' ', 1)[0].strip(), flag, j, pos, quality, str(len(fqOut[i][0])) + 'M', "rnext", "pnext", "tlen", fqOut[i][0], fqOut[i][1]))
                break
            else:
                #try reverse string
                leftpos = bwt.find(faDict[j], fqOut[i][0][::-1])
                
                # Perfect match on reverse
                if leftpos != None and len(leftpos) == 1:
                    flag = 16               # bitwise code for reverse
                    pos = leftpos[0] + 1    # +1 to be 1 index
                    foundMatch = True
                    zenith.write(sambuild.readToString(i.split(' ', 1)[0].strip(), flag, j, pos, quality, str(len(fqOut[i][0])) + 'M', "rnext", "pnext", "tlen", fqOut[i][0], fqOut[i][1]))
                    break
        # Could not find perfect match          
        if not foundMatch:
            flag = 4
            pos = 0
            zenith.write(sambuild.readToString(i.split(' ', 1)[0].strip(), flag, '*', pos, 0, '*', "*", 0, 0, fqOut[i][0], fqOut[i][1]))
    zenith.close()
    exit()

# SW 
for i in fqOut:
    bestAlignments = {}
    for j in faOut.keys():
        bestAlignments[j] = align.betterswalgo(str(faOut[j]), fqOut[i][0], match, mismatch, indel)
        # bestAlignments[j] = align.nwAlgo(j, faOut[j], i, fqOut[i], match, mismatch, indel, gapPenalty)

        # Not using quality score but matching score
        # Score max 60 min 0
        # bestAlignments[j][0]  = bestAlignments[j][0] #/ len(fqOut[i][0])
        if bestAlignments[j][0] < 0:
            bestAlignments[j][0] = 0
        if bestAlignments[j][0] > 60:
            bestAlignments[j][0] = 60
        # Note: flag = 0 since only detects reads
    # Find highest score alignment
    key, highest = align.maxAlignmentDict(bestAlignments)
    flag = 0
    quality = highest[0]

    zenith.write(sambuild.readToString(i.strip(), flag, key, highest[2], quality, highest[3], '*', 0, 0, fqOut[i][0], fqOut[i][1]))
zenith.close()
