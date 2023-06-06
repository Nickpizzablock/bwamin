import os, argparse, sys
from pyfaidx import Fasta
import align, sambuild, bwt

# python bwamin.py --index -A 1 -B 1 -O 1 -E 1 short.fa short.fq > output.txt
# python bwamin.py --index -A 1 -B 1 -O 1 -E 1 testfastring.fa testfqstring.fq > output.txt

# python bwamin.py --index --bwt -A 1 -B 1 -O 1 -E 1 .\benchmark\mydata\SRR10769501.fasta.fixed .\benchmark\mydata\SRR10769501.fastq

# Parser
parser = argparse.ArgumentParser(description='minimum bwa')

# Tool options
parser.add_argument('--index', action="store_true", help='index a fasta file')
parser.add_argument('--mem', action="store_true", help='index a fasta file')
parser.add_argument('--bwt', action="store_true", help='enables bwt search')


# My options
parser.add_argument('fasta', type=str, help="fasta input")
parser.add_argument('fastq', type=str, help="fastq input")

# Mem options (Found in bwa mem help)
# parser.add_argument('-t', type=int, default=1, help="number of threads [1]")
# parser.add_argument('-k', type=int,	default=19, help="Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19]")
# parser.add_argument('-w', type=int,	default=100, help="Band width. Essentially, gaps longer than INT will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. [100]")
# parser.add_argument('-d', type=int,	default=100, help="Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. [100]")
# parser.add_argument('-r', type=float, default=1.5, help="Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [1.5]")
# parser.add_argument('-c', type=int,	default=10000, help="Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000]")
# parser.add_argument('-P', help="In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.")
parser.add_argument('-A', type=int,	default=1, help="Matching score. [1]")
parser.add_argument('-B', type=int,	default=4, help="Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4]")
parser.add_argument('-O', type=int,	default=6, help="Gap open penalty. [6]")
parser.add_argument('-E', type=int,	default=1, help="Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1]")
# parser.add_argument('-L', type=int,	default=5, help="Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [5]")
# parser.add_argument('-U', type=int,	default=9, help="Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. [9]")
# parser.add_argument('-p', help="Assume the first input query file is interleaved paired-end FASTA/Q. See the command description for details.")
# parser.add_argument('-R', type=str, help="Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. [null]")
# parser.add_argument('-T', type=int,	default=30, help="Don’t output alignment with score lower than INT. This option only affects output. [30]")
# parser.add_argument('-a', help="Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.")
# parser.add_argument('-C', help="Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.")
# parser.add_argument('-H', help="Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.")
# parser.add_argument('-M', help="Mark shorter split hits as secondary (for Picard compatibility).")
# parser.add_argument('-v', type=int,	default=3, help="Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. When this option takes value 4, the output is not SAM. [3] ")
args = parser.parse_args()

# Check if only 1 option is selected (-h still works here)
optionsList = [args.index, args.mem]
if optionsList.count(True) != 1:
    print('ERROR: Only pick index or mem')
    raise False

# File checking
# Note: Should I add a .fa and .fq checker? bc txt should be fine
if not os.path.exists(args.fasta):
    print("ERROR: Fasta file invalid")
    raise OSError()
else:
    faFile = args.fasta

if not os.path.exists(args.fastq):
    print("ERROR: Fastq file invalid")
    raise OSError()
else:
    fqFile = args.fastq

if args.index and args.bwt:
    bwtindex = open(faFile + '.myindex','w')
    faOut = Fasta(faFile)
    for j in faOut.keys():
        bwtindex.write('>' + j + '\n')
        bwtindex.write(bwt.bwt(str(faOut[j])) + '\n')
    bwtindex.close()
    exit()
else:
    # Reading Fasta and fastq
    faOut = Fasta(faFile)
# for each in faOut.keys():
#     print(each)
#     print(faOut[each])

# exit()
fqOut = align.sortFqFile(fqFile)

# print(faOut)
# exit()
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
# print(sys.argc)
# exit()
bestAlignments = {}
# print(faOut.keys(faOut.keys[0]))
# for i in faOut.keys():
#     print(faOut[i])
# exit()
# For each read, look at each chromsome and find the best score



zenith = open('zenith.txt','w')

# Making the header on zenith
# https://samtools.github.io/hts-specs/SAMv1.pdf
for i in faOut.keys():
    zenith.write(sambuild.sq(i, len(faOut[i]), ''))
zenith.write(sambuild.hd('1.5', 'unsorted', 'query'))
command = 'python'
for i in sys.argv:
    command = command + ' ' + i
zenith.write(sambuild.pg('bwamin', 'bwamin', '0.001-alpha', command))

if args.bwt:
    # bwtindex = open(faFile + '.myindex','r')
    faDict = align.alignGenome(faFile + '.myindex')
    # print(faDict)

    # exit()
    for i in fqOut:
        for j in faDict:
            # print(fqOut[i])
            # exit()
            leftpos = []
            leftpos = bwt.find(faDict[j], fqOut[i][0])
            # if leftpos != None
            if leftpos != None and len(leftpos) == 1:
                flag = 0
                pos = leftpos[0]
                print('found at least 1 exact match')
                #exact matching only
                #flag = 0 because not checking reverse
                #i think leftpos is one off
                zenith.write(sambuild.readToString(i.split(' ', 1)[0].strip(), flag, j, pos+1, "quality", str(len(fqOut[i][0])) + 'M', "rnext", "pnext", "tlen", fqOut[i][0], fqOut[i][1])) # note: you need to put \n
            else:
                #try reverse string
                leftpos = bwt.find(faDict[j], reversed(fqOut[i][0]))
                if leftpos != None and len(leftpos) == 1:
                    flag = 16
                    pos = leftpos[0]
                    print('found at least 1 exact match reversed')
                    #exact matching only
                    #flag = 0 because not checking reverse
                    #i think leftpos is one off
                    zenith.write(sambuild.readToString(i.split(' ', 1)[0].strip(), flag, j, pos+1, "quality", str(len(fqOut[i][0])) + 'M', "rnext", "pnext", "tlen", fqOut[i][0], fqOut[i][1])) # note: you need to put \n
                else:
                    #try the reverse string
                    flag = 4
                    pos = 0
                    zenith.write(sambuild.readToString(i.split(' ', 1)[0].strip(), flag, '*', pos, 0, '*', "*", 0, 0, fqOut[i][0], fqOut[i][1])) # note: you need to put \n
    zenith.close()
    exit()

# Making the reads list on zenith
for i in fqOut:
    # genomeScores = []
    for j in faOut.keys():
        # genomeScores.append(align.alignThem(j, i))
        # if align.alignThem((j, faOut[j]), (i, fqOut[i])):
        #     print("there is a match")
        bestAlignments[j] = align.nwAlgo(j, faOut[j], i, fqOut[i], match, mismatch, indel, gapPenalty)

        # Note: flag = 0 since only single reads
        zenith.write(sambuild.readToString(i, 'flag', j, bestAlignments[j][2], "quality", bestAlignments[j][3], "rnext", "pnext", "tlen", fqOut[i][0], fqOut[i][1])) # note: you need to put \n
        # zenith.write('hi\thi')

        # exit()
    # bestScore = max(genomeScores)
    # samOut.append(bestScore)
# align.nwAlgo()
# Writing the contents out
# with open("output.sam", "w") as writer:
#     for i in samOut:
#         writer.write(i)
# print(optionsList.count(True))
# parser.add_argument("index", type=int, required=True, help="match val")
# parser.add_argument("mem", type=int, required=True, help="mismatch val")
zenith.close()
zenith = open('zenith.txt','r')
print('----------ZenithFile-----------')
print(zenith.read())
print('----------ZenithFile-----------')


print(bestAlignments)
print('hellowaorld')
