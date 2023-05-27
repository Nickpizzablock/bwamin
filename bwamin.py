"""
Notes:
Do mem portion. Prob not do index.

I need to have a package done by friday night


"""


import os, argparse
import align

# Parser
parser = argparse.ArgumentParser(description='minimum bwa')

# Tool options
parser.add_argument('--index', action="store_true", help='index a fasta file')
parser.add_argument('--mem', action="store_true", help='index a fasta file')

# My options
parser.add_argument('file', type=str, help="fasta input")
parser.add_argument('fastq', type=str, help="fastq input")



# Mem options (Found in bwa mem help)
# parser.add_argument('-t', type=int, default=1, help="number of threads [1]")
# parser.add_argument('-k', type=int,	default=19, help="Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19]")
# parser.add_argument('-w', type=int,	default=100, help="Band width. Essentially, gaps longer than INT will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. [100]")
# parser.add_argument('-d', type=int,	default=100, help="Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. [100]")
# parser.add_argument('-r', type=float, default=1.5, help="Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [1.5]")
# parser.add_argument('-c', type=int,	default=10000, help="Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000]")
# parser.add_argument('-P', help="In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.")
# parser.add_argument('-A', type=int,	default=1, help="Matching score. [1]")
# parser.add_argument('-B', type=int,	default=4, help="Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4]")
# parser.add_argument('-O', type=int,	default=6, help="Gap open penalty. [6]")
# parser.add_argument('-E', type=int,	default=1, help="Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1]")
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
    print('ERROR: Too many / few arguments')
    raise False

# File checking
faFile = ''
if not os.path.exists(args.file):
    print("ERROR: First file not found")
    raise False
else:
    faFile = args.file

faOut = align.alignGenome(faFile)
print(faOut)

fqOut = ''
if os.path.exists(args.fastq):
    fqFile = args.fastq
    fqOut = align.sortFqFile(fqFile)
    print(fqOut)


# print(optionsList.count(True))
# parser.add_argument("index", type=int, required=True, help="match val")
# parser.add_argument("mem", type=int, required=True, help="mismatch val")
print('hellowaorld')
