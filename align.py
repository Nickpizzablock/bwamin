# Note: should we have binary compression at least?
# 

def alignGenome(faFile):
    print("Enering align.py")

    # Reading the faFile to store chr in dict
    genome = {}
    seq = ''
    prevLine = ''
    with open(faFile) as reader:
        line = reader.readline()
        while line != '':                       # Until EOF
            line = line.strip()
            if '>' in line:                     # Detect chr, new dict entry
                line = line.strip('>')
                if prevLine != '':             
                    genome[prevLine] = seq
                seq = ''
                genome[line] = ''
                prevLine = line
            else:
                seq = seq + line
            line = reader.readline()
        genome[prevLine] = seq                  # Fill last dict entry
    return genome

def sortFqFile(fqFile):
    reads = {}
    name = ''
    seq = ''
    score = ''
    seqScore = []
    with open(fqFile) as reader:
        line = reader.readline().strip()
        while line != '':

            name = line.strip('@')
            line = reader.readline().strip()

            seq = line
            line = reader.readline()
            line = reader.readline().strip()

            score = line
            line = reader.readline()

            # seqScore.append(list(zip(seq, score)))
            # reads[name] = seqScore
            reads[name] = (seq,score)

    return reads
            
def alignThem(ref, read):
    for i in range(len(ref[1]) - len(read[1][0]) + 1):
        if ref[1][i:len(read[1][0])+i] == read[1][0]:
            return True
    # print(ref)
    # print(read)

# def align