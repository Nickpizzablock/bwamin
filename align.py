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

def printArray(array, title):
    print(title)
    for i in range(len(array)):
        print(array[i])

def printBackArray(backArray, title):
    print(title)
    for i in range(len(backArray)):
        for j in range(len(backArray[i])):
            if i == 0 or j == 0:
                backArray[i][j] = 'x'
            if backArray[i][j] == 0:
                backArray[i][j] = 'd'
            if backArray[i][j] == 1:
                backArray[i][j] = 'u'
            if backArray[i][j] == 2:
                backArray[i][j] = 'l'
            if backArray[i][j] == 3:
                backArray[i][j] = 'x'
        print(backArray[i])
    
def maxAlignment(align):
    maxVal = -1 # This could bean issue TODO: -------------------------------------
    maxAlign = []
    for i in align:
        # print(i[0])
        if i[0] > maxVal:
            maxVal = i[0]
            maxAlign = i[1]
    return [maxVal, maxAlign]

def nwAlgo(chr, refStage, id, readStage, match, mismatch, indel):
    # Note: read is a tuple, [0] is the letters
    # Set up array
    # print(read[0])
    ref = '-' + str(refStage)
    read = '-' + readStage[0]
    print(ref)
    print(len(ref))
    print(read)
    print(len(read))

    # Array prep
    array = [[(0) for h in range(len(ref))] for g in range(len(read))]
    backArray = [[(0) for h in range(len(ref))] for g in range(len(read))]

   
    for i in range(len(ref)):
        # array[0][i] = -i
        backArray[0][i] = 3
    for i in range(len(read)):
        # array[i][0] = -i
        backArray[i][0] = 3

    print(array)
    for i in range(1,len(read)):
        for j in range(1,len(ref)):
            # print('diag: ' + str(array[i-1][j-1]))
            # print('up: ' + str(array[i-1][j]))
            # print('left: ' + str(array[i][j-1]))
            print('comparing: ' + ref[j] +' and ' + read[i])
            if ref[j] == read[i]:
                diag = array[i-1][j-1] + match
            elif ref[j] == '-' or read[i] == '-':
                diag = array[i-1][j-1] - indel
            else:
                diag = array[i-1][j-1] - mismatch
            up = array[i-1][j] - indel
            left = array[i][j-1] - indel
            scores = [diag, up, left, 0]
            array[i][j] = max(scores)
            backArray[i][j] = scores.index(max(scores))
            
            print(scores)
            print('score: ' + str(array[i][j]))
            # exit()

    # Backtrack
    # Going backwards through reversed
    alignment = []
    for i in reversed(range(1,len(read))):
        for j in reversed(range(1,len(ref))):
            if backArray[i][j] != 3:
                score = 0
                string = ['',''] # index 0 is read 1 is ref
                prefix = ['','']
                suffix = ['','']
                isave = i
                jsave = j

                # find suffix
                if isave != len(read):
                    suffix[0] = read[isave+1:]
                if jsave != len(ref):
                    suffix[1] = ref[jsave+1:]

                while True:
                    score += array[isave][jsave]
                    if backArray[isave][jsave] == 0:
                        string[0] = read[isave] + string[0]
                        string[1] = ref[jsave] + string[1]
                        backArray[isave][jsave] = 3
                        isave -= 1
                        jsave -= 1
                    elif backArray[isave][jsave] == 1:
                        string[0] = read[isave] + string[0]
                        string[1] = '-' + string[1]
                        backArray[isave][jsave] = 3
                        isave -= 1
                    else:
                        backArray[isave][jsave] = 3
                        string[0] = '-' + string[0]
                        string[1] = ref[jsave] + string[1]
                        jsave -= 1
                    if backArray[isave][jsave] == 3:
                        break
                
                # find prefix
                if isave != 0:
                    prefix[0] = read[:isave]
                if jsave != 0:
                    prefix[1] = ref[:jsave]

                string[0] = prefix[0] + string[0] + suffix[0]
                string[1] = prefix[1] + string[1] + suffix[1]

                alignment.append([score, string])
                # print('dubz')
    # print(array)
    print(alignment)
    best = maxAlignment(alignment)
    print(best)
    # exit()
    printArray(array, 'array')
    printArray(backArray, 'backarray')
    # printBackArray(backArray, 'backarray')
    return best
    # exit()
# def align