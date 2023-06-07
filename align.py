# Note: should we have binary compression at least?
# 

# Deprecated - Old .fa parser
def alignGenome(faFile):
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
    
def maxAlignmentDict(align):
    maxVal = -1 # This could bean issue TODO: -------------------------------------
    # maxAlign = []
    # maxPos = -1
    maxKey = -1
    for i in align.keys():
        # print(i[0])
        if align[i][0] > maxVal:
            # maxVal = i[0]
            # maxAlign = i[1]
            # maxPos = i[2]
            maxVal = align[i][0]
            maxKey = i
    # return [maxVal, maxAlign, maxPos]
    # print('maxindex; ' + str(i))
    return maxKey, align[maxKey]

def maxAlignment(align):
    maxVal = -1 # This could bean issue TODO: -------------------------------------
    # maxAlign = []
    # maxPos = -1
    maxIndex = -1
    for i in range(len(align)):
        # print(i[0])
        if align[i][0] > maxVal:
            # maxVal = i[0]
            # maxAlign = i[1]
            # maxPos = i[2]
            maxVal = align[i][0]
            maxIndex = i
    # return [maxVal, maxAlign, maxPos]
    # print('maxindex; ' + str(i))
    return align[maxIndex]

def stringToCigar(string):
    #Note: do we need to consider ref and read deletions
    lastLetter = string[0]
    counter = 0
    output = ''
    for i in string:
        if i != lastLetter:
            if counter != 1:
                output = output + str(counter) + lastLetter
            else:
                output = output + lastLetter
            counter = 1
            lastLetter = i
        else:
            counter += 1
        # lastLetter = i
    output = output + str(counter) + lastLetter
    return output
    

def nwAlgo(chr, refStage, id, readStage, match, mismatch, indel, gapPenalty):
    # Note: read is a tuple, [0] is the letters
    # Set up array
    # print(read[0])
    ref = '-' + str(refStage)
    read = '-' + readStage[0]
    
    # variable print block
    # print(ref)
    # print(len(ref))
    # print(read)
    # print(len(read))

    # Array prep
    array = [[(0) for h in range(len(ref))] for g in range(len(read))]
    backArray = [[(0) for h in range(len(ref))] for g in range(len(read))]

   
    for i in range(len(ref)):
        # array[0][i] = -i
        backArray[0][i] = 3
    for i in range(len(read)):
        # array[i][0] = -i
        backArray[i][0] = 3

    # print(array) # Printing initial 0'ed array
    for i in range(1,len(read)):
        for j in range(1,len(ref)):

            # Debug block
            # print('diag: ' + str(array[i-1][j-1]))
            # print('up: ' + str(array[i-1][j]))
            # print('left: ' + str(array[i][j-1]))
            # print('comparing: ' + ref[j] +' and ' + read[i]) #compare letters
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
            
            # Print block updated score per 
            # print(scores)
            # print('score: ' + str(array[i][j]))
            # exit()

    # Backtrack
    # Going backwards through reversed
    # Note refs always bigger than read, fill in '-'

    #Note: deleted suffix and prefix to save space
    alignment = []
    for i in reversed(range(1,len(read))):
        for j in reversed(range(1,len(ref))):
            if backArray[i][j] != 3:
                score = 0
                string = ['',''] # index 0 is read 1 is ref
                prefix = ['','']
                suffix = ['','']
                cigar = ''
                gapCount = 0
                isave = i
                jsave = j

                # find suffix
                if isave != len(read):
                    suffix[0] = read[isave+1:]
                    # suffix[0] = read[jsave+1:] # because len(ref) > read
                if jsave != len(ref):
                    suffix[1] = ref[jsave+1:]

                while True:
                    score += array[isave][jsave] # old score but i think wrong
                    if backArray[isave][jsave] == 0:
                        string[0] = read[isave] + string[0]
                        string[1] = ref[jsave] + string[1]
                        backArray[isave][jsave] = 3
                        isave -= 1
                        jsave -= 1
                        score += match #score maybe
                        score -= gapCount * gapPenalty
                        gapCount = 0
                        cigar = 'M' + cigar
                    elif backArray[isave][jsave] == 1:
                        string[0] = read[isave] + string[0]
                        string[1] = '-' + string[1]
                        backArray[isave][jsave] = 3
                        isave -= 1
                        score -= indel #score maybe
                        gapCount += 1
                        cigar = 'I' + cigar
                    else:
                        backArray[isave][jsave] = 3
                        string[0] = '-' + string[0]
                        string[1] = ref[jsave] + string[1]
                        jsave -= 1
                        score -= indel #score maybe
                        gapCount += 1
                        cigar = 'D' + cigar
                    if backArray[isave][jsave] == 3:
                        score -= gapCount * gapPenalty
                        gapCount = 0
                        break
                
                # find prefix
                if isave != 0:
                    prefix[0] = read[:isave]
                if jsave != 0:
                    prefix[1] = ref[:jsave]

                string[0] = prefix[0] + string[0] + suffix[0]
                string[1] = prefix[1] + string[1] + suffix[1]
                diff = len(string[1]) - len(string[0]) #suffix match
                string[0] = string[0] + '-' * diff

                #LeftPos = isave
                # print(cigar)
                cigar = stringToCigar(cigar)
                alignment.append([score, string, isave, cigar])
                # print('dubz')
    # print(array)
    best = maxAlignment(alignment)
    # print(best)
    # print(alignment)

    # print(best)
    # exit()

    #printArray block for visualizing the score and backarray
    # printArray(array, 'array')
    # printArray(backArray, 'backarray')
    # printBackArray(backArray, 'backarray')
    
    return best
    # exit()
# def align