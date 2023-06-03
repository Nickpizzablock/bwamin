# Note how to detect integer or float to string?
def hd(formatVersion, sortingOrder, go):
    return '@HD\tVN:' + str(formatVersion) + '\tSO:' + sortingOrder + '\tGO:' + go + '\n'
def sq(refName, refLength, refAssembly):
    return '@SQ\tSN:' + refName + '\tLN:' + str(refLength) + '\tAS:' + refAssembly + '\n'
def rg(groupId, platform):
    return '@RG\tID:' + groupId + '\tPL:' + platform + '\n'
def pg(programId, programName, programVersion, args):
    return '@PG\ID:' + programId + '\tPN:' + programName + '\tVN:' + programVersion + '\tCL:' + args + '\n'
def co(text):
    return '@CO' + text + '\n'

#Note  will we need to format something?
def readToString(readName, flag, refName, leftPos, quality, cigar, nextReadName, nextReadPos, templateLength, seq, qual):
    out = ''
    for i in [readName, flag, refName, leftPos, quality, cigar, nextReadName, nextReadPos, templateLength, seq]:
        out = out + i +'\t'
    out = out + qual
    return out + '\n'