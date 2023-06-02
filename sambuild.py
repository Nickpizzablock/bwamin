# Note how to detect integer or float to string?
def hd(formatVersion, sortingOrder):
    return '@HD\tVN:' + formatVersion + '\tSO:' + sortingOrder
def sq(refName, refLength, refAssembly):
    return '@SQ\tSN:' + refName + '\tLN:' + refLength + '\tAS:' + refAssembly
def rg(groupId, platform):
    return '@RG\tID:' + groupId + '\tPL:' + platform
def pg(programId, programName, programVersion):
    return '@PG\ID:' + programId + '\tPN:' + programName + '\tVN:' + programVersion
def co(text):
    return '@CO' + text

#Note  will we need to format something?
def readToString(readName, flag, refName, leftPos, quality, cigar, nextReadName, nextReadPos, templateLength, seq, qual):
    out = ''
    for i in [readName, flag, refName, leftPos, quality, cigar, nextReadName, nextReadPos, templateLength, seq]:
        out = out + i +'\t'
    out = out + qual
    return out