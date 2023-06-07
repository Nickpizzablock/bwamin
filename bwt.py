# Use more compression align compress to cigar
import align

def bwt(string):
    """
    Convert a string to bwa encoded string

    Parameters
    ----------
    string: str
        String to be encoded
    
    Returns
    -------
    last: str
        Resultant BWT string
    """
    bwtlist = []
    string = string + '$'
    last = ''
    for i in string:                        # Rotating string
        bwtlist.append(string)
        string = string[-1] + string[:-1]
    bwtlist = sorted(bwtlist)               # Alphabetical sort
    for i in bwtlist:                       # Take last char
        last = last + i[-1]
    return last

def bwtEnhanced(string):
    """
    Convert a string to compressed bwa encoded string

    Parameters
    ----------
    string: str
        String to be encoded
    
    Returns
    -------
    align.stringToCigar(string): str
        Resultant compressed BWT string
    """
    string = bwt(string)
    return align.stringToCigar(string)

def addl2f(last):
    """
    Creates L2F (last to first) indexes from BWT string

    Parameters
    ----------
    string: str
        BWT string to be encoded
    
    Returns
    -------
    l2f: [int]
        List of ints that correlate to corresponding first index
    """
    if '$' not in last:                     # In case bwt not made   
        last = bwt(last)
    first = sorted(last)
    l2f = [-1 for i in range(len(first))]   # Allocate array space
    for i in range(len(first)):
        firstLetter = first[i]
        for j in range(len(last)):          # xth in first is xth in last
            if firstLetter == last[j] and l2f[j] == -1:
                l2f[j] = i                  # Record corresponding index
                break
    return l2f
        
def rebuildbwt(last):
    """
    Convert bwa encoded string to a string with $ at end

    Parameters
    ----------
    string: str
        BWT string to get reversed to normal
    
    Returns
    -------
    original: str
        Original BWT string
    """
    if '$' not in last:                     # In case bwt not made  
        last = bwt(last)
    l2f = addl2f(last)
    first = sorted(last)
    original = ''
    i = 0                                   # i is now index of first
    while True:
        i = l2f.index(i)                    # last to front index
        if i != 0:
            original += first[i]            # Printing in order
        else:
            original += '$'
            break
    return original

def indexbwt(string):
    """
    From bwa encoded string to get indexes of first

    Parameters
    ----------
    string: str
        BWT string to get first indexes
    
    Returns
    -------
    index: [int]
        List of indexes that map first to the original string index
    """
    if '$' not in last:                     # In case bwt not made  
        last = bwt(last)
    l2f = addl2f(string)            
    i = 0
    place = 0
    index =  [-1 for g in range(len(string))]
    while True:                             # Following l2f in order 
        i = l2f.index(i) 
        index[i] = place
        if i == 0:                          # Back to beginning
            break
        place += 1
    return index

def find(bwtstring, query): #w is query
    """
    Finds indexes of original string from bwt string and query

    Parameters
    ----------
    bwtstring: str
        BWT string
    w: string
        query string inside find

    Returns
    -------
    locals: [int]
        0 index; List of original string indexes from exact query
    """
    l2f = addl2f(bwtstring)
    first = sorted(bwtstring)

    n = len(bwtstring)
    # combine = []
    # for i in range(len(l2f)):
    #     combine.append([bwtstring[i], l2f[i]]) # index 0 = bwt, 1
    top = 0
    bottom = n - 1
    offset = 0
    # print(bwtstring[top:bottom])
    # print(n)
    for c in query[::-1]:
        # print(c)
        # print('look: ' + bwtstring[top:bottom + 1])
        # exit()
        if c not in bwtstring[top:bottom+ 1]:
            return None #w does not exist
        i = bwtstring[top:bottom+ 1].find(c) + offset
        j = len(bwtstring[top:bottom+ 1]) - 1 - bwtstring[top:bottom+ 1][::-1].find(c) + offset
        # print(i, j)
        # print(bwtstring[top:bottom][::-1])
        # print(l2f)
        # exit()
        top = l2f[i]
        bottom = l2f[j]
        offset = top
        # print(top, bottom)
    # return all positions in original string corresponding to first[top,bottom]
    # return first[top:bottom+1] #, top, bottom

    indexes = indexbwt(bwtstring)
    locals = []
    for i in range(top,bottom+1):
        locals.append(indexes[i])
    return locals


# b = input('Transform what?: ')
# b = "BANANA"
# c = bwt(b)
# # e = addl2f(c)
# # print(c)
# # print(e)
# f = rebuildbwt(c)
# print(f)
# print(find(c, 'A'))