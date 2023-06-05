import align
# Use more compression align vcompress to cigar
def bwt(string):
    bwtlist = []
    string = string + '$'
    last = ''
    first = ''
    for i in string:
        bwtlist.append(string)
        string = string[-1] + string[:-1]
    bwtlist = sorted(bwtlist)
    for i in bwtlist:
        last = last + i[-1]
    return last

def bwtEnhanced(string):
    string = bwt(string)
    return align.stringToCigar(string)


def addl2f(last):
    first = sorted(last)
    l2f = [-1 for i in range(len(first))]
    for i in range(len(first)):
        firstLetter = first[i]
        for j in range(len(last)):
            if firstLetter == last[j] and l2f[j] == -1:
                l2f[j] = i
                break
    return l2f
    # for i in first:
        
def rebuildbwt(last):
    l2f = addl2f(last)
    print(l2f)
    first = sorted(last)
    original = ''
    i = 0
    while True:
        i = l2f.index(i) # i is now index of first
        if i != 0:
            original += first[i]
        else:
            original += '$'
            break
        # for j in range(len(l2f)): # finding i in l2f
        #     if l2f[j] == i: # found i in l2f
        #         original = original + first[j] # print the original string
        #         break
        # i = j
        # if i == 0:
        #     break
    return original

def find(bwtstring, w): #w is query
    # print('bwtstring: ' + bwtstring)
    # print('query: ' + w)
    
    # bwtstring = bwtstring + '$'

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
    for c in w[::-1]:
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
    return first[top:bottom+1]


# b = input('Transform what?: ')
b = "BANANA"
c = bwt(b)
# e = addl2f(c)
# print(c)
# print(e)
f = rebuildbwt(c)
print(f)
print(find(c, 'Q'))