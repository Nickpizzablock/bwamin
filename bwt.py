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
        

b = input('Transform what?: ')
c = bwt(b)
e = addl2f(c)
print(c)
print(e)
