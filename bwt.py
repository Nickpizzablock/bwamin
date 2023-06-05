def bwt(string):
    bwtlist = []
    string = string + '$'
    output = ''
    for i in string:
        bwtlist.append(string)
        string = string[-1] + string[:-1]
    bwtlist = sorted(bwtlist)
    for i in bwtlist:
        output = output + i[-1]
    return output

# print(bwt(input('Transform what?: ')))