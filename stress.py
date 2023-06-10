import bwt

a = 5

def getBits(n):
    while(n):
        if (n & 1):
            print("1")
        else:
            print("0")
        n = n >> 1



a = 'AA'
# a = 'ATCGGCGATCGATCGTGGCTAG'
def toBits(string):
    i = 0
    ntdict = {'$': 8,
              'A': 9,
              'C': 10,
              'G': 11,
              'N': 12,
              'T': 13}
    for j in string:
        i = i + ntdict[j]
        i = i << 4
    i = i >> 4
    return i


def bitSort(intList):
    perfect = True
    temp = 0


print(toBits(a))
print(bin(toBits(a)))
print(len(bin(toBits(a)))-2)



# b = 9
# print(bin(b))
# b = b<<
exit()

# for i in range(10, 16):
#     print(i)
#     # input("press enter to continue")
#     bwt.bwt('A'* 2**i)
# print('done')

bwt.bwt('A' * 18959)
print('done')
# a = 1
# a = a << 4*10000

print(a)
# 
# 15 starts getting slow