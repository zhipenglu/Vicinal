# samsoftfilter.py
#
# Vicinal: accurate determination of ncRNA ends
# using chimeric reads from RNA-seq experiments
#
# Author: Zhipeng Lu
# UNC at Chapel Hill, NC 27599
# zhipengluchina@gmail.com
#
# This program filter sam files to look for partially mapped (softclipped) reads


import re, sys

if len(sys.argv) < 3 :
    print "Error: too few arguments"
    print "Usage: Python samsoftfilter.py file_local.sam file_soft.sam"
    sys.exit()
    

samfile = sys.argv[1]
outfile = sys.argv[2]

totalcount     = 0
softclipcount  = 0
softclipreads  = []
softclipmin    = 5     #minimum softclipped length
f = open(samfile, "r")
while True:
    line = f.readline()
    record = line.strip('\n').split()
    if not record: break
    elif (record[0][0] != "@") and ((record[1] == "0") or (record[1] == "16")):
        totalcount += 1
        softclips = re.findall('\d+S', record[5])
        softcliplen = []
        for i in softclips:  softcliplen.append(int(i.strip('S')))
        if (len(softcliplen) > 0) and (max(softcliplen) > softclipmin):
            softclipcount += 1
            softclipreads.append(line)
f.close()

sorted(softclipreads, key= lambda line: line.split()[2])
print "Total number of mapped reads:", totalcount
print "Total number of softclipped reads:", softclipcount

g = open(outfile, 'w')
for i in softclipreads: g.write(i)
g.close()
