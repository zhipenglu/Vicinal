# Vicinal
# Vicinal: accurate determination of ncRNA ends
# using chimeric reads from RNA-seq experiments
#
# Author: Zhipeng Lu
# UNC at Chapel Hill, NC 27599
# zhipengluchina@gmail.com

import sys, os, re, argparse

#requires python2.7 and higher
if sys.version_info < (2, 7):
    print "must use python 2.7 or greater"
    sys.exit()

parser = argparse.ArgumentParser()
parser.add_argument('--len', type=int, help='read length', dest='readlen', required=True)
parser.add_argument('--rad', type=int, help='mapping radius for the softclipped fragment', dest='radius', default=100)
parser.add_argument('refdir', help='directory for reference directory')
parser.add_argument('softclipsam', help='name of the softclipped samfile')
parser.add_argument('outname', help='base name of the output wig files')

args = parser.parse_args()
#if len(args)<4:
#    print "Usage: Python", sys.argv[0], " [options] refdir softclip outname"
#    sys.exit(1)

readlen  = args.readlen
radius   = args.radius
refdir   = args.refdir
softclipsam = args.softclipsam
outname  = args.outname
chimfile = outname + "_chim.sam"
wig1file = outname + "_1.wig"
wig2file = outname + "_2.wig"



def rev_complement(seq):
    basepairs = {"A":"T", "T":"A", "G":"C", "C":"G", "Y":"R", "R":"Y", "N":"N"}
    complement = ''
    for letter in seq: complement += basepairs[letter]
    return complement[::-1]  
       

def getchrseq(chrname):
    print "Processing chimeric reads on", chrname, "..."
    chrseq = '' #store chromosome sequence as a string
    filelist = os.listdir(refdir)
    filepath = ''
    for chrfile in filelist:
        if chrfile.split(".")[0] == chrname: 
            filepath = refdir + chrfile      
            f = open(filepath, "r")
            f.readline() #skip first line
            for line in f:
                chrseq += line.strip('\n')
            f.close()
    return chrseq


def initchrposv(chrname):
    #initialize a dictionary to store matched fragments, takes several G memory for a human chromosome
    chrposv = {}
    for i in range(1, chrinfo[chrname]+1):
        chrposv[i] = 0
    return chrposv
 


#1. get the chromosome size information from the header lines
#2. get the softclipped reads and sort by chromosome
chrinfo = {} #store chr info as a dictionary {chrname:chrsize}
chimericcount = 0
softclipreads = []
infile = open(softclipsam, "r")
for line in infile:
    record = line.strip('\n').split()
    if len(record) > 3 and line[0] != "@": softclipreads.append(line)
    elif len(record) > 1 and  record[0] == "@SQ":
        chrinfo[record[1][3:]] = int(record[2][3:])
    else: continue
softclipreads.sort(key=lambda line: line.split('\t')[2])
infile.close()


#this block matches the 2nd part to the vicinity of the 1st part on the other strand
wig1 = open(wig1file, 'w')
wig2 = open(wig2file, 'w')
chim  = open(chimfile, 'w')
read1 = softclipreads[1].split()
chimcount = 0
currentchr = ''
label = softclipsam.strip("_softclipped.sam")
trackdef1 = 'track type=wiggle_0 name="' + label + '_1" visibility=full alwaysZero=on\n'
trackdef2 = 'track type=wiggle_0 name="' + label + '_2" visibility=full alwaysZero=on\n'
wig1.write(trackdef1)
wig2.write(trackdef2)


for line in softclipreads:
    record = line.strip('\n').split()
    if currentchr != record[2]:
        if currentchr != '':   #write to disk reads that are mapped to the current chromosome
            wig1.write("variableStep chrom=" + currentchr + " span=1" + '\n')
            wig2.write("variableStep chrom=" + currentchr + " span=1" + '\n')
            for pos in sorted(chrposv1.keys()):
                if chrposv1[pos] == 0: continue
                wig1.write(str(pos) + '\t' + str(chrposv1[pos]) + '\n')
            for pos in sorted(chrposv2.keys()):
                if chrposv2[pos] == 0: continue
                wig2.write(str(pos) + '\t' + str(chrposv2[pos]) + '\n')
        currentchr = record[2] #keep track of chrname
        chrseq = getchrseq(currentchr)  #get chromosome sequence
        chrposv1 = initchrposv(currentchr)#initialize dictionary to store partially mapped read
        chrposv2 = initchrposv(currentchr)#initialize dictionary to store partially mapped read
    strand     = record[1]
    matchstart = int(record[3])
    readseq    = record[9]
    refseq = chrseq[matchstart-radius: matchstart+readlen+radius]
    cigars = re.findall('\d+[MIDNSPH=X]', record[5])                #split this read to three parts
    clip1seq = readseq[0: int(cigars[0][0:-1])]                        #clip1 + match + clip2
    matchseq = readseq[int(cigars[0][0:-1]): readlen-int(cigars[-1][0:-1])] #existence is determined later          #
    clip2seq = readseq[readlen-int(cigars[-1][0:-1]):readlen]          #
    clip1seqrc = rev_complement(clip1seq)                           #
    clip2seqrc = rev_complement(clip2seq)                           #
    if strand == '0':
        if (cigars[0][-1] == 'S') and (cigars[-1][-1] == 'S') and (int(cigars[0][0:-1]) == int(cigars[0][0:-1])):
            a = [m.start() for m in re.finditer(clip1seqrc, refseq)]        
            b = [m.start() for m in re.finditer(clip2seqrc, refseq)]
            if len(a+b) == 1: #not likely both will match
                chimcount += 1
                chim.write(line)
                if len(a) == 1: #left softclip is matched
                    for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                        chrposv1[i] += 1
                    for i in range(matchstart-radius+a[0], matchstart-radius+a[0]+int(cigars[0][0:-1])):
                        chrposv2[i] -= 1
                else: #right softclip is matched
                    for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                        chrposv1[i] += 1
                    for i in range(matchstart-radius+b[0], matchstart-radius+b[0]+int(cigars[-1][0:-1])):
                        chrposv2[i] -= 1
        elif (cigars[0][-1] == 'S') and (cigars[-1][-1] == 'S') and (int(cigars[0][0:-1]) > int(cigars[0][0:-1])):
            a = [m.start() for m in re.finditer(clip1seqrc, refseq)]
            if len(a) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                    chrposv1[i] += 1
                for i in range(matchstart-radius+a[0], matchstart-radius+a[0]+int(cigars[0][0:-1])):
                    chrposv2[i] -= 1          
        elif (cigars[0][-1] == 'S') and (cigars[-1][-1] != 'S'):
            a = [m.start() for m in re.finditer(clip1seqrc, refseq)]
            if len(a) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])):
                    chrposv1[i] += 1
                for i in range(matchstart-radius+a[0], matchstart-radius+a[0]+int(cigars[0][0:-1])):
                    chrposv2[i] -= 1                  
        elif (cigars[0][-1] == 'S') and (cigars[-1][-1] == 'S') and (int(cigars[0][0:-1]) < int(cigars[0][0:-1])):
            b = [m.start() for m in re.finditer(clip2seqrc, refseq)]
            if len(b) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                    chrposv1[i] += 1
                for i in range(matchstart-radius+b[0], matchstart-radius+b[0]+int(cigars[-1][0:-1])):
                    chrposv2[i] -= 1
        elif ((cigars[0][-1] != 'S') and (cigars[-1][-1] == 'S')): 
            b = [m.start() for m in re.finditer(clip2seqrc, refseq)]
            if len(b) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[-1][0:-1])):
                    chrposv1[i] += 1
                for i in range(matchstart-radius+b[0], matchstart-radius+b[0]+int(cigars[-1][0:-1])):
                    chrposv2[i] -= 1
    else:  #strand == '16'. Note when this happens, the read sequence is reverse_complemented in the output.
        if ((cigars[0][-1] == 'S') and (cigars[-1][-1] == 'S')) and (int(cigars[0][0:-1]) == int(cigars[0][0:-1])):
            a = [m.start() for m in re.finditer(clip1seqrc, refseq)]        
            b = [m.start() for m in re.finditer(clip2seqrc, refseq)]
            if len(a+b) == 1:
                chimcount += 1
                chim.write(line)
                if len(a) == 1: #left softclip is matched
                    for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                        chrposv2[i] -= 1
                    for i in range(matchstart-radius+a[0], matchstart-radius+a[0]+int(cigars[0][0:-1])):
                        chrposv1[i] += 1
                else: #right softclip is matched
                    for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                        chrposv2[i] -= 1
                    for i in range(matchstart-radius+b[0], matchstart-radius+b[0]+int(cigars[-1][0:-1])):
                        chrposv1[i] += 1
        elif (cigars[0][-1] == 'S') and (cigars[-1][-1] == 'S') and (int(cigars[0][0:-1]) > int(cigars[0][0:-1])):
            a = [m.start() for m in re.finditer(clip1seqrc, refseq)]
            if len(a) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                    chrposv2[i] -= 1
                for i in range(matchstart-radius+a[0], matchstart-radius+a[0]+int(cigars[0][0:-1])):
                    chrposv1[i] += 1            
        elif ((cigars[0][-1] == 'S') and (cigars[-1][-1] != 'S')):
            a = [m.start() for m in re.finditer(clip1seqrc, refseq)]
            if len(a) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])):
                    chrposv2[i] -= 1
                for i in range(matchstart-radius+a[0], matchstart-radius+a[0]+int(cigars[0][0:-1])):
                    chrposv1[i] += 1                  
        elif (cigars[0][-1] == 'S') and (cigars[-1][-1] == 'S') and (int(cigars[0][0:-1]) < int(cigars[0][0:-1])):
            b = [m.start() for m in re.finditer(clip2seqrc, refseq)]
            if len(b) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[0][0:-1])-int(cigars[-1][0:-1])):
                    chrposv2[i] -= 1
                for i in range(matchstart-radius+b[0], matchstart-radius+b[0]+int(cigars[-1][0:-1])):
                    chrposv1[i] += 1            
        elif ((cigars[0][-1] != 'S') and (cigars[-1][-1] == 'S')): 
            b = [m.start() for m in re.finditer(clip2seqrc, refseq)]
            if len(b) == 1:
                chimcount += 1
                chim.write(line)
                for i in range(matchstart, matchstart+readlen-int(cigars[-1][0:-1])):
                    chrposv2[i] -= 1
                for i in range(matchstart-radius+b[0], matchstart-radius+b[0]+int(cigars[-1][0:-1])):
                    chrposv1[i] += 1

wig1.write("variableStep chrom=" + currentchr + " span=1" + '\n')
wig2.write("variableStep chrom=" + currentchr + " span=1" + '\n')
for pos in sorted(chrposv1.keys()):
    if chrposv1[pos] == 0: continue
    wig1.write(str(pos) + '\t' + str(chrposv1[pos]) + '\n')
for pos in sorted(chrposv2.keys()):
    if chrposv2[pos] == 0: continue
    wig2.write(str(pos) + '\t' + str(chrposv2[pos]) + '\n')
        
wig1.close()
wig2.close()
chim.close()            
print "Total number of mapped chimeric reads:", chimcount


