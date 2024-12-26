#!/usr/bin/python
import sys
import os
import re

#Read in input files
read1 = sys.argv[1]
read2 = sys.argv[2]

#Decompress if necessary
os.system('pigz -d ' + read1)
os.system('pigz -d ' + read2)

#Change name if files were decompressed
read1 = re.sub('\.gz$','',read1)
read2 = re.sub('\.gz$','',read2)

#Open files for reading
fin1 = open(read1,'r')
fin2 = open(read2,'r')

#Make sure headers are not too different
count = 0
force = 0
for line1 in fin1:
    line2 = next(fin2)
    
    #Update counter
    count += 1
    
    #Check length
    if len(line1) != len(line2):
        print("ERROR [force_header]: Reads at line " + count + " have different length. Exiting.")
        print(line1)
        print(line2)
        os._exit(1)

    #Check similarity
    sim = sum([i != j for i,j in zip(line1,line2)])
    if sim == 1:
        force = max(force,1)
        if force == 1:
            off1 = line1
            off2 = line2
    if sim > 1:
        force = 2
        off1 = line1
        off2 = line2
        

    for i in range(3):
        next(fin1)
        next(fin2)
    count += 3
    
    if count > 400: break

#Close files
fin1.close()
fin2.close()

#Decide what to do
name1 = os.path.basename(read1)
name2 = os.path.basename(read2)
if force == 0:
    print("Headers for files " + name1 + " and " + name2 + " seem to match. Leaving them as are.")
    os._exit(0)
elif force == 1:
    print("Headers for files " + name1 + " and " + name2 + " are off by one character. Replacing them.\nExample of mismatched headers:\n" + off1 + off2 + '\n')
elif force == 2:
    print("Headers for files " + name1 + " and " + name2 + " are off by more than one character.\nReplacing headers but user should ensure files are not off sync.\nExample of mismatched headers:\n" + off1 + off2 + '\n')

#Print action
print('Merging headers between ' + name1 + ' and ' + name2)

#Reopen files
fin1 = open(read1,'r')
fin2 = open(read2,'r')
fout2 = open(read2 + '.temp','w')

#Replace headers
for line1 in fin1:
    line2 = next(fin2)

    #Write out headers
    fout2.write(line1)
    for i in range(3):
        next(fin1)
        fout2.write(next(fin2))

#Rewrite file 2
os.system("mv " + read2 + '.temp ' + read2)