#!/usr/bin/python3

import sys
import os
import re
import argparse
import numpy

######################## Get Inputs ##########################
parser = argparse.ArgumentParser(description='Make all four freq files (freq and indel/allreads and duplex).')
parser.set_defaults(s='samtools',Q=30,M=40,i='-f3 -F2304',c='')
parser.add_argument('bam', help='<file> Bam file to generate freq files.')
parser.add_argument('selector', help='<file> Selector space to consider for freq file.')
parser.add_argument('target', help='<file> Selector space to consider on-taget (no pad).')
parser.add_argument('genome', help='<file> Genome file used to generate bam.')
parser.add_argument('-s', metavar='', type=str, help='<string> Samtools version to use. Defaults to path. (' + str(parser.get_default('g')) + ')')
parser.add_argument('-Q', metavar='', type=int, help='<int> Phred score cutoff to determine high quality bases (' + str(parser.get_default('Q')) + ')')
parser.add_argument('-M', metavar='', type=int, help='<int> Minimum MAPQ score (' + str(parser.get_default('M')) + ')')
parser.add_argument('-i', metavar='', type=str, help='<str> Filter flags to select reads from bam (' + str(parser.get_default('i')) + ')')
parser.add_argument('-c', metavar='', type=str, help='<str> Chromosome to make freq file over. If empty will take all chromosomes (' + str(parser.get_default('c')) + ')')
args = parser.parse_args()
 
#Get required inputs
fl = args.bam #freq file
selector = args.selector
target = args.target
genome = args.genome
#Tailor file name
directory = os.path.dirname(os.path.abspath(fl))
fl = os.path.basename(fl)
print('Making freq files for ' + fl)
#Get optional inputs
samtools = args.s
quality = str(args.Q)
mapq = str(args.M)
insert = args.i
tgtchro = args.c

#See if it is a split selector
if re.match('[0-9]{6}',os.path.basename(selector)):
    selpad = re.match('[0-9]{6}',os.path.basename(selector))[0]
else:
    selpad = ''

#Get pad positions
if selector == target:
    dopad = False
else:
    dopad = True
    #Initialize pad dictionary
    paddict = {}
    #Parse pads only
    for line in os.popen('bedtools slop -b 1 -i {0} -g {2}.fai | bedtools subtract -a {1} -b -'.format(target,selector,genome)):
        #Split line
        tokens = line.strip().split('\t')
        #Add chromosome to dict if not already present
        if tokens[0] not in paddict: paddict[tokens[0]] = []
        #Add pad positions
        paddict[tokens[0]] += list(range(int(tokens[1]),int(tokens[2])+1))
    #Make set for efficiency
    for chro in paddict: paddict[chro] = set(paddict[chro])
    #Get edges
    edges = {}
    for line in os.popen('bedtools merge -i ' + target,'r'):
        tokens = line.strip().split('\t')
        if tokens[0] not in edges: edges[tokens[0]] = set()
        edges[tokens[0]].add(int(tokens[1]))
        edges[tokens[0]].add(int(tokens[2]))

#Get selector reference sequences for each tile
#Initialize positions file
fiposout = selpad + fl + tgtchro + '.positions'
fout = open(os.path.join(directory,fiposout),'w')
#Write out positions
for line in open(selector,'r'):
    #Break up line
    tokens = re.sub('\s+','\t',line)
    tokens = tokens.strip().split('\t')
    #Skip header if present
    if not tokens[1].isdigit() or not tokens[2].isdigit(): continue
    fout.write(tokens[0] + ':' + str(max(0, int(tokens[1])-10)) + '-' + str(int(tokens[2])+10) + '\n')
fout.close()

#Look up positions
allrefs = os.popen(samtools + ' faidx -r {0} {1}'.format(os.path.join(directory,fiposout),genome)).read().strip().split('>')[1:]
posvec = []

for ref in allrefs:    
    #Split by line    
    ref = ref.strip().split('\n')
    #print(ref)
    #Get chromosome positions 
    chr,pos = re.split(':',ref[0])    
    poss = [chr]
    poss += re.split('-',pos)
    #print(poss)
    poss[1] = int(poss[1])
    poss[2] = int(poss[2])
    #Make sequence
    seq = ''.join(ref[1:]).upper()
    poss.append(seq)
    #print(poss)
    #Append
    posvec.append(poss)

#Set current positions to first one
curpos = [0] + posvec[0]
#Remove positions file
os.system('rm {0}'.format(os.path.join(directory,fiposout))) 

#Get removed indels if available
fremoved = os.path.join(directory,re.sub('\.sorted\.bam$','.removed-indels.txt',fl))
if os.path.isfile(fremoved):
    doreminfo = True
    rem = {}
    rec = {}
    for line in open(fremoved,'r'):
        #Split line
        tokens = line.strip().split('\t')
        #See which dict to add to
        if tokens[3] == 'RM':
            dct = rem
        else:
            dct = rec
        #If chromosome not there add it
        if tokens[0] not in dct: dct[tokens[0]] = {}
        #Add positions
        dct[tokens[0]][tokens[1]] = tokens[2]
else:
    doreminfo = False

#Define positions to update curpos
def getcurpos(posvec,chro,pos,index):
    #See if just next tile (most likely)     
    if isinstance(index, int) and index+1 < len(posvec) and chro == posvec[index+1][0] and posvec[index+1][1]+5 < pos < posvec[index+1][2]-5:
        return [index+1] + posvec[index+1]
    else:
        for i in range(len(posvec)):
            if chro == posvec[i][0] and posvec[i][1]+5 < pos < posvec[i][2]-5:
                return [i] + posvec[i]        
        return ['','','','']

#Define viewcall
if 'dualindex-deduped' in fl or 'barcode-deduped' in fl:
    viewcall = ' '.join(['(',samtools,'view -H',os.path.join(directory,fl),';',samtools,'view',insert,os.path.join(directory,fl),tgtchro,'| cut -d ":" -f 1,3-)'])
    #disabling this since it is not informtive for cfRNA
    doduplex = False
    dupinsert = '--output-QNAME'
else:    
    viewcall = ' '.join([samtools,'view','-hb',insert,os.path.join(directory,fl),tgtchro,'|',samtools,'calmd','-',genome,'2> /dev/null'])
    doduplex = False
    #dupinsert = ''
    dupinsert = '--output-QNAME'

#Make pileup call
cmd = ' '.join([viewcall,'|',samtools,'mpileup --output-extra MD --output-BP',dupinsert,'-Q' + quality,'-q' + mapq, '-d10000000 -f',genome,'-l',selector, '-'])
print(cmd)

#Initialize output file names
if tgtchro == '':
    replace = "freq.paired.Q" + quality + '.txt'
    replaced = "freq.paired.duplex.Q" + quality + '.txt'
    replacei = "indels.paired.Q" + quality + '.txt'
    replacedi = "indels.paired.duplex.Q" + quality + '.txt'
else:
    replace = tgtchro + ".freq.paired.Q" + quality + '.txt'
    replaced = tgtchro + ".freq.paired.duplex.Q" + quality + '.txt'
    replacei = tgtchro + ".indels.paired.Q" + quality + '.txt'
    replacedi = tgtchro + ".indels.paired.duplex.Q" + quality + '.txt'

#Initialize output files
fout = selpad + re.sub('bam$',replace,fl)
fiout = selpad + re.sub('bam$',replacei,fl)
freq_output = open(os.path.join(directory,fout),'w')
if selpad == '': freq_output.write("CHR\tPOS\tDEPTH\tREF\tR+\tR-\tA+\tA-\tC+\tC-\tT+\tT-\tG+\tG-\tPAD\tMOTIF\n")
indel_output = open(os.path.join(directory, fiout),'w')
if selpad == '':
    if doreminfo:
        indel_output.write("CHR\tPOS\tDEPTH\tREF\tINDEL\tPLUS\tMINUS\tPAD\tMOTIF\tRECOVERED\tREMOVED\n")
    else:
        indel_output.write("CHR\tPOS\tDEPTH\tREF\tINDEL\tPLUS\tMINUS\tPAD\tMOTIF\n")
if doduplex:
    fdout = selpad + re.sub('bam$',replaced,fl)
    duplex_output = open(os.path.join(directory,fdout),'w')
    if selpad == '': duplex_output.write("CHR\tPOS\tDEPTH\tREF\tR+\tR-\tA+\tA-\tC+\tC-\tT+\tT-\tG+\tG-\tPAD\tMOTIF\n")
    fdiout = selpad + re.sub('bam$',replacedi,fl)
    duplex_indel_output = open(os.path.join(directory,fdiout),'w')
    if selpad == '':
        if doreminfo:
            duplex_indel_output.write("CHR\tPOS\tDEPTH\tREF\tINDEL\tPLUS\tMINUS\tPAD\tMOTIF\tRECOVERED\tREMOVED\n")
        else:
            duplex_indel_output.write("CHR\tPOS\tDEPTH\tREF\tINDEL\tPLUS\tMINUS\tPAD\tMOTIF\n")
            
#Call pileup
for line in os.popen(cmd):
    #Split line    
    var = line.split('\t')
    #Get variables we want
    chro = var[0]
    pos = var[1]
    posint = int(pos)
    ref = var[2]
    depth = var[3]
    calls = var[4] 
    start_pos = var[6].split(",")
    MD = var[8].split(",")
    
    ref_skips = calls.count("<") + calls.count(">")    
    depth = str(int(depth) - ref_skips)
    
    #Skip no depth
    if depth == "0": continue
    #Get pad information
    if dopad and chro in paddict and posint in paddict[chro]:
        pad = min([abs(posint-x) for x in edges[chro]]) 
    else:
        pad = 0
    #Get motif information
    if curpos[1] != chro or not isinstance(curpos[0], int) or not curpos[2]+5 < posint < curpos[3]-5: curpos = getcurpos(posvec,chro,posint,curpos[0])
    
    if isinstance(curpos[0], int) and (posint-curpos[2]-1) >= 0 and (posint-curpos[2]+4) <= len(curpos[4]):    
        motif = curpos[4][(posint-curpos[2]-1):(posint-curpos[2]+4)]
    else:
        motif = ""    
            
    #Tailor pileup calls to remove map quality and end information        
    calls = re.sub('\^.','',calls)
    calls = re.sub('\$','',calls)    
    
    #If we are doing duplex get read names
    if doduplex: names = var[7].split(',')
    #Filter and count indels
    allsizes = set(re.findall('[\+,-][0-9]+',calls))
    indels = []
    base = calls
    for sz in allsizes:        
        indels += re.findall('\\' + sz + '[ATCGNatcgn]{' + sz[1:] + '}',calls)
        base = re.sub('\\' + sz + '[ATCGNatcgn]{' + sz[1:] + '}','',base)
    #Make list of indels unique
    indels = list(set(indels))
            
    #Count frequencies of SNVs
    freq = {'.':0, ',':0, 'A':0, 'a':0, 'C':0, 'c':0, 'T':0, 't':0, 'G':0, 'g':0}
    if doduplex:
        #Initialize duplex frequencies
        dupfreq = {'.':0, ',':0, 'A':0, 'a':0, 'C':0, 'c':0, 'T':0, 't':0, 'G':0, 'g':0}
        dupdepth = 0
        for sym, nam in zip(base,names):
            #Normal
            try:
                freq[sym] += 1
            except:
                pass 
            #Duplex
            if nam.split(':')[1] == '2':                
                try:
                    dupfreq[sym] += 1
                    dupdepth += 1
                except:
                    pass
    else:
        for sym in base:
            try:
                freq[sym] += 1
            except:
                pass 
    
    #distances = {'A':[], 'C':[], 'T':[], 'G':[]}
    #if freq['.'] + freq[','] < int(depth):        
    #    for sym, spos, md in zip(base,start_pos,MD):
    #        if not sym  in "ATCGatcg":
    #            continue
    #        new_md = re.sub('\^[A-Za-z]+','MD',md)
    #        matches = re.findall("[0-9]+", new_md)
    #        mismatches = re.findall("[ATCGNatcgn]+", new_md)
    #        template_len = len(mismatches) + numpy.sum([int(x) for x in matches])                                    
    #        distances[sym.upper()].append(min(int(spos) - 1, template_len - int(spos)))
        
    #Output SNVs
    freq_output.write('\t'.join(map(str,[chro,pos,depth,ref,freq['.'],freq[','],freq['A'],freq['a'],freq['C'],freq['c'],freq['T'],freq['t'],freq['G'],freq['g'],pad,motif])) + '\n')
    if doduplex:
        duplex_output.write('\t'.join(map(str,[chro,pos,dupdepth,ref,dupfreq['.'],dupfreq[','],dupfreq['A'],dupfreq['a'],dupfreq['C'],dupfreq['c'],dupfreq['T'],dupfreq['t'],dupfreq['G'],dupfreq['g'],pad,motif])) + '\n')
    
    #Do indels
    if len(indels) > 0:        
        #Get reminfo
        if doreminfo:
            if chro in rem and pos in rem[chro]:
                remval = rem[chro][pos]
            else:
                remval = '0'
            if chro in rec and pos in rec[chro]:
                recval = rec[chro][pos]
            else:
                recval = '0'
            remstr = '\t' + recval + '\t' + remval
        else:
            remstr = ''
        #Remove n containing indels
        for indel in indels:
            if 'N' in indel or 'n' in indel:
                calls = re.sub('\\' + indel,'',calls)
        indels = [x for x in indels if 'N' not in x and 'n' not in x]
        #Get unique indels that show up
        uniqueindels = list(set([x.upper() for x in indels]))
        #Parse each unique indel
        for indel in uniqueindels:
            #Get lowrcase indel and indel size
            lindel = indel.lower()
            sz = re.match('^[\+,-][0-9]+',indel)[0]
            #Get original base and replace indel we are counting with number
            base = calls
            base = re.sub('.\\' + indel,'1',base)
            base = re.sub('.\\' + lindel,'2',base)
            #Replace all other indels
            for otherindel in uniqueindels:
                if indel == otherindel: continue
                base = re.sub('\\' + otherindel,'',base)
                base = re.sub('\\' + otherindel.lower(),'',base)
            #Debug for now            
            if len(base) - ref_skips != int(depth):
                print("indel no match")
                foutline = open('line-indel{0}.txt'.format(fl),'w')
                foutline.write(fl + '\n')
                foutline.write(line)
                foutline.close()
                os._exit(1)
            #Initialize counts
            count = [0,0]
            #Parse                        
            if doduplex:            
                #Initialize duplex frequencies
                duplexcount = [0,0]
                for sym, nam in zip(base,names):
                    #Normal
                    if sym == '1':
                        count[0] += 1
                    elif sym == '2':
                        count[1] += 1
                    else:
                        continue
                    #Duplex
                    if nam.split(':')[1] == '2':
                        if sym == '1':
                            duplexcount[0] += 1
                        elif sym == '2':
                            duplexcount[1] += 1
            else:
                for sym in base:
                    if sym == '1':
                        count[0] += 1
                    elif sym == '2':
                        count[1] += 1
                        
            indel_output.write('\t'.join([chro,pos,str(depth),ref.upper(),re.sub('[0-9]+','',indel),str(count[0]),str(count[1]),str(pad),motif]) + remstr + '\n')            
            if doduplex:            
                duplex_indel_output.write('\t'.join([chro,pos,str(dupdepth),ref.upper(),re.sub('[0-9]+','',indel),str(duplexcount[0]),str(duplexcount[1]),str(pad),motif]) + remstr + '\n')

freq_output.close()
indel_output.close()
if doduplex:
    duplex_output.close()
    duplex_indel_output.close()
