import itertools
import re
import os
from operator import itemgetter

def dedupe(readgroup,runrefseq,phredict,revphredict,int minPhred,float var_ratio,int min_var, consensus_required, key,samtools,genome,fai):        
    #See if all reads have the same length
    tgtlen = len(readgroup[0][9])
    for read in readgroup:
        if len(read[9]) != tgtlen:
            return '', '', '', ''
    #Gather sequences and qualities
    seqs = [read[9] for read in readgroup]
    quals = [read[10] for read in readgroup]
    #Get reference sequence
    offset = readgroup[0][3]-runrefseq[0]
    #Read in reference sequence and gather indels
    
    #if 'I' in readgroup[0][5] or 'D' in readgroup[0][5]:
    #    print("Indel!!!!!!!")    
    refseq, posoffset, indels = tailor_indel_seq(runrefseq,readgroup,len(readgroup[0][9]),samtools,genome,fai)
    #else:
    #    indels = set()
    #    refseq = runrefseq[2][offset:(offset + len(readgroup[0][9]))]
    #    posoffset = list(range(len(readgroup[0][9])))
    #Save family size
    fs = len(readgroup)
    #Initialize output
    writequal = ''
    writeref = ''
    snvs = set()
    #Parse through them
    for ref, p, ss, qq in zip(refseq,posoffset,itertools.zip_longest(*seqs),itertools.zip_longest(*quals)):
        #Get quality vector and see if we have any high quality bases
        hqc = hq = lqc = lq = 0
        hqb = []
        #Parse bases
        for s, q in zip(ss,qq):
            if phredict[q] < minPhred:
                lqc += 1
                lq += phredict[q]
            else:
                hqc += 1
                hq += phredict[q]
                hqb.append(s)
        #Decide what to do
        #If we don't have any high quality bases
        if hqc == 0:
            writeref += ref
            writequal += revphredict[round(lq/lqc)]
        else:
            #Get high quality basepairs
            bps = set(hqb)
            #If they all agree
            if len(bps) == 1:
                #If they are all reference
                if hqb[0] == '=':
                    writeref += ref
                    writequal += revphredict[round(hq/hqc)]
                #If it's part of an insertion
                elif p == -1:
                    writeref += list(bps)[0]
                    writequal += revphredict[round(hq/hqc)]
                #If not enough variant support squash it
                elif hqc/fs < var_ratio or hqc < min_var:
                    writeref += ref
                    writequal += '!'
                #Otherwise it's a true variant
                else:
                    writeref += hqb[0]
                    writequal += revphredict[round(hq/hqc)]
                    snvs.add(readgroup[0][2] + ':' + str(p) + ':' + ref + ':' + hqb[0])
            else:
                #If consensus is required and we don't have it
                if consensus_required:
                    writeref += ref
                    writequal += '!'
                else:
                    for bp in bps:
                        if ss.count(bp)/fs > var_ratio:
                            writeref += bp
                            writequal += revphredict[round((lq+hq)/fs)]
                            snvs.add(readgroup[0][2] + ':' + str(p) + ':' + ref + ':' + bp)
                            break
                    else:
                        writeref += ref
                        writequal += '!'
    return writeref, writequal, snvs, indels

#Update running read of reference sequence
def updaterefseq(runrefseq,pos,chr,samtools,genome,fai):  
    #print("updating genome", chr, pos)          
    runrefseq[0] = max(1, pos - 10000)
    end = min(fai[chr], pos + 1000000)
    runrefseq[1] = chr
    #print(' '.join([samtools,'faidx',genome,chr + ':' + str(runrefseq[0]) + '-' + str(end)]),'r',1)
    runrefseq[2] = os.popen(' '.join([samtools,'faidx',genome,chr + ':' + str(runrefseq[0]) + '-' + str(end)]),'r',1).read().split('\n',1)[1].replace('\n','').upper()
    #print(runrefseq[2])
    return runrefseq

#Tailor reference sequence to match reference
def tailor_indel_seq(refseq,readgroup,leng,samtools,genome,fai):
    #Get cigar
    cigar = readgroup[0][5]
    
    #Initialize position ofset
    #Break up CIGAR
    cigarvec = re.findall('[0-9]+[I,D,M,N]',cigar)
    #Initalize variables
    newseq = ''
    posoffset = []
    indels = set()
    readbps = 0

    if readgroup[0][3] < refseq[0] or readgroup[0][3] >= refseq[0] + len(refseq[2]):
        #print("I am outside the genome chunk loaded into memory (1). Reading in new genome chunk starting at: %s:%s"%(readgroup[0][2],readgroup[0][3]))                
        refseq = updaterefseq(refseq,readgroup[0][3],readgroup[0][2],samtools,genome,fai)

    genomebps = readgroup[0][3] - refseq[0]
       
    #Parse
    sm = 0  
    for cig in cigarvec:
        numbps = int(cig[:-1])
        if cig[-1] == 'M':
            #Add to reference sequence                        
            if genomebps + numbps > len(refseq[2]):
                old_pos = refseq[0] + genomebps
                #print("I am outside the genome chunk loaded into memory(2). Reading in new genome chunk starting at: %s:%s"%(readgroup[0][2],old_pos))                
                refseq = updaterefseq(refseq,refseq[0] + genomebps,readgroup[0][2],samtools,genome,fai)
                genomebps = old_pos - refseq[0]            
            newseq += refseq[2][genomebps:(genomebps+numbps)]
            posoffset += [x + refseq[0] for x in range(genomebps,(genomebps+numbps))]
            genomebps = genomebps+numbps
            #Account for read basepairs
            readbps += numbps    
            sm += numbps        
        elif cig[-1] == 'N':
            genomebps = genomebps+numbps                        
        elif cig[-1] == 'I':
            #print("Insert", readgroup)
            #Save insertion position
            insertpos = posoffset[-1]
            #Make reference sequence accounting for insertion
            newseq += ''.join('N' for i in range(numbps))
            posoffset += [-1 for i in range(numbps)]
            #See if they all agree
            inserts = list(set([read[9][readbps:(readbps+numbps)] for read in readgroup]))
            if len(inserts) == 1:
                indels.add(':'.join(['I',readgroup[0][2],str(insertpos),inserts[0]]))
            else:
                tmpinsert = ''.join(['N' for x in range(numbps)])
                indels.add(':'.join(['I',readgroup[0][2],str(insertpos),tmpinsert]))
            #Account for read basepairs
            readbps += numbps
            sm += numbps        
        elif cig[-1] == 'D':
            #print("Delete", readgroup)
            #Remove from reference            
            genomebps = genomebps+numbps                        
            #Add to set
            indels.add(':'.join(['D',readgroup[0][2],str(posoffset[-1]),str(numbps)]))
    #Tailor size    
    #if sm != len(newseq):
    #    print(readgroup, leng)
    #    print(newseq, posoffset, indels)
        
    return newseq, posoffset, indels

def get_molecule_end(genomebps, cigar):
    cigarvec = re.findall('[0-9]+[I,D,M,N]',cigar)
    intronsig="i"
    for cig in cigarvec:
        numbps = int(cig[:-1])
        if cig[-1] == "N":                    
            intronsig += str(genomebps) + "+" + cig 
        if cig[-1] in ['M', 'N', "D"]:                    
            genomebps = genomebps+numbps

    return genomebps, intronsig

def get_indel_from_cigar(readgroup):
    indels = set()
    for read in readgroup:
        cigar = read[5]
        if 'I' not in cigar and 'D' not in cigar: continue
        curpos = read[3]
        cigarvec = re.findall('[0-9]+[I,D,M,N]',cigar)
        for cig in cigarvec:
            numbps = int(cig[:-1])
            if cig[-1] == 'M' or cig[-1] == 'N':
                #Update current position
                curpos += numbps
            elif cig[-1] == 'I' or cig[-1] == 'D':
                indels.add(curpos-1)
    return indels

def get_indel_signature(curpos, cigar, mate_curpos, mate_cigar):        
    indels = ""    
    
    if 'I' in cigar or 'D' in cigar:        
        cigarvec = re.findall('[0-9]+[I,D,M,N]',cigar)
        for cig in cigarvec:
            numbps = int(cig[:-1])
            if cig[-1] == 'M' or cig[-1] == 'N':
                #Update current position
                curpos += numbps
            elif cig[-1] == 'I' or cig[-1] == 'D':
                indels = indels + (str(curpos-1) + "gp" + cig)
    
    indels += "-"

    cigar = mate_cigar
    curpos = mate_curpos
    if 'I' in cigar or 'D' in cigar:        
        cigarvec = re.findall('[0-9]+[I,D,M,N]',cigar)
        for cig in cigarvec:
            numbps = int(cig[:-1])
            if cig[-1] == 'M' or cig[-1] == 'N':
                #Update current position
                curpos += numbps
            elif cig[-1] == 'I' or cig[-1] == 'D':
                indels = indels + (str(curpos-1) + "gp" + cig)    
    return indels

#Filter discrepant families

def filter_discr_cigars(readgroup,readinfo):
    #Initialize    
    cigarfamilies = {}
    perfcigars = []
    indelcigars = []        
    readbarcode = {}  
    indelsignatures = {}  
    
    #Parse each read
    for read in readgroup:
        chro = read[2]
        pos1 = read[3]
        cigar = read[5]
        pos2 = read[7]
        mate_cigar = get_mate_cigar(read)
        
        indelsig = get_indel_signature(pos1, cigar, pos2, mate_cigar)
        barcode = ":".join([chro, str(pos1), cigar, str(pos2), mate_cigar])        
        readbarcode[read[0]] = (barcode, indelsig)

        try:
            indelsignatures[indelsig] += 1
        except:
            indelsignatures[indelsig] = 1
        
        try:
            cigarfamilies[barcode] += 1
        except:
            cigarfamilies[barcode] = 1
    
    if len(indelsignatures) == 1: 
        if len(cigarfamilies) == 1:            
            return readgroup, readinfo, set()
        else:            
            or_x = sorted(cigarfamilies.items(), key=lambda x: (x[1],x[0]), reverse=True)
            allcount = len(readgroup)            
            for index in reversed(range(allcount)):
                k = readgroup[index][0]
                if readbarcode[k][0] != or_x[0][0]:                                        
                    del readinfo[1][index]
                    del readgroup[index]                                    
                    try:                        
                        del cigarfamilies[readbarcode[k][0]]                                    
                    except:
                        pass
            readinfo[0] = readgroup[0][0]
            return readgroup, readinfo, set()
    else:                    
        allcount = len(readgroup)
        indels = get_indel_from_cigar(readgroup) 
        for index in reversed(range(allcount)):
            k = readgroup[index][0]
            if readbarcode[k][1] != "-":                           
                del readinfo[1][index]
                del readgroup[index]                                    
                try:                        
                    del cigarfamilies[readbarcode[k][0]]
                except:
                    pass

        if len(cigarfamilies) == 0:             
            return [], [], indels
        
        if len(cigarfamilies) == 1:
            readinfo[0] = readgroup[0][0]                        
            return readgroup, readinfo, indels
        
        allcount = len(readgroup)        
        or_x = sorted(cigarfamilies.items(), key=lambda x: (x[1],x[0]), reverse=True)
        allcount = len(readgroup)                
        for index in reversed(range(allcount)):
            k = readgroup[index][0]
            if readbarcode[k][0] != or_x[0][0]:                 
                del readinfo[1][index]
                del readgroup[index]                                    
                try:                        
                    del cigarfamilies[readbarcode[k][0]]                                    
                except:
                    pass
    
        readinfo[0] = readgroup[0][0]
        return readgroup, readinfo, indels

    ##If they are all indels
    #if allcount == len(indelcigars):
    #    #If they all agree
    #    if len(list(set(indelcigars))) == 1:
    #        return readgroup, readinfo, set()
    #    #If they disagree
    #    else:
    #        print cigarfamilies
    #        indels = get_indel_from_cigar(readgroup) 
    #        #mate_indels = get_indel_from_mate_cigar(readgroup) 
    #        #print indels, mate_indels          
    #        return [], [], indels
    #
    #return [], [], set()
    #print cigarfamilies
    #
    ##If there is a combination of perfect and indels
    #indels = get_indel_from_cigar(readgroup)    
    #for index in reversed(range(allcount)):
    #    if 'I' in readgroup[index][5] or 'D' in readgroup[index][5]:
    #        del readgroup[index]
    #        del readinfo[1][index]
    #return readgroup, readinfo, indels

def get_mate_cigar(line):
    
    for el in line[11:]:
       field = el.split(':')
       if field[0] == "MC":
         return field[2]
    print("Error: BAM file doesn't seem to have the 'MC' field, containing the cigar string of the mate read! For RNA-seq deduping, this field is required for defining the barcode family!")
    print(line)
    exit(0)
 
def remove_indels(pos, cigar, pos2):
    if 'D' in cigar:        
        cigarvec = re.findall('(([0-9]+)[M]([0-9]+)[D]([0-9]+)[M])',cigar)        
        pr = False
        while len(cigarvec) > 0:                        
            delbp  = int(cigarvec[0][2])

            if pos2 + delbp < pos:               
                #print("\nCase 1", pos,  pos2)                    
                newcigarsub = str(int(cigarvec[0][1]) + int(cigarvec[0][3])) + "M"            
                pos = pos + delbp                         
                #print(cigar, pos, pos2)                    
            elif pos + delbp < pos2:               
                #print("\nCase 2", pos,  pos2)                    
                newcigarsub = str(int(cigarvec[0][1]) + int(cigarvec[0][3])) + "M"            
                #print(cigar, pos, pos2)                    
            else:
                if pos != pos2:
                    pr = True
                    print("\nCase 3", pos,  pos2)                    
                    newcigarsub = str(int(cigarvec[0][1]) + int(cigarvec[0][3]) + delbp + pos2 - pos) + "M"            
                    print(cigar, pos, pos2)                 
                else:
                    newcigarsub = str(int(cigarvec[0][1]) + int(cigarvec[0][3]) + delbp + pos2 - pos) + "M"            
                
            cigar = cigar.replace(cigarvec[0][0], newcigarsub, 1)            
            if pr: 
                print(cigar)
            cigarvec = re.findall('(([0-9]+)[M]([0-9]+)[D]([0-9]+)[M])',cigar)

    if 'I' in cigar:        
        cigarvec = re.findall('(([0-9]+)[M]([0-9]+)[I]([0-9]+)[M])',cigar)
        cnt = 0
        while len(cigarvec) > 0:                        
            newcigarsub = str(int(cigarvec[0][1]) + int(cigarvec[0][3])) + "M"            
            cigar = cigar.replace(cigarvec[0][0], newcigarsub, 1)            
            #if update_pos:                
            insbp  = int(cigarvec[0][2])
            pos = pos - insbp
            cigarvec = re.findall('(([0-9]+)[M]([0-9]+)[I]([0-9]+)[M])',cigar)            
            cnt = cnt + 1
    return pos, cigar

def get_code(barcode, chro, pos1, pos2, cigar, mate_cigar, fraglength):
    
    if fraglength > 0:
        pos1_tmp, cigar_tmp = remove_indels(pos1, cigar, pos2)    
        pos2_tmp, mate_cigar_tmp = remove_indels(pos2, mate_cigar, pos1)    
    else:
        pos1_tmp, cigar_tmp = remove_indels(pos1, cigar, pos2)    
        pos2_tmp, mate_cigar_tmp = remove_indels(pos2, mate_cigar, pos1)    
    
    return barcode + ':' + chro + ':' + str(pos1) + ':' + str(pos2) + ':' + cigar + ':' + mate_cigar + ":" + str("-" if int(fraglength) < 0 else "+")
