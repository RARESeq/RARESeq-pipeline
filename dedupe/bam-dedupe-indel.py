#PURPOSE: Dedupe barcoded reads in position-sorted indexed paired-end BAM file(s) (multithreaded)
#AUTHOR: Bogdan Luca (2021)
#Original AUTHOR2: Andre Schultz (2019)
#ORIGINAL AUTHOR: Aaron M Newman (2015) (This code was translated from perl to python. The original script was bam-dedup-params.pl)

import sys
import os
import re
from operator import itemgetter
import time
import itertools 
import argparse
from dedupe import dedupe
from dedupe import tailor_indel_seq
from dedupe import filter_discr_cigars
from dedupe import get_mate_cigar
from dedupe import updaterefseq
from dedupe import get_molecule_end

######################## Define Main functions ######################
def run(myfile,chro):
    #Define run constants
    bc_type = adapterdict[myfile] #Adapter type
    readgroup = {} #Reads in groups to be deduped
    readinfo = {} #Dict to save read header and other useful info
    singletons = {}
    prev = -1 #Initialize previously deduped position
    dedupelag = 300 #Lag for position to dedupe (in case of adjustment)
    runrefseq = [-1000000,'',''] #Vector for reference sequence [pos,chor,seq]
    allsnvs = set()
    allindels = set()
    filteredkeys = set()
    rmindelcountA = {}
    rmindelcountM = {}

    #Open files to write out
    outbam = open(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.sam',myfile)),'w')      
    outdup = open(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.duplex.sam',myfile)),'w')
    outcop = open(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.copies.txt',myfile)),'w')
    outremindel = open(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.removed-indels.txt',myfile)),'w')

    inputfiltsam = os.path.join(bamdir,re.sub('.bam$','.' + chro + '.filt.sam',myfile))    
    inputfiltbam = os.path.join(bamdir,re.sub('.bam$','.' + chro + '.filt.bam',myfile))
    
    #Define running command
    
    cmd = ' '.join([samtools,'view','-b',dedupe_flags,'-L',targetdict[myfile],os.path.join(bamdir,myfile),chro,'|',samtools,'view','-h','-','| awk \'substr($1,1,1) == \"@\" || $6 !~ \"S|H|P\"\'','|',samtools,'view','-Shb','-o',inputfiltbam,'-'])
    #cmd = ' '.join([samtools,'view',dedupe_flags,'-L',targetdict[myfile],os.path.join(bamdir,myfile),chro,'| awk \'substr($1,1,1) == \"@\" || $6 !~ \"S|H|P\"\'',">",inputfiltbam])
    os.system(cmd)    
    cmd = ' '.join([samtools,'view',inputfiltbam,'-o',inputfiltsam])
    os.system(cmd)   
    
    cmd = ' '.join([samtools,'view','-h',inputfiltbam,'|',samtools,'calmd','-e','-',genome,'2> /dev/null'])
    #cmd = ' '.join([samtools,'view','-b',dedupe_flags,'-L',targetdict[myfile],os.path.join(bamdir,myfile),chro,'|',samtools,'calmd','-e','-',genome,'2> /dev/null',"| awk '$6 !~ \"S|H|P\"'"])
        
    #Parse
    for line in os.popen(cmd,'r',-1):
        #Skip header
        if line[0] == '@': continue
        #breakup line
        var = line.strip().split('\t')
        #Filter read with too many SNVs
        seq = var[9]
        if len(seq) - seq.count('=') > snvfilter: continue
        #Convert to numeric
        for i in [1,3,7,8]: var[i] = int(var[i])
        #Check fragment length
        if abs(var[8]) < dedupe_fraglength_min: continue #or abs(var[8]) > dedupe_fraglength_max: continue
        #Update previous one
        if prev < 0: prev = var[3]

        #If position is new possible dedupe
        if prev != var[3]:            
            #Update prev
            prev = var[3]
            #Process families
            bctodel = set() #Barcodes to delete
            for key in readgroup:
                #Skip if too close
                if readgroup[key][0] > var[3] - dedupelag:
                    #print("Skipping family", readgroup[key][0])
                    continue
                #Otherwise process
                else:
                    bctodel.add(key)
                #Save for later if singleton
                if len(readgroup[key]) == 1:
                    singletons[key] = readgroup[key]
                    continue
                #Delete first element (position)
                curpos = readgroup[key][0]
                del readgroup[key][0]
                
                #Account for mismatched barcodes
                readgroup[key], readinfo[key], rmindels = filter_discr_cigars(readgroup[key],readinfo[key])

                #Save mismatched indels
                if len(readgroup[key]) == 0:
                    for indel in rmindels:
                        try:
                            rmindelcountA[indel] += 1
                        except:
                            rmindelcountA[indel] = 1
                else:
                    for indel in rmindels:
                            try:
                                rmindelcountM[indel] += 1
                            except:
                                rmindelcountM[indel] = 1

                #See if disagreeing family
                if len(readgroup[key]) == 0: continue

                #Read in more sequence if we have to (if we're within 500bp from the end read in +- 1mb)
                if runrefseq[0] + 999500 < readgroup[key][0][3] or runrefseq[0] > readgroup[key][0][3] or runrefseq[1] != readgroup[key][0][2]:
                    runrefseq = updaterefseq(runrefseq,readgroup[key][0][3],readgroup[key][0][2],samtools,genome,fai)

                #Otherwise process it                
                fsnvs, findels = process_family(readgroup[key],key,runrefseq,readinfo[key],outbam,outcop,outdup)
                #Add snvs and indels to list
                allsnvs |= fsnvs
                allindels |= findels

            #Delete processed families
            for key in bctodel:
                del readgroup[key]
                del readinfo[key]
            
            #Process singletons
            bctodel = set() #Barcodes to delete
            for key in singletons:
                #Skip if too close
                if singletons[key][0] > var[3] - (2*dedupelag):
                    continue
                else:
                    bctodel.add(key)

                #Remove position
                del singletons[key][0]

                #Read in more sequence if we have to (if we're within 500bp from the end read in +- 10kb)
                if runrefseq[0] + 999500 < singletons[key][0][3] or runrefseq[0] > singletons[key][0][3] or runrefseq[1] != singletons[key][0][2]:                    
                    runrefseq = updaterefseq(runrefseq,singletons[key][0][3],singletons[key][0][2],samtools,genome,fai)

                #Process singleton
                writeref, writequal = process_singleton(singletons[key][0],runrefseq,allindels,allsnvs)

                #Write out
                if len(writeref) > 0:
                    singletons[key][0][0] = makeheader(singletons[key][0][0],'1','1')
                    singletons[key][0][9] = writeref
                    singletons[key][0][10] = writequal
                    outbam.write('\t'.join([str(x) for x in singletons[key][0]]) + '\n')
                    outcop.write('\t'.join([singletons[key][0][5]] + key.split(':') + ['1', '1']) + '\n')
            
            #Delete processed families
            for key in bctodel:
                del singletons[key]

        #Gather read data
        var2 = var[0].split(':')
        samFlag = var[1]
        if samFlag & samFlags['suppAlignment'] == samFlags['suppAlignment']: continue
        fraglength = var[8]
        chro = var[2]
        pos1 = var[3]
        cigar = var[5]
        pos2 = var[7]
        qual = var[10]
        mate_cigar = get_mate_cigar(var)

        #Get barcode and fragment ID
        if bc_type == 'index':
            barcode = var2[-1]
            if len(var2) == 14:
                barcode = var2[-3]
            if ((samFlag & r1revV1) == r1revV1 or
                (samFlag & r1revV2) == r1revV2):
                strand = 1
            else:
                strand = 2
        else:
            if bc_type == 'flex':
                bc1 = var2[-2][1:3]
                bc2 = var2[-1][1:3]
            elif bc_type == 'tandem':
                bc1 = var2[-2]
                bc2 = var2[-1]
            if ((samFlag & r1revV1) == r1revV1 or
                (samFlag & r1revV2) == r1revV2):
                barcode = bc2 + bc1
                strand = 1
            else:
                barcode = bc1 + bc2
                strand = 2
                        
        end1, intronsig1 = get_molecule_end(pos1, cigar)
        end2, intronsig2 = get_molecule_end(pos2, mate_cigar)
        intronsig = [intronsig1, intronsig2]
        intronsig.sort()
        
        code = barcode + ':' + chro + ':' + str(min(pos1, pos2)) + ':' + str(max(end1, end2)) + ':' + ":".join(intronsig) + ":" + str("-" if int(fraglength) < 0 else "+")
        codemate = barcode + ':' + chro + ':' + str(min(pos1, pos2)) + ':' + str(max(end1, end2)) + ':' + ":".join(intronsig) + ":" + str("+" if int(fraglength) < 0 else "-")        

        #Save read id based on code
        if code not in readinfo:
            readinfo[code] = [var[0], []]
        if codemate not in readinfo:
            readinfo[codemate] = [var[0], []] #Read header, strand, read
        
        #Save reads
        try:
            readgroup[code].append(var)
        except:
            readgroup[code] = [pos1,var]
        
        #Save read info
        readinfo[code][1].append(strand)
    
    #Process remaining families
    for key in readgroup:
        #Save for later if singleton
        if len(readgroup[key]) == 1:
            try:
                singletons[key] = readgroup[key]
            except:
                singletons = {}
                singletons[key] = readgroup[key]
            continue
        #Delete first element (position)
        curpos = readgroup[key][0]
        del readgroup[key][0]
        #Account for mismatched barcodes
        readgroup[key], readinfo[key], rmindels = filter_discr_cigars(readgroup[key],readinfo[key])

        #Save mismatched indels
        if len(readgroup[key]) == 0:
            for indel in rmindels:
                try:
                    rmindelcountA[indel] += 1
                except:
                    rmindelcountA[indel] = 1
        else:
            for indel in rmindels:
                    try:
                        rmindelcountM[indel] += 1
                    except:
                        rmindelcountM[indel] = 1
        
        #See if disagreeing family
        if len(readgroup[key]) == 0: continue

        #Read in more sequence if we have to (if we're within 500bp from the end read in +- 10kb)
        if runrefseq[0] + 999500 < readgroup[key][0][3] or runrefseq[0] > readgroup[key][0][3] or runrefseq[1] != readgroup[key][0][2]:
            #runrefseq = updaterefseq(runrefseq,readgroup[key])
            runrefseq = updaterefseq(runrefseq,readgroup[key][0][3],readgroup[key][0][2],samtools,genome,fai)

        #Otherwise process it        
        fsnvs, findels = process_family(readgroup[key],key,runrefseq,readinfo[key],outbam,outcop,outdup)
        #Add snvs and indels to list
        allsnvs |= fsnvs
        allindels |= findels
    
    #Process singletons
    for key in singletons:
        #Remove position
        del singletons[key][0]
        #Read in more sequence if we have to (if we're within 500bp from the end read in +- 10kb)
        if runrefseq[0] + 999500 < singletons[key][0][3] or runrefseq[0] > singletons[key][0][3] or runrefseq[1] != singletons[key][0][2]:
            #runrefseq = updaterefseq(runrefseq,singletons[key])
            runrefseq = updaterefseq(runrefseq,singletons[key][0][3],singletons[key][0][2],samtools,genome,fai)
        #Process singleton
        writeref, writequal = process_singleton(singletons[key][0],runrefseq,allindels,allsnvs)
        #Write out
        if len(writeref) > 0:
            singletons[key][0][0] = makeheader(singletons[key][0][0],'1','1')
            singletons[key][0][9] = writeref
            singletons[key][0][10] = writequal
            outbam.write('\t'.join([str(x) for x in singletons[key][0]]) + '\n')
            outcop.write('\t'.join([singletons[key][0][5]] + key.split(':') + ['1', '1']) + '\n')
    
    #Close out files
    outbam.close()
    outdup.close()
    outcop.close()

    #Write out removed indels
    for indel in rmindelcountM:
        outremindel.write(chro + '\t' + str(indel) + '\t' + str(rmindelcountM[indel]) + '\t' + 'REC' + '\n')
    for indel in rmindelcountA:
        outremindel.write(chro + '\t' + str(indel) + '\t' + str(rmindelcountA[indel]) + '\t' + 'RM' + '\n')
    outremindel.close()

#Process a read family
def process_family(readgroup,key,runrefseq,readinfo,outbam,outcop,outdup):
    #Defupe family    
    writeseq, writequal, snvs, indels = dedupe(readgroup,runrefseq,phredict,revphredict,int(minPhred),var_ratio,min_var,consensus_required, key,samtools,genome,fai)
    #Check if they were all the same length. Should be but collisions may lead them to not be
    if writeseq == '':
        return set(), set()
    #Duplex support
    dplx = str(len(set(readinfo[1])))
    #Make read header
    header = makeheader(readinfo[0],dplx,str(len(readgroup)))
    #Find index to write out
    if dplx == "2":
        #Find index of a read 2
        ind2 = readinfo[1].index(1)
        ind1 = readinfo[1].index(2)
        #Write out the duplex
        read = readgroup[ind2]
        read[0] = header
        read[9] = writeseq
        read[10] = writequal
        outdup.write('\t'.join([str(x) for x in read]) + '\n')
    else:
        ind1 = 0
    #Write out the read
    read = readgroup[ind1]
    read[0] = header
    read[9] = writeseq
    read[10] = writequal
    outbam.write('\t'.join([str(x) for x in read]) + '\n')
    #Write out copies
    outcop.write('\t'.join([readgroup[ind1][5]] + key.split(':') + [str(len(readgroup)), dplx]) + '\n')
    #Return snvs and indels
    return snvs, indels

#Function to dedupe barcode family. Assumes format define in the run function (processed already)
def process_singleton(singleton,runrefseq,allindels,allsnvs):
    #Get reference sequence
    offset = singleton[3]-runrefseq[0]
    #Read in reference sequence and gather indels
    if 'I' in singleton[5] or 'D' in singleton[5]:
        refseq, posoffset, indels = tailor_indel_seq(runrefseq,[singleton],len(singleton[9]),samtools,genome,fai)
        #If indels are not contained, return
        if not indels.issubset(allindels):
            return '', ''
    else:
        refseq, posoffset, indels = tailor_indel_seq(runrefseq,[singleton],len(singleton[9]),samtools,genome,fai)
        
    #Initialize output
    writequal = ''
    writeref = ''
    snvs = set()
    #Parse through them
    for ref, p, ss, qq in zip(refseq,posoffset,singleton[9],singleton[10]):
        #If base is reference or low quality or part of insertion
        if ss == '=' or phredict[qq] < int(minPhred) or p == -1:
            writeref += ref
            writequal += qq
        #Otherwise we have a variant
        else:
            #Make snv
            snv = singleton[2] + ':' + str(singleton[3] + p) + ':' + ref + ':' + ss
            #See if it has family support
            if snv in allsnvs or keepsingletonvariants:
                writeref += ss
                writequal += qq
            else:
                writeref += ref
                writequal += '!'
    #Return
    return writeref, writequal

#Make header function
def makeheader(header,dplx,leng):
    #Make read header
    header = header.split(':')
    header.insert(1,dplx)
    header.insert(1,leng)
    return ':'.join(header)


def wrapup(myfile,filemode):
    #Get name insert
    if adapterdict[myfile] == 'index':
        insert = 'singleindex-deduped'
    else:
        insert = 'dualindex-deduped'
    #Get chromosomes
    chros = targetchros[targetdict[myfile]]
    #Get mode
    if filemode == 0:
        outbam = os.path.join(bamdir,re.sub('.bam$','.' + insert + '.bam',myfile))
        outbamsor = os.path.join(bamdir,re.sub('.bam$','.' + insert + '.sorted.bam',myfile))
        outinpfilt = os.path.join(bamdir,re.sub('.bam$','.filt.bam',myfile))
        outinpfiltsor = os.path.join(bamdir,re.sub('.bam$','.filt.sorted.bam',myfile))        
        outsamheader = os.path.join(bamdir,re.sub('.bam$','.' + insert + '.header.sam',myfile))
        os.system(samtools + ' view -H ' + os.path.join(bamdir,myfile) + ' > ' + outsamheader)
        insam = []
        for chro in chros:
            print (chro)
            insam.append(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.sam',myfile)))        
        insam = ' '.join(insam)        
        os.system('cat ' + outsamheader + ' ' + insam + ' | ' + samtools + ' view -b -o ' + outbam)        
        os.system(samtools + ' sort -o ' + outbamsor + ' ' + outbam)        
        os.system(samtools + ' index ' + outbamsor + ' ' + outbamsor + '.bai')
        os.system('rm ' + insam + ' ' + outbam)
        
        #write input filt bam
        insam = []
        inbam = []
        for chro in chros:
            print (chro)
            insam.append(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.filt.sam',myfile)))
            inbam.append(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.filt.bam',myfile)))
        insam = ' '.join(insam)        
        inbam = ' '.join(inbam)        
        os.system('cat ' + outsamheader + ' ' + insam + ' | ' + samtools + ' view -b -o ' + outinpfilt)        
        os.system(samtools + ' sort -o ' + outinpfiltsor + ' ' + outinpfilt)        
        os.system(samtools + ' index ' + outinpfiltsor + ' ' + outinpfiltsor + '.bai')
        os.system('rm ' + outsamheader + ' ' + insam + ' ' + inbam + ' ' + outinpfilt)
    elif filemode == 1:
        outdup = os.path.join(bamdir,re.sub('.bam$','.' + insert + '.duplex.bam',myfile))
        outdupsor = os.path.join(bamdir,re.sub('.bam$','.' + insert + '.duplex.sorted.bam',myfile))
        outdupheader = os.path.join(bamdir,re.sub('.bam$','.' + insert + '.duplex.header.sam',myfile))
        os.system(samtools + ' view -H ' + os.path.join(bamdir,myfile) + ' > ' + outdupheader)
        indup = []
        for chro in chros:
            print (chro)
            indup.append(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.duplex.sam',myfile)))
        indup = ' '.join(indup)
        os.system('cat ' + outdupheader + ' ' + indup + ' | ' + samtools + ' view -b -o ' + outdup)
        os.system(samtools + ' sort -o ' + outdupsor + ' ' + outdup)
        os.system(samtools + ' index ' + outdupsor + ' ' + outdupsor + '.bai')
        os.system('rm ' + outdupheader + ' ' + indup + ' ' + outdup)
    elif filemode == 2:
        outcop = os.path.join(bamdir,re.sub('.bam$','.' + insert + '.copies.txt',myfile))
        incop = []
        for chro in chros:
            print (chro)
            incop.append(os.path.join(bamdir,re.sub('.bam$','.' + chro + '.copies.txt',myfile)))
        incop = ' '.join(incop)
        os.system('cat ' + incop + ' > ' + outcop)
        os.system('rm ' + incop)

def forkpool(files,threads,mode):
    jobs = threads #Maximum number of jobs
    pids = {} #Process ids
    #Do deduping
    if mode == 1:
        #Loop through them
        for fl in files:
            #Define target chromosomes for sample
            targetfiles = targetchros[targetdict[fl]]
            for tgt in targetfiles:
                #If max number of jobs running wait for one to end
                while len(pids) == jobs:
                    for child in pids:
                        if child == os.waitpid(child, os.WNOHANG)[0]:
                            del pids[child]
                            break
                        time.sleep(0.5) #Small pause before checking if process has ended
                #Start new fork
                child = os.fork()
                if child > 0:
                    pids[child] = child
                else:
                    run(fl,tgt)
                    os._exit(0)
    #Wrap up
    if mode == 2:
        #Loop through them
        for fl in files:
            for filemode in range(3):
                #If max number of jobs running wait for one to end
                while len(pids) == jobs:
                    for child in pids:
                        if child == os.waitpid(child, os.WNOHANG)[0]:
                            del pids[child]
                            break
                        time.sleep(0.5) #Small pause before checking if process has ended
                #Start new fork
                child = os.fork()
                if child > 0:
                    pids[child] = child
                else:
                    wrapup(fl,filemode)
                    os._exit(0)
    #Make freq files
    if mode == 3:
        #Loop through them
        for fl in files:
            #Define target chromosomes for sample
            targetfiles = targetchros[targetdict[fl]]
            for tgt in targetfiles:
                #If max number of jobs running wait for one to end
                while len(pids) == jobs:
                    for child in pids:
                        if child == os.waitpid(child, os.WNOHANG)[0]:
                            del pids[child]
                            break
                        time.sleep(0.5) #Small pause before checking if process has ended
                #Start new fork
                child = os.fork()
                if child > 0:
                    pids[child] = child
                else:
                    makefreq(fl,tgt)
                    os._exit(0)
    #Removed families and freq files
    if mode == 4:
        #Loop through them
        for fl in files:
            for filemode in range(3,5):
                #If max number of jobs running wait for one to end
                while len(pids) == jobs:
                    for child in pids:
                        if child == os.waitpid(child, os.WNOHANG)[0]:
                            del pids[child]
                            break
                        time.sleep(0.5) #Small pause before checking if process has ended
                #Start new fork
                child = os.fork()
                if child > 0:
                    pids[child] = child
                else:
                    wrapup(fl,filemode)
                    os._exit(0)
        
    #Wait for all forks to finish
    for child in pids:
        os.waitpid(child, 0)

######################## Get Inputs ##########################
parser = argparse.ArgumentParser(description='Indel friendly barcode deduper.')
parser.set_defaults(s='samtools',g='/indexes/hg19.fa',e='add500bp',d='-F2060 -f3',t=12,m=1,f=0.5,Q=30,l=60,v=7,c=True,S=False,F=True)
parser.add_argument('bamdir', help='<dir> Directory of bams to be deduped. Will look for "sorted.bam$" files.')
parser.add_argument('sample2bc', help='<file> Sample to barcode file. Will be used to define adapter type and selector for each sample.')
parser.add_argument('-s', metavar='', type=str, help='<string> Samtools version to use. Defaults to path. (' + str(parser.get_default('s')) + ')')
parser.add_argument('-g', metavar='', type=str, help='<string> Genome with which samples were mapped. (' + str(parser.get_default('g')) + ')')
parser.add_argument('-e', metavar='', type=str, help='<string> Extendion to find padded selector. Use "." to run script without padding (' + str(parser.get_default('e')) + ')')
parser.add_argument('-d', metavar='', type=str, help='<string> Dedupe flags. Will dedupe reads matching this filter. Default includes paired and matched in proper pair, excluding unmapped read, read with mate unmapped, and supplemental alignment (' + str(parser.get_default('d')) + ')')
parser.add_argument('-t', metavar='', type=int, help='<int> Number of threads for parallelization (' + str(parser.get_default('t')) + ')')
parser.add_argument('-m', metavar='', type=int, help='<int> Minimum number of reads to call a variant (' + str(parser.get_default('m')) + ')')
parser.add_argument('-f', metavar='', type=float, help='<float> Minimum fraction of reads to call a variant (' + str(parser.get_default('f')) + ')')
parser.add_argument('-Q', metavar='', type=int, help='<int> Phred score cutoff to determine high quality bases (' + str(parser.get_default('Q')) + ')')
parser.add_argument('-l', metavar='', type=int, help='<int> Exclude fragments shorted than this length (' + str(parser.get_default('l')) + ')')
parser.add_argument('-v', metavar='', type=int, help='<int> Filter reads with more than these many non-reference bases. Will consider SNVs and insertions. (' + str(parser.get_default('v')) + ')')
parser.add_argument('-c', action='store_false' if parser.get_default('c') else 'store_true', help='Requires consensus between variants. (' + str(parser.get_default('c')) + ')')
parser.add_argument('-S', action='store_false' if parser.get_default('S') else 'store_true', help='Keep singleton variants. (' + str(parser.get_default('S')) + ')')
args = parser.parse_args()

# Define input.
bamdir = args.bamdir
sample2bc = args.sample2bc
samtools = args.s
genome = args.g
bedadd = args.e
dedupe_flags = args.d
threads = args.t
min_var = args.m
var_ratio = args.f
minPhred = str(args.Q)
dedupe_fraglength_min = args.l
snvfilter = args.v
consensus_required = args.c
keepsingletonvariants = args.S

######################## Run main routine ####################
#Read in sample to barcode file
samplebcs = open(sample2bc,"r").read().split('\n')
while samplebcs[-1] == '': samplebcs = samplebcs[:-1]
samplebcs = [x.split('\t') for x in samplebcs]

#Define phred dictionary
phredict = {}
for p in range(43):
    phredict[chr(p+33)] = p
revphredict = {}
for v in phredict:
    revphredict[phredict[v]] = v

# SAM FORMAT FLAGS
samFlags = {
    'paired': 1,
    'properPaired': 2,
    'unmapped': 4,
    'mateUnmapped': 8,
    'reverse': 16,
    'mateReverse': 32,
    'r1': 64,
    'r2': 128,
    'secondaryAlignment': 256,
    'qcFailed': 512,
    'duplicate': 1024,
    'suppAlignment': 2048
}

#Sam Flags for reverse strand
r1revV1 = samFlags['r1'] + samFlags['reverse']
r1revV2 = samFlags['r2'] + samFlags['mateReverse']

files = list(filter(re.compile('sorted.bam$').search, os.listdir(bamdir)))
files = [x for x in files if not re.search('-deduped',x)]
files = [x for x in files if not re.search('sorted.filt.sorted.bam',x)]

#Print parameters
print('No. threads:\t' + str(threads) + '\n')

#Get unique targets and match bed files
targets = []        #Different Selectors used
targetdict = {}     #Selector by sample
adapterdict = {}    #Adapter by sample
for fl in files:    
    sample = "_".join(fl.split("__")[0].split('_')[1:-1])
    print(sample, flush = True)

    found = False
    for vec in samplebcs:
        if vec[0] == sample: 
            found = True
            break
    if not found: 
        print("Error: I could not find sample '" + sample + "' in the s2b file. This error should not occur normally. If it does, there is a problem with the deduper!", flush = True)
        os._exit(1)
        
    #Get full name of bed and extended bed
    bed = os.path.basename(vec[4].strip())
    if bedadd != '.':        
        bedext = os.path.join('/selectors/',re.sub('\.bed$','.' +  bedadd + '.bed',bed))
    else: 
        bedext = os.path.join('/selectors/',bed)
    bed = os.path.join('/selectors/',bed)
    #Get unique bed files
    if bedext not in targets: targets.append(bedext)
    #Make sample dictionary
    targetdict[fl] = bedext
    print("Using selector '" + bedext + "' for deduping sample: '" + sample + "'")
    #Get adapter type
    adapterdict[fl] = vec[3].lower()
    if adapterdict[fl] == 'tanbar': adapterdict[fl] = 'tandem'

#read chromosome sizes
fai = {}
for line in open(genome + ".fai"):
    l = line.strip().split("\t")
    fai[l[0]] = int(l[1])
#print(fai)

# Get unique chromosomes
targetchros = {}    #Chromosomes in the Selector
for target in targets:
    targetfiles = []
    with open(target,'r') as f:
        for line in f:
            tmp = line.split('\t')
            if tmp[0] not in targetfiles:
                targetfiles.append(tmp[0])
    targetchros[target] = targetfiles

#Run deduping part
forkpool(files,threads,1)

#Run wrap-up for bams and copies
forkpool(files,threads,2)

#Combine removed indels files
for myfile in files:
    #Get name insert
    if adapterdict[myfile] == 'index':
        insert = 'singleindex-deduped'
    else:
        insert = 'dualindex-deduped'
    #initialize command
    cmd = []
    #Define target chromosomes for sample
    targetfiles = targetchros[targetdict[fl]]
    for tgt in targetfiles:
        cmd.append(os.path.join(bamdir,re.sub('.bam$','.' + tgt + '.removed-indels.txt',myfile)))
    #Make rm command
    rmcmd = 'rm ' + ' '.join(cmd)
    #Add last file
    cmd.append('>')
    cmd.append(os.path.join(bamdir,re.sub('.bam$','.' + insert + '.removed-indels.txt',myfile)))
    os.system('cat ' + ' '.join(cmd))
    os.system(rmcmd)
