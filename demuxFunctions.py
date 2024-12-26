import os, time, re, datetime
from operator import itemgetter
import time

def pct(enu,den,digits=1):
    if den == 0:
        return "NA"
    else:
        return ("{0:5."+str(digits)+"f}%").format(100.0 * enu / den)

def forkpool(script,cmds,threads):
    jobs = min(len(cmds),threads) #Maximum number of jobs
    pids = {} #Process ids

    #Loop through them
    for job in range(len(cmds)):
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
            os.execvp(script,cmds[job])

    for child in pids:
        os.waitpid(child, 0)

def forkpoollang(lang,cmds,threads):
    jobs = min(len(cmds),threads) #Maximum number of jobs
    pids = {} #Process ids
    #Loop through them
    for job in range(len(cmds)):
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
            os.system(lang + ' ' + ' '.join(cmds[job]))
            os._exit(0)

    for child in pids:
        os.waitpid(child, 0)

def splitFileFork(script_dir,fname,nlines,R,abs_path):
    child = os.fork()
    if child > 0:
        print('Splitting ' + R + ' files')
        return child
    else:
        script = os.path.join(script_dir,'demux','splitFile.sh')
        os.execvp(script,[script,fname,str(nlines),abs_path,R])

def demuxFork(threads,scratchfolder,script_dir,sample2bc,read1,read2,abs_path,nlines):
    #Start two processes to split the two reads
    readchildren = []
    readchildren.append(splitFileFork(script_dir,read1,nlines,'R1',abs_path))
    readchildren.append(splitFileFork(script_dir,read2,nlines,'R2',abs_path))

    #Variables
    pids = {} #Process ids
    runn = {} #Processes running
    iter = 0
    maxiter = 1e7
    rc0 = False
    rc1 = False
    #Loop through them
    while iter < maxiter:
        #If max number of jobs running wait for one to end
        if len(pids) == threads:
            print("Processes full, waiting for one to end")
        while len(pids) == threads:
            for child in pids:
                if child == os.waitpid(child, os.WNOHANG)[0]:
                    del pids[child]
                    del runn[child]
                    break
            time.sleep(0.5)

        #See if done writing
        if maxiter == 1e7:
            if rc0 == False:
                rc0 = readchildren[0] == os.waitpid(readchildren[0], os.WNOHANG)[0]
                if rc0: threads += 1
            if rc1 == False:
                rc1 = readchildren[1] == os.waitpid(readchildren[1], os.WNOHANG)[0]
                if rc1: threads += 1
            #If done writing all the files (splitting processes both ended)
            if rc0 and rc1:
                #If so find number of files written
                fls = os.listdir(scratchfolder)
                fls = [int(re.match('^R[1,2]([0-9]{6})',x).groups()[0]) for x in fls]
                maxiter = max(fls)

        #See if next file is ready, if so start new process
        fnext1 = 'R1' + str(iter+1).zfill(6)
        fnext2 = 'R2' + str(iter+1).zfill(6)
        fls = os.listdir(scratchfolder)
        #If iter+1 files exist, then iter files are done writing and we can start demuxing them
        if fnext1 in fls and fnext2 in fls:
            fname1 = 'R1' + str(iter).zfill(6) #Define R1 name
            fname2 = 'R2' + str(iter).zfill(6) #Define R2 name
            #Start new process
            child = os.fork()
            #If in the parent
            if child > 0:
                #Save process and current file
                pids[child] = child
                runn[child] = str(iter)
                #Print output
                #print(str(len(pids)) + ' processes running: ' + ', '.join(runn.values().sort()))
                #Update iterator
                iter += 1
            else:
                print('Demuxing ' + str(iter))
                script = os.path.join(script_dir,'demux','demuxFile.sh')
                os.execvp(script,['python',os.path.join(script_dir,'demux','decode_demultiplex.py'),
                    os.path.join(scratchfolder,fname1),os.path.join(scratchfolder,fname2),
                    sample2bc,scratchfolder])

    #Check again if done writing (loop can exit before last file is done writing)
    while not (rc0 and rc1):
        if rc0 == False:
            rc0 = readchildren[0] == os.waitpid(readchildren[0], os.WNOHANG)[0]
        if rc1 == False:
            rc1 = readchildren[1] == os.waitpid(readchildren[1], os.WNOHANG)[0]
    fname1 = 'R1' + str(maxiter).zfill(6) #Define R2 name
    fname2 = 'R2' + str(maxiter).zfill(6) #Define R2 name
    #Start new process
    child = os.fork()
    #If in the parent
    if child > 0:
        #Save process
        pids[child] = child
        #Update iterator
        iter += 1
    else:
        print('Demuxing last file')
        script = os.path.join(script_dir,'demux','demuxFile.sh')
        os.execvp(script,['python',os.path.join(script_dir,'demux','decode_demultiplex.py'),os.path.join(scratchfolder,fname1),os.path.join(scratchfolder,fname2),sample2bc,scratchfolder])

    #Wait for demuxing to end
    for child in pids:
        os.waitpid(child, 0)
    print('Done demuxing')

#Reads in separate stats files for demuxed files and combines data appropriately
def readStatsFiles(scratchfolder):
    #Get all files in scratch
    files = os.listdir(scratchfolder)
    files = [x for x in files if re.search('^R[0-9]*_demux_stats\.txt$',x)]
    #order them
    index = [int(x[2:8]) for x in files]
    index = sorted(range(len(index)), key=lambda k: index[k])
    files = [files[i] for i in index]
    #Read them all
    allvec = []
    vec = []
    for f in files:
        with open(os.path.join(scratchfolder,f),'r') as r:
            allvec.append(r.read().split('\n'))
    for i in range(len(allvec[0])):
        vec.append('')
    #Start time and version
    vec[0] = allvec[0][0]
    vec[1] = allvec[0][1]
    #variables to add across all files
    for i in [2,3,4,8,9,10,15,16,17,18,20,21,22,23,24,26,27,28,29,30,32,33,34,35,36,37,39,40,41,42,43,44,45,46]:
        count = 0
        for j in range(len(allvec)):
            count += int(allvec[j][i])
        vec[i] = count
    #Variables that are numerical vectors
    i6 = [0,0]
    i7 = [0,0]
    for j in range(len(allvec)):
        allvec[j][6] = allvec[j][6].split('\t')
        allvec[j][7] = allvec[j][7].split('\t')
        i6[0] += int(allvec[j][6][0])
        i6[1] += int(allvec[j][6][1])
        i7[0] += int(allvec[j][7][0])
        i7[1] += int(allvec[j][7][1])
    vec[6] = i6
    vec[7] = i7
    #Fragments per sample
    vec[5] = break3(allvec,5)
    #sample barcodes
    tmp = {}
    cols = []
    rows = []
    for j in range(len(allvec)):
        allvec[j][12] = allvec[j][12].split('\t')[:-1]
        allvec[j][13] = allvec[j][13].split('\t')[:-1]
        for s in allvec[j][12]:
            if s not in cols:
                cols.append(s)
        for s in allvec[j][13]:
            if s not in rows:
                rows.append(s)
    vec[12] = cols
    vec[13] = rows
    for c in cols:
        tmp[c] = {}
        for r in rows:
            tmp[c][r] = 0
    for j in range(len(allvec)):
        allvec[j][14] = allvec[j][14].split(';')
        for i in range(len(allvec[j][14])):
            allvec[j][14][i] = allvec[j][14][i].split('\t')[:-1]
        for c in range(len(allvec[j][12])):
            for r in range(len(allvec[j][13])):
                tmp[allvec[j][12][c]][allvec[j][13][r]] += int(allvec[j][14][c][r])
    vec[14] = tmp
    #Vectors to sort by four to make table
    for i in [11,19,25,31,38]:
        vec[i] = break4(allvec,i)
    return vec

#suport function for readStatsFiles
def break4(allvec,idx):
    tmp = {}
    for j in range(len(allvec)):
        x = allvec[j][idx].split('\t')[:-1]
        for i in range(0,len(x),4):
            try:
                tmp[x[i]]['f'] += int(x[i+1])
            except:
                tmp[x[i]] = {}
                tmp[x[i]]['k'] = x[i]
                tmp[x[i]]['f'] = 0
                tmp[x[i]]['s'] = x[i+2]
                tmp[x[i]]['t'] = x[i+3]
                tmp[x[i]]['f'] += int(x[i+1])
    tmp = sorted(tmp.values(), key=itemgetter('f'),reverse = True)
    idx = min(40,len(tmp))
    x = []
    for j in range(idx):
        x = x + [tmp[j]['k'],tmp[j]['f'],tmp[j]['s'],tmp[j]['t']]
    return x

def break3(allvec,idx):
    tmp = {}
    for j in range(len(allvec)):
        x = allvec[j][idx].split('\t')[:-1]
        for i in range(0,len(x),3):
            try:
                tmp[x[i]]['f'] += int(x[i+1])
            except:
                tmp[x[i]] = {}
                tmp[x[i]]['k'] = x[i]
                tmp[x[i]]['f'] = 0
                tmp[x[i]]['s'] = x[i+2]
                tmp[x[i]]['f'] += int(x[i+1])
    tmp = sorted(tmp.values(), key=itemgetter('f'),reverse = True)
    x = []
    for j in range(len(tmp)):
        x = x + [tmp[j]['k'],tmp[j]['f'],tmp[j]['s']]
    return x

#Writes out combined stats files
def writeStatsFiles(out,scratchfolder):
    """
    vec in order:
    0: start time (taken from first file)
    1:version (taken from first file)
    2: total fragment (added from all files) (sum_frags)
    3: demultiplexed fragments (added from all files) (sum_demux)
    4: unidentified fragments (added from all files) (unid_frags)
    5: vector, every three define (1) sample barcode, (2) read from that sample,
        and (3) sample name.
    6: vector: correct first base 1/2 [correct_first_base_i1,correct_first_base_i2]
    7: vector: correct punctuation 1/2 [correct_pm_i1,correct_pm_i2]
    8: Phix frags (phix_frags)
    9: Invalid index (invalid_index)
    10: Invalid UMI (invalid_umis)
    11: vector, present if unidentified barcodes are present
    12: vector, expected_i1_species
    13: vector, expected_i2_species
    14: matrix of values related to above
    15: total_frags
    16: expected_frags
    17: swapped_frags
    18: sum_u1
    19: tab delimited umis1. Every 4 define a umi.
    20: perfect_match1
    21: corrected1
    22: n_corrected1
    23: non_decodable1
    24: sum_u2
    25: tab delimited umis2. Every 4 define a umi.
    26: perfect_match2
    27: corrected3
    28: n_corrected2
    29: non_decodable2
    30: sum_all
    31: tab delimited I1. Evey four.
    32: perfect_match1
    33: corrected1
    34: n_corrected1
    35: other_barcode1
    36: non_decodable1
    37: phix1
    38: tab delimited I2. Every four.
    39: perfect_match2
    40: corrected2
    41: n_corrected2
    42: other_barcode2
    43: non_decodable2
    44: phix2
    45: error_corrected_frags
    46: other
    47: stop_time
    """
    vec = readStatsFiles(scratchfolder)
    out = open(os.path.join(out,'demux_stats.txt'),'w')
    out.write("Demultiplexing Statistics\n")
    out.write("-------------------------\n")
    out.write("Date: {0}\n".format(vec[0]))
    out.write("Version: {0}\n".format(vec[1]))
    out.write("\n")
    out.write("Total fragments:\t{0}\n".format(vec[2]))
    out.write("Demultiplexed:\t\t{0}\t{1}\n".format(vec[3], pct(vec[3], vec[2])))
    out.write("Unidentified:\t\t{0}\t{1}\n".format(vec[4], pct(vec[4], vec[2])))
    if vec[4] > vec[3]:
        out.write("\nWARNING: The majority of reads could not be demultiplexed. This suggests some problem with the data or sample definitions.\n")
    out.write("\n# Demultiplexed fragments\n\n")
    out.write("Total demultiplexed:\t{0}\n\n".format(vec[3]))
    out.write("## Fragments per sample\n")
    for i in range(0,len(vec[5]),3):
        out.write("{0}\t{1}\t{2}\t{3}\n".format(vec[5][i], vec[5][i+1], pct(vec[5][i+1], vec[3]),vec[5][i+2]))
    out.write("\n## Error rates\t\t    R1\t    R2\n")
    out.write("First base\t\t{0}\t{1}\n".format(pct(vec[3]-vec[6][0],vec[3]),pct(vec[3]-vec[6][1],vec[3])))
    out.write("Punctuation\t\t{0}\t{1}\n".format(pct(vec[3]-vec[7][0],vec[3]),pct(vec[3]-vec[7][1],vec[3])))
    out.write("\n# Unidentified fragments\n\n")

    out.write("Total unidentified:\t{0}\n".format(vec[4]))
    out.write("Suspected PhiX:\t\t{0}\t{1}\n".format(vec[8], pct(vec[8], vec[4])))
    out.write("Invalid Index:\t\t{0}\t{1}\t(including swapping, see below)\n".format(vec[9], pct(vec[9], vec[4])))
    out.write("Invalid UMI:\t\t{0}\t{1}\n".format(vec[10], pct(vec[10], vec[4])))
    for i in range(0, len(vec[11]), 4):
        out.write("WARNING: unexpected barcode {0} with high read count {1} ({2}+{3})\n".format(vec[11][i],vec[11][i+1],vec[11][i+2],vec[11][i+3]))
    out.write("\n# Swapping\n\n")
    out.write("## Counts\n")
    out.write("I1/I2\t")
    for i2 in vec[13]:
        out.write("\t{0}".format(i2))
    out.write("\n")
    for i1 in vec[12]:
        out.write("{0}".format(i1))
        for i2 in vec[13]:
            out.write("\t{0:8}".format(vec[14][i1][i2]))
        out.write("\n")
    out.write("\n## Percent\n")
    out.write("I1/I2\t")
    for i2 in vec[13]:
        out.write("\t{0}".format(i2))
    out.write("\n")
    for i1 in vec[12]:
        out.write("{0}".format(i1))
        for i2 in vec[13]:
            out.write("\t{0:8.3f}".format(100.0 * vec[14][i1][i2] / vec[15]))
        out.write("\n")
    out.write("\n")
    out.write("Total fragments:\t{0}\t(all expected I1 or I2, incl. invalid UMI)\n".format(vec[15]))
    out.write("Expected fragments:\t{0}\t{1}\n".format(vec[16], pct(vec[16], vec[15])))
    out.write("Swapped fragments:\t{0}\t{1}\n".format(vec[17], pct(vec[17], vec[15])))
    out.write("\n# UMIs\n\n")
    out.write("Total fragments:\t{0}\n\n".format(vec[18]))
    out.write("## UMI1 (Top40)\n")
    for idx in range(0,len(vec[19]),4):
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t\t{5}\n".format(idx/4+1, vec[19][idx], vec[19][idx+1], pct(vec[19][idx+1], vec[18], 2), vec[19][idx+2], vec[19][idx+3]))
    out.write("Perfect match:\t\t{0}\t{1}\n".format(vec[20], pct(vec[20], vec[18],2)))
    out.write("Error corrected:\t{0}\t{1}\n".format(vec[21], pct(vec[21], vec[18],2)))
    out.write("N base corrected:\t{0}\t{1}\t(included in above)\n".format(vec[22], pct(vec[22], vec[18],2)))
    out.write("Non-decodable:\t\t{0}\t{1}\n".format(vec[23], pct(vec[23], vec[18],2)))
    if vec[23] > vec[18] / 2:
        out.write("\nWARNING: The majority of UMI sequences in read 1 could not be identified. This suggests problems with the data or UMI definitions.\n")
    out.write("\n## UMI2 (Top40)\n")
    for idx in range(0,len(vec[25]),4):
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t\t{5}\n".format(idx/4+1, vec[25][idx], vec[25][idx+1], pct(vec[25][idx+1], vec[24], 2), vec[25][idx+2], vec[25][idx+3]))
    out.write("Perfect match:\t\t{0}\t{1}\n".format(vec[26], pct(vec[26], vec[24],2)))
    out.write("Error corrected:\t{0}\t{1}\n".format(vec[27], pct(vec[27], vec[24],2)))
    out.write("N base corrected:\t{0}\t{1}\t(included in above)\n".format(vec[28], pct(vec[28], vec[24],2)))
    out.write("Non-decodable:\t\t{0}\t{1}\n".format(vec[29], pct(vec[29], vec[24],2)))
    if vec[23] > vec[18] / 2:
        out.write("\nWARNING: The majority of UMI sequences in read 2 could not be identified. This suggests problems with the data or UMI definitions.\n")
    out.write("\n# Index barcodes\n\n")
    out.write("Total fragments:\t{0}\n\n".format(vec[30]))
    out.write("## I1 (Top40)\n")
    for idx in range(0,len(vec[31]),4):
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t\t{5}\n".format(idx/4+1, vec[31][idx], vec[31][idx+1], pct(vec[31][idx+1], vec[30], 2), vec[31][idx+2], vec[31][idx+3]))
    out.write("\nPerfect match:\t\t{0}\t{1}\n".format(vec[32], pct(vec[32], vec[30],2)))
    out.write("Error corrected:\t{0}\t{1}\n".format(vec[33], pct(vec[33], vec[30],2)))
    out.write("N base corrected:\t{0}\t{1}\t(included in above)\n".format(vec[34], pct(vec[34], vec[30],2)))
    out.write("Other barcodes: \t{0}\t{1}\n".format(vec[35], pct(vec[35], vec[30],2)))
    out.write("Non-decodable:\t\t{0}\t{1}\n".format(vec[36], pct(vec[36], vec[30],2)))
    out.write("Suspected PhiX:\t\t{0}\t{1}\t(included in above)\n".format(vec[37], pct(vec[37], vec[30],2)))
    out.write("\n")
    out.write("## I2 (Top40)\n")
    for idx in range(0,len(vec[38]),4):
        out.write("{0}\t{1}\t{2}\t{3}\t{4}\t\t{5}\n".format(idx/4+1, vec[38][idx], vec[38][idx+1], pct(vec[38][idx+1], vec[30], 2), vec[38][idx+2], vec[38][idx+3]))
    out.write('\n')
    out.write("\nPerfect match:\t\t{0}\t{1}\n".format(vec[39], pct(vec[39], vec[30],2)))
    out.write("Error corrected:\t{0}\t{1}\n".format(vec[40], pct(vec[40], vec[30],2)))
    out.write("N base corrected:\t{0}\t{1}\t(included in above)\n".format(vec[41], pct(vec[41], vec[30],2)))
    out.write("Other barcodes: \t{0}\t{1}\n".format(vec[42], pct(vec[42], vec[30],2)))
    out.write("Non-decodable:\t\t{0}\t{1}\n".format(vec[43], pct(vec[43], vec[30],2)))
    out.write("Suspected PhiX:\t\t{0}\t{1}\t(included in above)\n".format(vec[44], pct(vec[44], vec[30],2)))
    out.write("\n")

    out.write("# Summary\n\n")

    out.write("Total fragments:    \t{0}\n\n".format(vec[30]))
    out.write("  Demultiplexed:    \t{0}\t{1}\n".format(vec[3], pct(vec[3], vec[30])))
    out.write("    Perfect match:  \t{0}\t{1}\n".format(vec[3]-vec[45], pct(vec[3]-vec[45], vec[30])))
    out.write("    Error corrected:\t{0}\t{1}\n\n".format(vec[45], pct(vec[45], vec[30])))
    out.write("  Unidentified:     \t{0}\t{1}\n".format(vec[4], pct(vec[4], vec[30])))
    out.write("    Suspected PhiX: \t{0}\t{1}\n".format(vec[8], pct(vec[8], vec[30])))
    out.write("    Swapped:        \t{0}\t{1}\n".format(vec[17], pct(vec[17], vec[30])))
    out.write("    Invalid Index:  \t{0}\t{1}\n".format(vec[46], pct(vec[46], vec[30])))
    out.write("    Invalid UMI:    \t{0}\t{1}\n".format(vec[10], pct(vec[10], vec[30])))
    out.write("\n")
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    diff = datetime.datetime.strptime(now, "%Y-%m-%d %H:%M:%S") - datetime.datetime.strptime(vec[0], "%Y-%m-%d %H:%M:%S")
    out.write("Runtime: {0}\n\n".format(str(diff)))
    #"%Y-%m-%d %H:%M:%S"
