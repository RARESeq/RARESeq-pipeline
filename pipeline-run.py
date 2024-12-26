#Import dependencies (if not standard define why dependency was added)
import os
import sys
import time
import re
import math
import demuxFunctions as demux

#Define directories we will be working with
script_dir = os.environ['pipeline_dir']
abs_path = os.path.dirname(os.path.abspath(sys.argv[1]))
timesfolder = os.path.join(abs_path,'timestamps')
if not os.path.isdir(timesfolder): 
    os.mkdir(timesfolder)
dryfile = os.path.join(timesfolder,'checks.log')
os.system("> " + dryfile)

#Define error check function
def doubleprint(msg):
    #print(msg)
    os.system("echo \"" + msg + "\" | tee -a " + dryfile)
    
def checkfile(filename):
    if not os.path.exists(filename):
        doubleprint("ERROR: " + filename + " does not exist. Exiting.")
        return False
    else:
        doubleprint("Exists: " + filename)
        return True
def checksoftware(software):
    if os.system("which " + software + " >/dev/null 2>&1") != 0:
        doubleprint("ERROR: " + software + " does not exist. Exiting.")
        return False
    else:
        doubleprint("Exists: " + software)
        return True

def make_folder(*argv):
    path = os.path.join(*argv)
    try:
        os.makedirs(path)
    except: 
        pass
    return path

#Define function to return selector name
def getbed(fl): 
    return re.search('(Sample_.*?_[a-zA-Z]*)\.',fl).group(1)

def getBedSize(bed):
    count = 0
    sz = 0
    for line in os.popen('bedtools merge -i {0}'.format(bed)):
        line = re.sub('\s+','\t',line)
        line = line.strip().split('\t')
        count += 1
        #Skip header if present
        if count == 1 and (not line[1].isdigit() or not line[2].isdigit()):
            continue
        sz += (int(line[2]) - int(line[1]) + 1)
    return sz

def breakbed(bed):
    count = 0
    sz = 0
    fout = open(os.path.join(abs_path,str(count).zfill(6) + os.path.basename(bed)),'w')
    for line in os.popen('bedtools merge -i {0}'.format(bed),'r'):
        line = re.sub('\s+','\t',line)
        vars = line.strip().split('\t')
        #Skip header if present
        if count == 0 and (not vars[1].isdigit() or not vars[2].isdigit()):
            continue
        #Update
        fout.write("\t".join(vars) + '\n')
        sz += (int(vars[2]) - int(vars[1]) + 1)
        if sz > config['bed_break']:
            sz = 0
            count += 1
            fout.close()
            fout = open(os.path.join(abs_path,str(count).zfill(6) + os.path.basename(bed)),'w')
    fout.close()

def check_and_load_sample2bc():
    #Define sample barcodes
    global checks, samplebcs, adapter, selector, samplenames, keyfolders, barcodes
    
    samplebcs = open(sample2bc,"r").read().split('\n')
    while samplebcs[-1] == '': 
        samplebcs = samplebcs[:-1]
    samplebcs = [x.split('\t') for x in samplebcs]

    #Define sample names
    samplenames = {}
    keyfolders = []
    adapter = {}
    barcodes = {}
    for vec in samplebcs:
        #Get sample type
        samplekey = 'Sample_' + vec[0] + '_' + vec[2].lower()

        if len(vec) < 4:
            doubleprint('ERROR: Not all rows in sample2barcode file have 4 columns. This can happen when other sperator than tab is used (e.g. multiple spaces).')
            checks = False

        if vec[2].lower() not in keywords:
            doubleprint('ERROR: Invalid keyword: ' + vec[2] + '. The values on the third columns of sample2barcode file are restricted to the following keywords: ' + ", ".join(keywords))
            checks = False

        if samplekey in samplenames:
            doubleprint('ERROR: Non-unique sample key: ' + adapter[samplekey] + '. The string obtained by concatenating the first and third columns of sample2barcode file need to be unique.')
            checks = False

        if vec[1] in barcodes:
            doubleprint('ERROR: Non-unique barcodes: ' + barcodes[vec[1]] + '. The values on the second columns of sample2barcode file need to be unique.')
            checks = False

        samplenames[samplekey] = vec[2].lower()
        keyfolders.append(vec[2].lower())
        barcodes[vec[1]] = samplekey 
        #Get selector and polishing files
        
        selector[samplekey] = os.path.join('/selectors/',os.path.basename(vec[4].strip()))
        checks = checks and checkfile(selector[samplekey])

        #Define adapter type
        adapter[samplekey] = vec[3].lower()
        if adapter[samplekey] not in ['flex','tandem','index']:
            doubleprint('ERROR: Adapter ' + adapter[samplekey] + ' not valid. Value must be "flex", "index", or "tandem".')
            checks = False
        else:
            doubleprint("Adapter type: " + adapter[samplekey])


        if config['bed_ext'] == '':
            selector_ext[samplekey] = selector[samplekey]
        else:
            selector_ext[samplekey] = re.sub('\.bed$','.' + config['bed_ext'] + '.bed',selector[samplekey])

        selector_freq[samplekey] = selector_ext[samplekey]
        check = checkfile(selector_freq[samplekey])
        if check:
            selector_freq_size[selector_freq[samplekey]] = getBedSize(selector_freq[samplekey])
            if selector_freq_size[selector_freq[samplekey]] > config['bed_limit']:
                if config['bed_force_freq']:
                    doubleprint('WARNING: Will generate freq files for large selector' + selector_freq[samplekey])
                else:
                    doubleprint('WARNING: Will not generate freq files for ' + selector_freq[samplekey])
        else:
            checks = False

    keyfolders = list(set(keyfolders))

#Initialize checks
checks = True
dryrun = False

#Check inputs
for fl in sys.argv:
    if fl != "dryrun":
        checks = checks and checkfile(fl)

#doubleprint
doubleprint("Read 1: " + sys.argv[1])
doubleprint("Read 2: " + sys.argv[2])
doubleprint("S2B: " + sys.argv[3])
if len(sys.argv) > 4:
    if sys.argv[4] == "dryrun":
        dryrun = True
        doubleprint("No user-specific configurations file. Will run all defaults.")        
    else:
        doubleprint("Config file: " + sys.argv[4])        
else:
    doubleprint("No user-specific configurations file. Will run all defaults.")
if len(sys.argv) > 5 and sys.argv[5] == "dryrun":
    dryrun = True

#Define input variables
read1 = sys.argv[1]
read2 = sys.argv[2]
sample2bc = sys.argv[3]
selector = {}
selector_ext = {}
selector_freq = {}
selector_freq_size = {}
adapter = {}
samplenames = {}
barcodes = {}

#See whether run is container run or not
try:
    os.environ['DOCKERRUN']
except:
    docker = False
else:
    docker = True

#Read in default configs
if docker:    
    f = open(os.path.join(script_dir,'configs/config_docker.py'),'r')
    os.system("cp " + os.path.join(script_dir,'configs/config_docker.py') + " " + os.path.join(abs_path, "configuration_used.py"))
else:    
    f = open(os.path.join(script_dir,'configs/config.py'),'r')
tmp = re.sub('\r','',f.read())
exec(tmp)

#Save original
okeys = config.keys()

#Define variables from config file
if len(sys.argv) > 4 and sys.argv[4] != "dryrun":
    f = open(sys.argv[4],'r')
    os.system("cp " + sys.argv[4] + " " + os.path.join(abs_path, "configuration_used.py"))
    tmp = re.sub('\r','',f.read())
    exec(tmp)

#Output configurations
for key in sorted(config):
    if key not in okeys:
        doubleprint("ERROR: '" + str(key) + "' is not a valid key. Exiting")
        checks = False
    doubleprint("CONFIG " + key + ": " + str(config[key]))

#Define variables
nlines = config['demux_split_lines']
avoid = config['avoid']
avoid2 = config['avoid2']

for i in avoid:
    if i not in [-1,0,1,2,3,4,5,6,7,8]:
        doubleprint("WARNING: Invalid avoid option " + str(i))

for i in avoid2:
    if i not in [-1]:
        doubleprint("WARNING: Invalid avoid2 option " + str(i))

if config['threads'] == 'auto':
    if docker:
        threads = 20
    else:
        threads = 12
else:
    threads = int(config['threads'])
doubleprint("Will run on " + str(threads) + " threads.")

selectors_dir = config['selectors_dir']
if docker:
    genome = os.path.join('/indexes',config['fa'])
    gtf = os.path.join('/indexes',config['gtf'])
    star_index = os.path.join('/indexes',config['star_index'])
    protein_coding_genes = os.path.join('/indexes',config['protein_coding_genes'])
    rnaseqc_gtf = os.path.join('/indexes',config['rnaseqc_gtf'])
else:
    genome = config['fa']
    gtf = config['gtf']
    star_index = config['star_index']
    protein_coding_genes = config['protein_coding_genes']
    rnaseqc_gtf = config['rnaseqc_gtf']

STAR = config['STAR']
RSEM_prepare_reference = config['RSEM_prepare_reference']
rnaseqc_cmd = config['rnaseqc_cmd']
sambaba_cmd = config['sambaba_cmd']
STAR_Fusion_CTAT_folder = config['STAR_Fusion_CTAT_folder']
downsampling_breaks = config['downsampling_breaks']
DNA_contamination_exon_positions = config['DNA_contamination_exon_positions']
downsampling_exon_positions = config['downsampling_exon_positions']
rnaseqc_subsample_limit = int(config['rnaseqc_subsample_limit'])
bed_ext = config['bed_ext']
dedupe_params = config['dedupe_params']
freq_params = config['freq_params']
freq_nondeduped = config['freq_nondeduped']
freq_deduped = config['freq_deduped']

#Check that files exist
checks = checks and checkfile(genome)

#Check software
checks = checks and checksoftware(STAR) and checksoftware(RSEM_prepare_reference)

#Define keywords
keywords = ['cfrna', 'cfRNA', 'tumor', 'Tumor', 'normal', 'Normal', 'pbmc', 'PBMC']

check_and_load_sample2bc()

#Define subthreads and subjobs
if config['subjobs'] == 'auto' or config['subthreads'] == 'auto':
    if docker: # If running docker version hard-code 20 threads
        if len(samplenames) in list(range(1,13)):
            #no. samples: 1, 2,3,4,5,6,7,8,9,10,11,12      
            sjvec =  [ 1, 2,3,4,5,6,4,4,5, 5, 6, 6]
            stvec =  [20,10,6,5,4,3,5,5,4, 4, 3, 3]
            subjobs = sjvec[len(samplenames)-1]
            subthreads = str(stvec[len(samplenames)-1])
        else:
            subjobs = 5
            subthreads = '4'
    else: #Otherwise hard-code 12
        if len(samplenames) in list(range(1,13)):
            #no. samples: 1,2,3,4,5,6,7,8,9,10,11,12       
            sjvec =  [ 1,2,3,4,5,6,4,4,3, 5, 6, 4]
            stvec =  [12,6,4,3,2,2,3,3,4, 2, 2, 3]
            subjobs = sjvec[len(samplenames)-1]
            subthreads = str(stvec[len(samplenames)-1])
        else:
            subjobs = 4
            subthreads = '3'
else:
    subjobs = int(config['subjobs'])
    subthreads = config['subthreads']

doubleprint('Will run ' + str(subjobs) + ' subjobs with ' + subthreads + ' subthreads for parallelization')

#Define output folders
scratchfolder = make_folder(abs_path,'scratch')
demuxfolder = make_folder(abs_path,'demultiplexed')
exp_folder = make_folder(abs_path, "exp_analysis")
genotyping_folder = make_folder(abs_path, "genotyping")
fastp_good_folder = make_folder(exp_folder, 'fastp', 'good')
fastp_qc_folder = make_folder(exp_folder, 'fastp', 'QC')
fastp_bad_folder = make_folder(exp_folder, 'fastp', 'bad')
mapping_folder = make_folder(exp_folder, 'STAR')
first_pass = os.path.join(mapping_folder, "first_pass")
second_pass = os.path.join(mapping_folder, "second_pass")
first_pass_unimapping = os.path.join(mapping_folder, "first_pass_unimapping")
second_pass_unimapping = os.path.join(mapping_folder, "second_pass_unimapping")
second_pass_genome = os.path.join(mapping_folder, "second_pass_genome")
QC_prededupe_folder = os.path.join(exp_folder, "QC_prededupe")
dedupe_folder = make_folder(exp_folder,'deduped')
barcode_dedupe_input = make_folder(dedupe_folder, "input", "barcode-deduped")
barcode_dedupe_output = make_folder(dedupe_folder, "output", "barcode-deduped")
exp_postdedupe_folder = make_folder(exp_folder,'expression', "barcode-deduped")
downsampling_pre_folder = make_folder(exp_folder,'downsampling', "barcode-deduped", "pre-ontarget-filtering")
STAR_fusion_root_folder = make_folder(genotyping_folder, 'STAR_Fusion')
#gatk_dedupfolder = make_folder(exp_folder,'deduped','gatk')

#See if we are running the job
if dryrun:
    if checks:
        doubleprint("Passed all checks. Job should run properly.")
    else:
        doubleprint("Did not pass all checks. Job will most likely not run properly. Will exit")
    os._exit(0)
else:
    if not checks:
        doubleprint("Job did not pass tests. Exiting.")
        os._exit(1)

#Initialize timestamps file
timestampfile = os.path.join(timesfolder,'timestamps.txt')
os.system("echo >> " + timestampfile)

print(avoid)

##################################Run Stage 1#######################
if -1 not in avoid:
    if 0 not in avoid:
        #####################Demux#################################
        os.system('echo "Starting demultiplexing: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
        #Start scratch folder
        if not os.path.isdir(scratchfolder): os.mkdir(scratchfolder)
        #Run demultiplexing
        demux.demuxFork(threads-2,scratchfolder,script_dir,sample2bc,read1,read2,abs_path,nlines)
        #Create demultiplexed directory if not already created
        if not os.path.isdir(demuxfolder): os.mkdir(demuxfolder)
        #Combine files
        files = os.listdir(scratchfolder) #Get names
        files = [re.sub('^R[1-2][0-9]{6}_','',x) for x in files]
        files = list(set(files))
        files = [x for x in files if x != 'demux_stats.txt']
        #Combine them
        script = os.path.join(script_dir,'demux/combineFiles.sh')
        cmds = [[script,abs_path,x] for x in files]
        demux.forkpool(script,cmds,threads)
        #Combine stats files
        print('Done combining files, combining stats files')
        demux.writeStatsFiles(demuxfolder,scratchfolder)
        for f in os.listdir(scratchfolder):
            os.unlink(os.path.join(scratchfolder,f))
        print('Done combining stats files, demuxing step done')
    
        #Get files
        files = os.listdir(demuxfolder) #Get names
        files = [x for x in files if re.search('_R1_',x)]
        files = list(set(files))
        files.sort()

        for samp in samplenames.keys():
            if not os.path.isdir(os.path.join(demuxfolder,samp)):
                make_folder(demuxfolder,samp)
        #Parse through files
        
        for f1 in files:
            #Make file names
            barcode = re.match(r".+?_S[0-9]+?_R1_(.+?).fastq", f1)
            if barcode == None:
                continue
            barcode = barcode.groups()[0]
            if not barcode in barcodes:
                doubleprint("Barcode " + barcode + "not found in the sample2barcode file. Skipping it")
                continue
            samp = barcodes[barcode]
            #Match output folder to file
            f2 = f1.replace('_R1_','_R2_')
            otpfolder = os.path.join(demuxfolder,samp)
            make_folder(otpfolder)
            #Move files
            os.system('mv ' + os.path.join(demuxfolder,f1) + ' ' + os.path.join(otpfolder,samp + "_R1.fastq"))
            os.system('mv ' + os.path.join(demuxfolder,f2) + ' ' + os.path.join(otpfolder,samp + "_R2.fastq"))
           
    #Run in forking
    if 1 not in avoid:
        os.system('echo "Starting fastp processing: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
        #Get files
        folders = os.listdir(demuxfolder) #Get names
        folders = [x for x in folders if x != 'barcode-deduped' and x != 'demux_stats.txt']
        folders = list(set(folders))
        folders.sort()
        #Initialize
        cmds = []
        script = os.path.join(script_dir,'fastp/fastp.sh')
        #Parse through files
        for fol in folders:
            if not re.search("^Sample_", fol):
                continue
            #Get full folder
            otpfolder = os.path.join(demuxfolder,fol)
            #Get files
            files = os.listdir(otpfolder)
            f1 = next(filter(re.compile('_R1').search,files))
            f2 = next(filter(re.compile('_R2').search,files))
            #Make file names
            f = f1.replace('_R1','')
            f = f.replace('.fastq','')
            fb = f1.replace('.fastq','.bad.fq')
            fb = fb.replace('_R1','')

            #Run with single output as bad
            cmds.append([script,os.path.join(otpfolder,f1),os.path.join(otpfolder,f2),os.path.join(fastp_qc_folder,f),
                config['afterp_overlap_len_require'],\
                config['afterp_overlap_diff_limit'],config['afterp_length_required'],config['afterp_qualified_quality_phred'],config['afterp_unqualified_percent_limit'],\
                config['afterp_trim_r1_front'], config['afterp_trim_r1_tail'], config['afterp_trim_r2_front'], config['afterp_trim_r2_tail'],
                os.path.join(fastp_bad_folder,fb), os.path.join(fastp_good_folder,fol), subthreads])
        #Run forkpool
        demux.forkpool(script,cmds,subjobs)
        os.system("Rscript %s/fastp/fastp_stats.R %s"%(script_dir, fastp_qc_folder)) 

    if 2 not in avoid:
        ##############################Mapping############################
        os.system('echo "Starting mapping: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
       
        #Initialize times file
        mapstampfile = os.path.join(timesfolder,'map_timestamps.txt')
        open(mapstampfile, 'w').close() #Clear this file
        
        #Generate STAR and RSEM indices if they don't exist
        if not os.path.exists(os.path.join(star_index, "SAindex")):
            doubleprint("The STAR index doesn't exist. Generating it now. If the job doesn't have access to >64GB RAM and >30GB space on disk, it might fail. This process might take up to 2-3 hours.")
            os.system("sh %s/mapping/mapping_generate_STAR_index.sh %s %s %s %s %s %s"%(script_dir, STAR, RSEM_prepare_reference, star_index, genome, gtf, threads))
        
         #Parse through samples
        for samp in samplenames:
            #Get folder
            otpfolder = os.path.join(fastp_good_folder,samp)
            os.system('echo "First-pass mapping file ' + samp + ": " + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
            os.system("sh %s/mapping/mapping_STAR_1pass.sh %s %s %s %s %s"%(script_dir, STAR, otpfolder, star_index, first_pass, threads))
        
        os.system("sh %s/mapping/mapping_generate_STAR_index_pass2.sh %s %s %s %s %s %s %s"%(script_dir, STAR, otpfolder, first_pass, second_pass_genome, gtf, genome, threads))
        
        #Parse through samples
        for samp in samplenames:
            #Get folder
            otpfolder = os.path.join(fastp_good_folder,samp)
            os.system('echo "Second-pass mapping file ' + samp + ": " + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
            os.system("sh %s/mapping/mapping_STAR_2pass.sh %s %s %s %s %s"%(script_dir, STAR, otpfolder, second_pass, second_pass_genome, threads))
        
        cmds = []
        script = os.path.join(script_dir,'mapping/uniquely_mapping_reads.sh')
        
        os.system('echo "Starting extracting unimapping reads: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
        
        for samp in samplenames:
            os.system('echo "Extracting uniquely mapping reads for file ' + samp + ": " + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
            cmds.append([script, os.path.join(second_pass, samp + ".Aligned.sortedByCoord.out.bam"), second_pass_unimapping])
        demux.forkpoollang("sh",cmds,threads)  
        os.system("Rscript %s/mapping/STAR_tabular_logs.R %s"%(script_dir, second_pass))  
        os.system("Rscript %s/mapping/STAR_plot_tabular_logs.R %s"%(script_dir, os.path.dirname(second_pass)))
        
    if 3 not in avoid:
        os.system('echo "Starting pre-dedupe QC analysis: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)  
        
        script = os.path.join(script_dir,'bashes/downsample_bam_sambaba.sh')
        cmds = []
        n_reads_d = {}
        for samp in samplenames:    
            n_reads = rnaseqc_subsample_limit
            n_reads_d[samp] = n_reads
            cmds.append([script, 
                sambaba_cmd,
                os.path.join(second_pass_unimapping, samp + ".bam"), 
                str(n_reads),
                str(subthreads),
                os.path.join(QC_prededupe_folder, samp)])
        demux.forkpoollang("sh", cmds, subjobs) 

        script = os.path.join(script_dir,'QC_prededupe/RNASeQC.sh')
        cmds = []
        for samp in samplenames:
            n_reads = n_reads_d[samp]
            cmds.append([script, 
                rnaseqc_cmd, 
                rnaseqc_gtf, 
                os.path.join(QC_prededupe_folder, samp, samp + "__" + str(n_reads) + ".sorted.bam"), 
                os.path.join(QC_prededupe_folder, samp)
                ])
        demux.forkpoollang("sh", cmds, threads) 

        script = os.path.join(script_dir,'QC_prededupe/RNASeQC_stats.R')
        os.system("Rscript " + script + " " + QC_prededupe_folder)

        script = os.path.join(script_dir,'QC_prededupe/bam_input_molecule_length.py')
        cmds = []
        for samp in samplenames:
            n_reads = n_reads_d[samp]
            cmds.append([script, 
                os.path.join(QC_prededupe_folder, samp, samp + "__" + str(n_reads) + ".sorted.bam"),
                os.path.join(QC_prededupe_folder, samp)
                ])
        demux.forkpoollang("python", cmds, threads) 

        script = os.path.join(script_dir,'QC_prededupe/bam_input_molecule_length.R')
        os.system("Rscript " + script + " " + QC_prededupe_folder)

        script = os.path.join(script_dir,'QC_prededupe/DNA_contamination.sh')
        cmds = []
        for samp in samplenames:
            n_reads = n_reads_d[samp]
            cmds.append([script, 
                os.path.join(QC_prededupe_folder, samp, samp + "__" + str(n_reads) + ".sorted.bam"),
                DNA_contamination_exon_positions, 
                os.path.join(QC_prededupe_folder, samp)
                ])
        demux.forkpoollang("sh", cmds, threads) 

        script = os.path.join(script_dir,'QC_prededupe/DNA_contamination_plot.R')
        os.system("Rscript " + script + " " + QC_prededupe_folder)

    if 4 not in avoid:
        os.system('echo "Starting to prepare input for deduping: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)     
        cmds = []
        script = os.path.join(script_dir,'dedupe_input/dedupe_prepare_input.sh')
        for samp in samplenames:
            #Get folder 
            cmds.append([script, os.path.join(second_pass_unimapping, samp + ".bam"), genome + ".fai", barcode_dedupe_input])        
        demux.forkpoollang("sh", cmds, threads) 
         
        os.system('echo "Starting deduping: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)     
        cmds = []
        for samp in samplenames:
            #Get folder 
            make_folder(os.path.join(barcode_dedupe_output, samp))
            os.system('echo "Barcode deduping ' + samp + ": " + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
            os.system('ln -sf ' + os.path.join(barcode_dedupe_input, samp + '.sorted.bam') + ' ' + os.path.join(barcode_dedupe_output, samp, samp + '.sorted.bam'))
            os.system('ln -sf ' + os.path.join(barcode_dedupe_input, samp + '.sorted.bam.bai') + ' ' + os.path.join(barcode_dedupe_output, samp, samp + '.sorted.bam.bai'))
            cmds.append(["%s/dedupe/bam-dedupe-indel.py"%(script_dir), "-g", genome, "-e", bed_ext, "-t", subthreads, dedupe_params, os.path.join(barcode_dedupe_output, samp), sample2bc])
        demux.forkpoollang("python", cmds, subjobs) 

        cmds = []
        script = os.path.join(script_dir,'bashes/bam_get_number_reads.sh')        
        for samp in samplenames:
            cmds.append([script, 
                os.path.join(barcode_dedupe_output, samp, samp + ".sorted.bam")])
            cmds.append([script, 
                os.path.join(barcode_dedupe_output, samp, samp + ".sorted.filt.sorted.bam")])
        demux.forkpoollang("sh", cmds, threads)

        
    if 5 not in avoid:
        os.system('echo "Starting to quantify expression: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)     
       
        script = os.path.join(script_dir,'bashes/recover_deduped_read_ids.sh') 
        cmds = []
        for samp in samplenames:
            cmds.append([script, os.path.join(barcode_dedupe_output, samp, samp + '.sorted.dualindex-deduped.sorted.bam')])
        demux.forkpoollang("sh",cmds,threads) 
        
        script = os.path.join(script_dir,'bashes/downsample_bam_defined_reads.sh') 
        cmds = []
        for samp in samplenames:
            cmds.append([script, 
                os.path.join(second_pass, samp + ".Aligned.toTranscriptome.out.bam"), 
                os.path.join(barcode_dedupe_output, samp, samp + ".sorted.dualindex-deduped.sorted.bam.read_ids"), 
                exp_postdedupe_folder])
        demux.forkpoollang("sh",cmds,threads)  
        
        script = os.path.join(script_dir,'expression_quantification/expression_quantification_RSEM.sh') 
        cmds = []
        cmds2 = []
        for samp in samplenames:
            os.system('echo "Quantifying expression for sample ' + samp + ': ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)                 
            cmds.append([script, 
                os.path.join(exp_postdedupe_folder, samp + ".Aligned.toTranscriptome.out.bam"),
                os.path.join(star_index, "RSEM"),
                exp_postdedupe_folder, str(threads)])
            cmds2.append(["%s/expression_quantification/expression_add_gene_symbol.R %s %s"%(script_dir, os.path.join(exp_postdedupe_folder, samp + ".genes.results"), protein_coding_genes)])
        demux.forkpoollang("sh",cmds,threads)  
        demux.forkpoollang("Rscript",cmds2,threads)  
        os.system("Rscript %s/expression_quantification/expression_tabular.R %s %s"%(script_dir, exp_postdedupe_folder, exp_postdedupe_folder))
        
    if 6 not in avoid:
        os.system('echo "Starting the downsampling analysis: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)     

        script = os.path.join(script_dir,'bashes/downsample_bam_sambaba.sh')
        cmds = []
        n_reads_d = {}
        for samp in samplenames:    
            for fract in downsampling_breaks:
                n_reads = int(fract)
                available_reads = int(open(os.path.join(barcode_dedupe_output, samp, samp + ".sorted.nr_reads"), "r").read())
                n_reads_d[samp] = available_reads
                if n_reads < n_reads_d[samp]:
                    cmds.append([script, 
                        sambaba_cmd,
                        os.path.join(barcode_dedupe_output, samp, samp + ".sorted.bam"), 
                        str(n_reads),
                        str(subthreads),
                        os.path.join(downsampling_pre_folder, samp + "__" + str(n_reads))])
        demux.forkpoollang("sh", cmds, subjobs) 

        cmds = []
        for samp in samplenames:
            for fract in downsampling_breaks:
                n_reads = int(fract)
                if n_reads < n_reads_d[samp]:
                    os.system('echo "Barcode deduping downsamples file ' + samp + " - " + str(n_reads) + ": " + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)                    
                    cmds.append(["%s/dedupe/bam-dedupe-indel.py -g %s -e %s -t %s %s %s %s"%(script_dir, genome, bed_ext, subthreads, dedupe_params, os.path.join(downsampling_pre_folder, samp + "__" + str(n_reads)), sample2bc)])        
        demux.forkpoollang("python", cmds, subjobs) 
        
        cmds = []
        script = os.path.join(script_dir,'bashes/bam_coverage_exons.sh')
        for samp in samplenames:
            for fract in downsampling_breaks:
                n_reads = int(fract)
                if n_reads < n_reads_d[samp]:
                    cmds.append([script, 
                        os.path.join(downsampling_pre_folder, samp + "__" + str(n_reads), samp + "__" + str(n_reads) + ".sorted.dualindex-deduped.sorted.bam"),
                        downsampling_exon_positions,
                        os.path.join(downsampling_pre_folder, samp + "__" + str(n_reads))])
            cmds.append([script, os.path.join(barcode_dedupe_output, samp, samp + '.sorted.dualindex-deduped.sorted.bam'),
                downsampling_exon_positions,
                os.path.join(barcode_dedupe_output, samp)
                ])
        demux.forkpoollang("sh", cmds, threads)  
        
        cmds = []
        script = os.path.join(script_dir,'bashes/bam_get_number_reads.sh')
        for samp in samplenames:
            for fract in downsampling_breaks:
                n_reads = int(fract)
                if n_reads < n_reads_d[samp]:
                    cmds.append([script, 
                        os.path.join(downsampling_pre_folder, samp + "__" + str(n_reads), samp + "__" + str(n_reads) + ".sorted.bam")])
                    cmds.append([script, 
                        os.path.join(downsampling_pre_folder, samp + "__" + str(n_reads), samp + "__" + str(n_reads) + ".sorted.dualindex-deduped.sorted.bam")])
            cmds.append([script, os.path.join(barcode_dedupe_output, samp, samp + '.sorted.bam')])
            cmds.append([script, os.path.join(barcode_dedupe_output, samp, samp + '.sorted.dualindex-deduped.sorted.bam')])
        demux.forkpoollang("sh", cmds, threads)  
        os.system("Rscript %s/downsampling/downsample_plot.R %s %s"%(script_dir, downsampling_pre_folder, barcode_dedupe_output))

    if 7 not in avoid: 
        os.system('echo "Starting to generate freq files: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)     
        doubleprint('echo "Starting to generate freq files: "')

        cmds = []
        script = os.path.join(script_dir,'freq/bam-snvindelfreq.py') 
        for samp in samplenames:
            fls = []
            if freq_deduped:
                fls.append(samp + '.sorted.dualindex-deduped.sorted.bam')
            if freq_nondeduped:                
                fls.append(samp + '.sorted.filt.sorted.bam')
            for fl in fls:
                if re.search('sorted\.bam$',fl) and not re.search('\.duplex\.',fl) and (selector_freq_size[selector_freq[getbed(fl)]] < config['bed_limit'] or config['bed_force_freq']):
                    #Break selector if it's loong
                    if selector_freq_size[selector_freq[getbed(fl)]] > config['bed_break'] and not os.path.exists(os.path.join(abs_path,'000000' + os.path.basename(selector_freq[getbed(fl)]))):
                        breakbed(selector_freq[getbed(fl)])
                    #If we are working with a broken selector
                    if selector_freq_size[selector_freq[getbed(fl)]] > config['bed_break']:
                        allbeds = os.listdir(abs_path)
                        allbeds = [x for x in allbeds if re.match('^[0-9]{6}',x) and os.path.basename(selector_freq[getbed(fl)]) in x]
                        print(allbeds)
                        for sel in allbeds:
                            cmds.append([script,
                                freq_params,
                                os.path.join(barcode_dedupe_output, samp, fl),
                                os.path.join(abs_path,sel),
                                selector[getbed(fl)],
                                genome
                                ])
                    #Otherwise
                    else:
                        cmds.append([script,
                            freq_params,
                            os.path.join(barcode_dedupe_output, samp, fl),
                            selector_freq[getbed(fl)],
                            selector[getbed(fl)],
                            genome
                            ])
        #Execute        
        demux.forkpoollang("python",cmds,threads)
        
        # Merge if we split them up
        for samp in samplenames:
            bamdir = os.path.join(barcode_dedupe_output, samp)
            allbeds = os.listdir(bamdir)
            allsamples = list(set([x[6:] for x in allbeds if re.match('[0-9]{6}',x)]))
            print(allsamples)
            for sample in allsamples:
                if re.search('\.indels\.',sample):
                    os.system('echo "CHR\tPOS\tDEPTH\tREF\tINDEL\tPLUS\tMINUS\tPAD\tMOTIF\tRECOVERED\tREMOVED" > ' + os.path.join(bamdir,sample))
                else:
                    os.system('echo "CHR\tPOS\tDEPTH\tREF\tR+\tR-\tA+\tA-\tC+\tC-\tT+\tT-\tG+\tG-\tPAD\tMOTIF" > ' + os.path.join(bamdir,sample))
                os.system('cat ' + os.path.join(bamdir,'0*' + sample) + ' >> ' + os.path.join(bamdir,sample))
                os.system('rm ' + os.path.join(bamdir,'0*' + sample))
           
        # Remove split selectors. 
        splitbeds = os.listdir(abs_path)
        splitbeds = [os.path.join(abs_path,x) for x in splitbeds if re.match('^[0-9]{6}',x)]
        if len(splitbeds) > 0: os.system('rm ' + ' '.join(splitbeds))
        
    if 8 not in avoid:
        os.system('echo "Starting STAR-Fusion: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)                 
        cmds = []
        for samp in samplenames:            
            otpfolder = os.path.join(fastp_good_folder,samp)
            os.system('echo "STAR-Fusion on ' + samp + ": " + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
            cmds.append(["%s/mapping/mapping_STAR_Fusion.sh"%(script_dir), otpfolder, STAR_Fusion_CTAT_folder, STAR_fusion_root_folder, subthreads])            
        demux.forkpoollang("sh", cmds, subjobs)    
    
os.system('echo "Done: ' + time.ctime(time.time()) + ' (epoch: ' + str(time.time()) + ')" | tee -a ' + timestampfile)
