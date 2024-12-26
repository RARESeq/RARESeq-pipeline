config = {
    # The demultiplexing step splits the parent FASTQ into muliple files in order to 
    # multithread and accelerate the process. The parent fastq is split into files with
    # this number of lines per file
    'demux_split_lines': 10000000,

    # Steps to skip in Stage 1. Define as vector integres
    # e.g. to skip steps 0 and 1 define 'avoid: [0,1]'
    # (-1) Skip all
    # () Skip none (default)
    # (0) Demultiplex
    # (1) fastp
    # (2) mapping
    # (3) pre-dedupe QC (DNA contamination, molecule length, read stats)
    # (4) deduping
    # (5) expression quantification
    # (6) downsampling analysis
    # (7) freq file generation 
    # (8) fusion calling   
    
    'avoid': [], 

    # Placeholder for stage 2. Not used yet.
    'avoid2': [],

    # fastp parameters. Please refer to https://github.com/OpenGene/fastp for
    'afterp_overlap_len_require': '30',
    'afterp_overlap_diff_limit': '5',
    'afterp_length_required': '35',
    'afterp_qualified_quality_phred': '30',
    'afterp_unqualified_percent_limit': '20',
    'afterp_trim_r1_front' : '10',
    'afterp_trim_r1_tail' : '10',
    'afterp_trim_r2_front' : '10',
    'afterp_trim_r2_tail' : '10',

    #rnaseqc params
    'rnaseqc_cmd' : 'rnaseqc.v2.3.5.linux',
    'rnaseqc_gtf' : 'gencode.v24lift37.annotation.plus.ERCC.forRNASeQC.gtf',
    'rnaseqc_subsample_limit' : '1000000',
    'DNA_contamination_exon_positions' : '/indexes/exonic_positions/exons_75bp.bed',
    
    # Path to reference genome used in mapping process and STAR index
    'fa': 'hg19_ERCC.fa',
    'gtf': 'gencode.v24lift37.annotation.plus.ERCC.gtf',
    'protein_coding_genes': 'gencode.v24lift37.annotation.plus.ERCC_protein_coding.txt',
    #if the STAR index doesn't exist, it will be created by the pipeline (which might take up to 1h)
    'star_index': 'star_index',
    'STAR_Fusion_CTAT_folder': '/indexes/CTAT_folder/GRCh37_gencode_v19_CTAT_lib_Mar272019.source/ctat_genome_lib_build_dir/',

    # Path to tools in specified vesions. Be aware that modifying these may lead to
    # errors as flags and other tings may change from version to version. 
    'samtools': 'samtools',
    'bedtools': 'bedtools',
    'STAR': "STAR",
    'RSEM_prepare_reference' : "rsem-prepare-reference",
    'sambaba_cmd' : "sambamba-0.7.1-linux-static",
    
    'bed_ext': 'add500bp',    
    'bed_force_freq': False,
    'bed_limit': 100000000,
    'bed_break': 5000000,

    #Selector parameters
    'selectors_dir' :  "/selectors/",

    # Additional command line parameters passed to the deduping script
    'dedupe_params' : "",

    # Additional command line parameters passed to the script generating freq files
    'freq_params': "-i '-f3 -F2304'",
    'freq_nondeduped': False,
    'freq_deduped': True,

    # Some steps in the pipeline have internal parallelization (e.g. fastqc). These parameters
    # split these processes into a given number of jobs with a given number of threads. For example, subjobs=2 
    # and subthreads=6 will run two processes at a time with 6 threads each. These parameters would be most 
    # beneficial for a lane with two samples. If either parameter is set to 'auto' these will be calculated
    # based on the number of samples in the lane.
    'subjobs': 'auto',
    'subthreads': 'auto',

    # Maximum number of threads used. Define as string.
    'threads': 'auto',

    #Downsampling parameters
    'downsampling_breaks' : [x for x in [1e6, 2.5e6, 5e6, 1e7, 2.5e7, 5e7, 1e8, 2.5e8, 5e8]],
    'downsampling_exon_positions' : '/indexes/exonic_positions/exons_merged.bed'
}
