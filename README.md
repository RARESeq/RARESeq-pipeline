## 1. Introduction

This tutorial is aimed at RARE-Seq pipeline users. It is structured in
the following sections:

-   [1. Introduction](#1_introduction)
-   [2. Instructions for Pipeline
    Users](#2_instructions_for_pipeline_users)
    -   [2.1 Overview of the RARE-Seq
        Pipeline](#21_overview_of_the_rare_seq_pipeline)
        -   [Stage 1: Expression analysis](#stage-1-expression-analysis)
        -   [Stage 2: Genotyping Analysis](#stage-2-genotyping-analysis)
    -   [2.2 The Pipeline Input](#22_the_pipeline_input)
    -   [2.3 Running the pipeline](#23_the_pipeline_run)
    -   [2.4 The Pipeline Output](#24_the_pipeline_output)

## 2. Instructions for Pipeline Users

### 2.1 Overview of the RARE-Seq Pipeline

The RARE-Seq pipeline is organized into the steps described below:

1.  **Demultiplexing**: Splits the parent FASTQ files into individual
    FASTQ files for each sample.
2.  **Fastp Processing**: Performs adapter sequence trimming, quality
    control and filtering of the sequencing reads.
3.  **Mapping**: Aligns the reads to the reference genome with STAR,
    using a standard 2-pass mapping strategy, and extracts uniquely
    mapping reads.
4.  **Pre-dedupe QC**: Performs various quality checks, including
    estimating the DNA contamination rate, rRNA percent estimation, and
    other stats.
5.  **Deduping**: Filters out off-target reads and collapses duplicated
    on-target reads into barcode families. A barcode family is defined
    as a set of reads that share the same UMI, start positions for them
    and their respective read mates, as well as the splicing signature
    for them and their respective mates.
6.  **Expression Quantification**: Applies RSEM on the deduped data to
    quantify gene expression.
7.  **Downsampling Analysis**: Runs a downsampling analysis to estimate
    the library complexity. This analysis includes the following steps:
    -   downsampling of the pre-deduping **on- and off-target** reads
        from each sample to predefined depth levels.
    -   deduping of the downsampled reads.
    -   calculating a downsampling curve with the number of input reads
        on the x-axis and the number of on-target barcode families on
        the y-axis.
8.  **Freq file generation**: Outputs freq files that can be used for
    variant calling. This step is disabled by default for whole-exome
    selectors, but can be enabled as described below.
9.  **Fusion Calling**: Runs a separate mapping of the FASTQ files with
    STAR-Fusion to identify fusion candidates.

### 2.2 The Pipeline Input

The example input data from this tutorial can be downloaded from the
[RARE-seq website](https://rareseq.stanford.edu/).

The pipeline needs 3 mandatory and 1 optional input files:

1.  FASTQ R1 file (required): The input FASTQ R1 file, in .fastq.gz
    format. This file can (and most of the times should) be a symlink:

        ls ../HiSeq1574_dowsampled/*1.fq.gz

        ## ../HiSeq1574_dowsampled/Hi1574_CKDL200169161-1a_HFJ7VBBXX_L5_1.fq.gz

2.  FASTQ R2 file (required): The input FASTQ R2 file, in .fastq.gz
    format. This file can (and most of the times should) be a symlink:

        ls ../HiSeq1574_dowsampled/*2.fq.gz

        ## ../HiSeq1574_dowsampled/Hi1574_CKDL200169161-1a_HFJ7VBBXX_L5_2.fq.gz

3.  Sample2barcode file (required): A tab-delimited file with the
    following columns:

    1.  Sample name.
    2.  Sample barcode.
    3.  Sample type. Currently restricted to: cfRNA (or cfrna), Tumor
        (or tumor), Normal (or normal).
    4.  The barcode scheme. Currently restricted to: Flex (or flex).
    5.  The selector bed file.  
        **Please note that:**

    -   **this file should not have any headers.**
    -   **columns should be separated by a <u>single</u> tab character,
        not space, or other separators.**
    -   **the combination of sample names/sample types should be unique
        across rows (i.e., column 1 + column 3 should be unique)**.  
    -   **the name of the sample2barcode file should include the “s2b”
        or “sample2barcode” strings.**  
        Here is an example sample2barcode file:

    <!-- -->

        s2b = read.delim("../HiSeq1574_dowsampled/sample2barcode_hiseq1574.txt", header = F)
        cat(apply(head(s2b), 1, paste, collapse = "\t"), sep = "\n")

        ## CTR73.C1.C5  AGGACATT-AGGACATT   cfrna   Flex    CFRNA_V1_withERCC.bed
        ## LUP341.C1.C2 CGGATTCA-CGGATTCA   cfrna   Flex    CFRNA_V1_withERCC.bed
        ## LUP1003.C1.C2    TATGGCAC-TATGGCAC   cfrna   Flex    CFRNA_V1_withERCC.bed
        ## CTR45.C1.C3  GATGAGGG-GATGAGGG   cfrna   Flex    CFRNA_V1_withERCC.bed
        ## CTR135.C2.C5 CCATAGGA-CCATAGGA   cfrna   Flex    CFRNA_V1_withERCC.bed
        ## LUP782.C1.C2 CTAGTAGC-CTAGTAGC   cfrna   Flex    CFRNA_V1_withERCC.bed

4.  Configuration file (optional). This is a configuration file
    specifying the pipeline parameters. If not provided, the default
    configuration output below is used:

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

    To provide a custom configuration to the pipeline, you need to copy
    the default configuration file from the following path:
    `configs/config_docker.py`into the folder containing the lane’s
    inputs.  
    Here we illustrate several parameters that you can modify in the
    configuration file to adjust the behavior of the pipeline. Only
    modify the rest if you know what you are doing, as it can lead to
    the pipeline crashing:

    -   Pipeline steps to skip:

    <!-- -->

         # Steps to skip. Define as vector integres
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

        'avoid': []

    -   Extension of the selector bed file. The extension of the
        selector file used for filtering for on-target reads. By default
        a 500bp padding is added around the target selector positions.
        If set to ’’, no padding is used. If other padding is set, you
        need to make sure that the corresponding bed file is available
        in the folder containing the selectors.

    <!-- -->

        'bed_ext': 'add500bp',

    -   Parameters controlling the generation of freq files. By default,
        freq files for samples whose selector is larger than
        *bed\_limit* base-pairs are not generated. You can set
        *bed\_force\_freq* to “True” to force the freq file generation
        in such cases. Please not that doing so will significantly
        increase the running time of the lanes that contain samples
        captured with large selectors:

    <!-- -->

        'bed_force_freq': False, 
        'bed_limit': 100000000,

    -   The downsampling levels:

    <!-- -->

        #Downsampling parameters
        'downsampling_breaks' : [x for x in [1e6, 2.5e6, 5e6, 1e7, 2.5e7, 5e7, 1e8, 2.5e8, 5e8]],

All these files should be available in the same folder,
e.g. `/silo5/cfRNA/v3.3/HiSeq1805/`.

### 2.3 Running the Pipeline

Run the pipeline using this command, where each of the 4 inputs are
described in the previous section:

    python pipeline-run.py $1 $2 $3 $4

### 2.4 The Pipeline Output

The output of the pipeline should include the following files inside the
root directory, e.g `../HiSeq1574_dowsampled/`:

-   RSEM outputs containing expression quantification in
    `exp_analysis/expression/barcode-deduped/*.genes.results`, e.g:

        cat ../HiSeq1574_dowsampled/exp_analysis/expression/barcode-deduped/Sample_LUP1004.C1.C2_cfrna.genes.results | head

        ## gene_id  transcript_id(s)    length  effective_length    expected_count  TPM FPKM
        ## C-ERCC-1 C-ERCC-1    93.00   0.00    0.00    0.00    0.00
        ## C-ERCC-10    C-ERCC-10   144.00  7.27    52.00   140325.59   286885.47
        ## C-ERCC-11    C-ERCC-11   145.00  7.50    0.00    0.00    0.00
        ## C-ERCC-12    C-ERCC-12   151.00  9.01    1.00    2177.62 4451.98
        ## C-ERCC-13    C-ERCC-13   159.00  11.25   1.00    1743.47 3564.39
        ## C-ERCC-14    C-ERCC-14   126.00  3.73    4.00    21026.06    42986.26
        ## C-ERCC-16    C-ERCC-16   202.00  27.53   2.00    1424.88 2913.06
        ## C-ERCC-17    C-ERCC-17   211.00  31.80   4.00    2466.96 5043.51
        ## C-ERCC-18    C-ERCC-18   202.00  27.53   87.00   61982.23    126718.16

-   The deduped bam files in
    `exp_analysis/deduped/output/barcode-deduped/Sample*/Sample*.sorted.dualindex-deduped.sorted.bam`,
    e.g.:

        samtools view ../HiSeq1574_dowsampled/exp_analysis/deduped/output/barcode-deduped/Sample_LUP1004.C1.C2_cfrna/Sample_LUP1004.C1.C2_cfrna.sorted.dualindex-deduped.sorted.bam | head

        ## K00124:1:1:637:HFJ7VBBXX:5:1101:8197:14635:I14+0+0:J14+0+0:U12+0+0+0:V22+0+0+0   163 chr1    879752  255 24M294N96M  =   880930  1813    GAGGTGGTGGAAGGGGCCAGGGGCTGCCTCAGTCGTCCTCTGAGAGCTGCAGATCCTCCAGCTCGTCCTCCGGCCCCTGGGCCAGCTGCTGCAGCTCCCCAGGGGCCAGCCCCGCCTCTG    JJJJJJJJJJJJJJJJJFJJF<JJFJJJJJJJJJJFFJJJJJJJJJJJ<JJJJJJJJJJJJFJFJJJJJJJJJJJJAJJJJJJJJJJJJJFJJJJJJJJJ7JJJJJJFJJAJFFJFF<FJ    NH:i:1  HI:i:1  AS:i:231    nM:i:0  MC:Z:104M519N14M    NM:i:0  MD:Z:120
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:17279:2276:I14+0+0:J14+0+0:U22+0+0+0:V30+0+0+0   163 chr1    879928  255 120M    =   881631  1938    ATGCTGAAAAATAAATAATAAAGCCTGTCCCGTGTCTACTGCCTCCCCCAACTGCACAGACGCCAGCCTCTAGGCCTGACTGCCAGGGAGGTGGAAACACTGGCCACCAGCCCGGCAGCC    JFJJJJJJJJJJJJJJFJJFJJFFAA<JFJFAAFJJJJFJJJJJJJFJFJJJJJJJJFJJJ7JJJJJJJFJJJJJJJJAFJJAFF-AA-AF<F<AFJAFJJJAFAFJJJ)77<A<FJJJ-    NH:i:1  HI:i:1  AS:i:239    nM:i:0  MC:Z:36M115N84M NM:i:0  MD:Z:120
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:13321:7655:I14+0+0:J14+0+0:U26+0+0+0:V06+0+0+0   163 chr1    879948  255 120M    =   880131  544 AAGCCTGTCCCGTGTCTACTGCCTCCCCCAACTGCACAGACGCCAGCCTCTAGGCCTGACTGCCAGGGAGGTGGAAACACTGGCCACCAGCCCGGCAGCCCCTACAGGCCCCCCAGATGG    JJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFJJJJJJJJJ<JJJ    NH:i:1  HI:i:1  AS:i:240    nM:i:0  MC:Z:50M241N70M NM:i:0  MD:Z:120
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:5690:3401:I14+0+0:J14+0+0:U11+0+0+0:V11+0+0+0    99  chr1    880093  255 88M256N31M  =   880110  392 TGCAGATCCTCCAGCTCGTCCTCCGGCCCCTGGGCCAGCTGCTGCAGCTCCCCAGGGGCCAGCCCCGCCTCTGCGTCTGGGTCTCCATCCTCCGAGTTGCTGCTGTCCTCCTCGCCCTC JAJJJJJJJJJJJJJJJJJJJAJJJFJJFJJJJJJJAJJJJJFJJJJF7A-FAJJJ<JJJJJJJJJJJJJJJJJJJJ<JJJJJFJJJFJJJ<JJJJJJFJAJJJJJJJFFAFF<JJJFF NH:i:1  HI:i:1  AS:i:240    nM:i:0  MC:Z:71M256N48M NM:i:0  MD:Z:119
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:5690:3401:I14+0+0:J14+0+0:U11+0+0+0:V11+0+0+0    147 chr1    880110  255 71M256N48M  =   880093  -392    GTCCTCCGGCCCCTGGGCCAGCTGCTGCAGCTCCCCAGGGGCCAGCCCCGCCTCTGCGTCTGGGTCTCCATCCTCCGAGTTGCTGCTGTCCTCCTCGCCCTCCTCCTCGTCCTCTTCAT <JJFFJJJJFJJJJJFFJJJJJJJJJJJJJJFJJJJJJJAJJJJJJJJJJJJJJJJJJJJJJFJJJJJJFAJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJ NH:i:1  HI:i:1  AS:i:240    nM:i:0  MC:Z:88M256N31M NM:i:0  MD:Z:119
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:13321:7655:I14+0+0:J14+0+0:U26+0+0+0:V06+0+0+0   83  chr1    880131  255 50M241N70M  =   879948  -544    CTGCTGCAGCTCCCCAGGGGCCAGCCCCGCCTCTGCGTCTGGGTCTCCATCCCAAGACCATTCACCCTCCGAGTTGCTGCTGTCCTCCTCGCCCTCCTCCTCGTCCTCTTCATCGTCTTC    JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJFJJJJJJJJJJJJJ    NH:i:1  HI:i:1  AS:i:240    nM:i:0  MC:Z:120M   NM:i:0  MD:Z:120
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:2950:4491:I14+0+0:J14+1+0:U19+0+0+0:V18+0+0+0    163 chr1    1956958 255 120M    =   1956984 145 AGTACACCATGACGGTGTTCCTGCACCAGAGCTGGCGGGACAGCAGGCTCTCCTACAACCACACCAACGAGACCCTGGGTCTGGACAGCCGCTTCGTGGACAAGCTGTGGCTGCCCGACA    F<AJF<FJJ-AJJJJAJJJJJFJJJJJJJJFJ-AFF<JJJJJJJJJJJ-FJJJFFFJFJJFJJJJJJJFJ<FJJJAJFJA7-A-AFA<AJFJJAF-7-7FFA<-A<<J7FF)7AJJF7J<    NH:i:1  HI:i:1  AS:i:237    nM:i:0  MC:Z:119M   NM:i:0  MD:Z:120
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:2950:4491:I14+0+0:J14+1+0:U19+0+0+0:V18+0+0+0    83  chr1    1956984 255 119M    =   1956958 -145    CAGAGCTGGCGGGACAGCAGGCTCTCCTACAACCACACCAACGAGACCCTGGGTCTGGACAGCCGCTTCGTGGACAAGCTGTGGCTGCCCGACACCTTCATCGTGAACGCCAAGTCGGC FAFAA<FJJFFJFJFJJFJFA-F7JJFJJJFFAA-JJJJFFFF7FFJJJJJJF<7JJAJF-JFJFAJJFFJJJJJJJJ7FJJJJAFJJJJJFFJFFFJFAJJJFJ<FJJJFJJFFFAA7 NH:i:1  HI:i:1  AS:i:237    nM:i:0  MC:Z:120M   NM:i:0  MD:Z:119
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:8724:2258:I14+0+0:J14+0+0:U11+0+0+0:V11+0+0+0    83  chr1    1956990 255 103M    =   1957000 -93 TTGCGGGACAGCAGGCTCTCCTACAACCACACCAACGAGACCCTGGGTCTGGACAGCCGCTTCGTGGACAAGCTGTGGCTGCCCGACACCTTCATCGTGAACG JJJJJJJJJJJJJJJJ7AJJJJFJJJJJJJJFJFJJJJJJJFJJJFJFJFFFFJFJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJFF7FJJJJJJJJJF7-J NH:i:1  HI:i:1  AS:i:202    nM:i:0  MC:Z:103M   NM:i:1  MD:Z:1G101
        ## K00124:1:1:637:HFJ7VBBXX:5:1101:8724:2258:I14+0+0:J14+0+0:U11+0+0+0:V11+0+0+0    163 chr1    1957000 255 103M    =   1956990 93  GCAGGCTCTCCTACAACCACACCAACGAGACCCTGGGTCTGGACAGCCGCTTCGTGGACAAGCTGTGGCTGCCCGACACCTTCATCGTGAACGCCAAGTCGGC JJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJ<JJ<F<FFJJJJJ<FJF<JJJJJJJ<FFJJJJ<JFAJJJJJFJJFFJAJFFJJFJJFJAFF-<F7F NH:i:1  HI:i:1  AS:i:202    nM:i:0  MC:Z:103M   NM:i:0  MD:Z:103

-   The freq and indel files in
    `exp_analysis/deduped/output/barcode-deduped/Sample*/Sample*.sorted.dualindex-deduped.sorted.freq.paired.Q30.txt`
    and
    `exp_analysis/deduped/output/barcode-deduped/Sample*/Sample*.sorted.dualindex-deduped.sorted.indels.paired.Q30.txt`,
    respectively, e.g.:

        cat ../HiSeq1574_dowsampled/exp_analysis/deduped/output/barcode-deduped/Sample_LUP1004.C1.C2_cfrna/Sample_LUP1004.C1.C2_cfrna.sorted.dualindex-deduped.sorted.freq.paired.Q30.txt | head
        cat ../HiSeq1574_dowsampled/exp_analysis/deduped/output/barcode-deduped/Sample_LUP1004.C1.C2_cfrna/Sample_LUP1004.C1.C2_cfrna.sorted.dualindex-deduped.sorted.indels.paired.Q30.txt | head

        ## CHR  POS DEPTH   REF R+  R-  A+  A-  C+  C-  T+  T-  G+  G-  PAD MOTIF
        ## chr1 879752  1   G   1   0   0   0   0   0   0   0   0   0   219 AGAGG
        ## chr1 879753  1   A   1   0   0   0   0   0   0   0   0   0   220 GAGGT
        ## chr1 879754  1   G   1   0   0   0   0   0   0   0   0   0   221 AGGTG
        ## chr1 879755  1   G   1   0   0   0   0   0   0   0   0   0   222 GGTGG
        ## chr1 879756  1   T   1   0   0   0   0   0   0   0   0   0   223 GTGGT
        ## chr1 879757  1   G   1   0   0   0   0   0   0   0   0   0   224 TGGTG
        ## chr1 879758  1   G   1   0   0   0   0   0   0   0   0   0   225 GGTGG
        ## chr1 879759  1   T   1   0   0   0   0   0   0   0   0   0   226 GTGGA
        ## chr1 879760  1   G   1   0   0   0   0   0   0   0   0   0   227 TGGAA
        ## CHR  POS DEPTH   REF INDEL   PLUS    MINUS   PAD MOTIF   RECOVERED   REMOVED
        ## chr1 20915589    10  T   -C  4   1   33  TTCCC   0   0
        ## chr1 24841007    30  C   -T  0   1   0   ACTTT   0   0
        ## chr1 114429186   9   G   -T  1   0   0   CGTTT   0   0
        ## chr1 151819884   8   A   +A  0   1   0   GAAAA   0   0
        ## chr1 153346251   161 G   -AAA    1   0   51  AGAAA   0   0
        ## chr1 156715104   29  A   +A  0   1   0   GAAAA   0   0
        ## chr10    17271931    175 C   -G  0   1   0   GCGTC   0   0
        ## chr10    17275871    198 C   -A  1   0   0   GCAAT   0   0
        ## chr10    17277260    188 T   -GA 1   0   0   ATGAG   0   0
