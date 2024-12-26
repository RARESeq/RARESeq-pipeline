#!/usr/bin/python
"""
Demultiplex CAPP-Seq paired-end FastQ files using error correcting codes.

Input:
  - R1.fastq[.gz]
  - R2.fastq[.gz]
  - I1.fastq[.gz]
  - I2.fastq[.gz]
  - sample2barcode.txt (tab delimited file with sample to adapter associations (one per line))
Output:
  1. demultiplexed fastq files (R1/R2) for each barcode given in the sample-barcode file
  2. demultiplexed fastq files (R1/R2) with all unidentified reads
  3. text file with statistics about error correction
Notes:
- If input file names end in .gz, assume gzipped input
- If index read files are not present, the script will attempt to read the barcodes from the read headers
- Reads are demultiplexed only if sample could be identified and UMIs are valid (after error correction)
- For demultiplexed reads, the stubs with the UMIs are trimmed to the length specified in the umis.txt file
- For the unidentified reads, the maximum stub length is trimmed off and the trimmed sequences are written to the header
Dependencies:
- ParseFastQ.py - library for reading FastQ files
- umis.txt - text file with list of insert UMIs (columns: number, umi, stag, punct, full_insert, length)
- codetable_umis.txt - text file with umi codewords (columns: message, codeword_number, hamming_distance, number_of_Ns, rev_comp_msg)
- codetable_samples.txt - text file with index/sample codewords (column: message, codeword_number, hamming_distance, number_of_Ns, rev_comp_msg)
See also:
- minimal_decode_demultiplex.py - minimal version with same input/output but without error statistics (faster)
- parallel_decode_demultiplex.py - multi-threaded implementation based on minimal version (i.e. no error stats)
- generate_code_table.py - script for generating code table files

Version history:
2017-07-09 v1.0 <stehr@stanford.edu>
2017-08-17 v1.1 added support for MiSeq/NextSeq data where I2 is not reverse complemented

TODO:
- error correct phix signatures
- add parameter datafiles_dir
- optimize runtime using profiling
"""

import sys, os, re, gzip, datetime
import ParseFastQ as fq
from operator import itemgetter

# ----------- Parameters -----------

# settings
i2_reverse_complement_auto_detect = True # if true, auto-detect from FastQ files, otherwise use setting below
i2_reverse_complement = True # if not auto-detect: whether I2 is reverse complement (e.g. HiSeq4000) or not (e.g. MiSeq/NextSeq)
i2_auto_detect_reads = 100 # if auto-detect: number of reads to check in R2 or I2

demultiplex = True   # if true, write demultiplexed R1,R2 fastq files for each sample, otherwise count only
write_unidentified = False # if true, and demultiplex=True, write unidentified reads to separate files (see below)
force_uncompressed = True # if true, write uncompressed FastQs (much faster), otherwise write gzipped output if input file has .gz extension
write_stats_file = True   # if true, write a file with demultiplexing statistics (see stats_filename below), otherwise write to stdout
verbosity = 0 # set level of verbosity in console output (0=silent, 1=minimal, 2=verbose, 3=debug)

# error stats
SHOW_MAX_UMIS = 40 # number of umi1/umi2 to show in stats table
SHOW_MAX_I1I2 = 40 # number of i1/i2 to show in stats table

# output files
stats_filename = "demux_stats.txt"                  # statistics log file
demux_file_template = "{0}_{1}_S{2:02}_{3}_{4}.fastq"   # demultiplexed fastq files
# {0}=sample_name                                   {1}=sample_number
# {2}=R1/R2/I1/I2                                   {3}=barcode
unid_file_template = "Unidentified_{0}.fastq"       # catch-all file for non-demultiplexed reads
# {0}=R1/R2/I1/I2
fastq_header_demux = "I{0:02}+{1}+{2}:J{3:02}+{4}+{5}:U{6:02}+{7}+{8}+{12}:V{9:02}+{10}+{11}+{13}" # this will be appended to the R1/R2 headers of demultiplexed reads
# {0}=I1 codeword number                            {1}=I1 number of corrected errors
# {2}=I1 number of corrected Ns                     {3}=I2 codeword number
# {4}=I2 number of corrected errors                 {5}=I2 number of corrected Ns
# {6}=UMI1 codeword number                          {7}=UMI1 number of corrected errors
# {8}=UMI1 number of corrected Ns                   {9}=UMI2 codeword number
# {10}=UMI2 number of corrected errors              {11}=UMI2 number of corrected Ns
# {12}=I1 errors in punctuation mark                {13}=I2 errors in punctuation mark
fastq_header_unid = "{0}+{1}:{2}:{3}" # this will be appended to the R1/R2 headers of unidentified reads
# {0}=I1 index sequence                             {1}=I1 index sequence
# {2}=R1 trimmed stub sequence                      {3}=R2 trimmed stub sequence

# required input files
#adapter_file = "umis.txt"
#codetable_samples_file = "codetable_samples.txt"
#codetable_umis_file = "codetable_umis.txt"

envar_name = 'pipeline_dir' # environment variable pointing to parent directory of datafiles
adapter_file = os.path.join(os.environ[envar_name],"configs", "umis.txt")
codetable_samples_file = os.path.join(os.environ[envar_name],"configs", "codetable_samples.txt")
codetable_umis_file = os.path.join(os.environ[envar_name],"configs", "codetable_umis.txt")

# PhiX index sequences
PHIX_I1 = "TCGAATGA" # specific to our adapter design
PHIX_I2 = "CGCATGTC" # off sequencer sequence

# ----------- Constants -----------

VERSION = "1.1" # version of demultiplexing script
CODETABLE_HEADER = "#CAPP-Seq code table v1.0" # expected header in codetable files
UMI_HEADER = "#CAPP-Seq insert UMIs v1.0" # expected header in umi file
ADAPTER_TYPE = "Flex" # expected adapter type in column 4 of sample2barcode file
INDEX_BARCODE_LENGTH = 8 # length of single index barcode, for index extraction from header
MAX_STUB_LENGTH = 12 # length of longest stub sequence (for trimming unidentified reads)
DUAL_INDEX_DELIMITER = '-' # delimiter in sample2barcode file, e.g. 'ACGTGATC-ACGTGATC'
CODEWORD_NOT_FOUND = None # value for decoded codeword number if message could not be decoded

# ----------- Functions -----------

"""Prints error message to stderr and exits"""
def error_and_exit(msg, exit_code = 1):
  sys.stderr.write("ERROR: {0}\n".format(msg))
  sys.exit(exit_code)

"""Prints a non-critical warning to stderr"""
def warning(msg):
  sys.stderr.write("WARNING: {0}\n".format(msg))

"""Returns a percentage as a formated string"""
def pct(enu,den,digits=1):
    if den == 0:
      return "NA"
    else:
      return ("{0:5."+str(digits)+"f}%").format(100.0 * enu / den)

"""Returns the complementary DNA sequence"""
def complement(seq):
    """Returns the complement DNA sequence (allowing for IUPAC ambiguity codes)"""
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))

"""Returns the reverse-complement DNA sequence"""
def reverse_complement(seq):
    """Returns the reverse complement DNA sequence (allowing for IUPAC ambiguity codes)"""
    return complement(seq)[::-1]

"""Returns the items in a dictionary as a list of tuples sorted by descending value"""
def dict_items_by_val(d):
  return sorted(d.items(), key = lambda item: item[1], reverse=True)

"""Prints error and exits if file is not readable"""
def check_infile(f):
  if not (os.path.isfile(f) and os.access(f, os.R_OK)):
    sys.stderr.write("ERROR: Can not read from file %s\n" %f)
    sys.exit(2)

"""Prints error and exits if directory is not writable"""
def check_outdir(d):
  if not (os.path.isdir(d) and os.access(d, os.W_OK)):
    sys.stderr.write("ERROR: Can not write to directory %s\n" %d)
    sys.exit(3)

# Output file names

"""Returns the name of the stats output file"""
def get_stats_filename(out_dir, file_name):
  return os.path.join(out_dir, file_name[-14:-6] + '_' + stats_filename)

"""Returns the name for a demultiplexed sample fastq output file.
   Can be customized by the constant demux_file_template.
   read_id = {R1,R2,I1,I2}
"""
def get_demux_filename(out_dir, file_name, sample_name, sample_number, read_id, barcode):
  return os.path.join(out_dir, demux_file_template.format(file_name.replace('.fastq',''), sample_name, sample_number, read_id, barcode))

"""Returns the name of the unidentitied reads fastq output file.
   read_id = {R1,R2,I1,I2}
"""
def get_unidentified_filename(out_dir, read_id):
  return os.path.join(out_dir, unid_file_template.format(read_id))

# ------------ Classes ------------

"""Error-correcting insert-UMI stubs"""
class Stub:
  def __init__(self, num, umi, stag, punct, full_insert, length):
    self.num = int(num) # umi number
    self.name = "stub{0:02}".format(num)
    self.umi = umi # 6bp UMI sequence
    self.stag = stag # staggered bases (0-3)
    self.punct = punct # punctuation mark bases (GT or CT)
    self.full_insert = full_insert # full insert sequence
    self.length = int(length) # length of full insert sequence
    assert self.full_insert[1:] == self.umi + self.stag + self.punct, "Error in UMI definitions: full_insert != umi+stag+punct"
    assert len(full_insert) == self.length, "Error in UMI definitions: lenght does not match full_insert sequence"

"""List of stubs with insert UMIs"""
class Stubs:

  """reads the stub definitions from the UMI file"""
  def __init__(self, adapter_file):
    self.stubs = []
    self.umi2stub = {}
    i = 1
    with open(adapter_file) as f:
      header = f.readline().rstrip()
      if not header == UMI_HEADER:
        sys.stderr.write("WARNING: Header in UMI file does not match. Expected '{0}' but found '{1}'\n".format(UMI_HEADER, header))
      cols = f.readline().rstrip().split()
      for line in f:
        fields = line.rstrip().split('\t')
        num = int(fields[0])
        assert num == i, "Warning: Unexpected stub number"
        umi = fields[1]
        stag = fields[2]
        punct = fields[3]
        full_insert = fields[4]
        length = int(fields[5])
        s = Stub(num, umi, stag, punct, full_insert, length)
        self.stubs.append(s)
        self.umi2stub[umi] = s
        i += 1

  """Returns a stub object by UMI sequence or None if sequence is invalid"""
  def getByUmi(self, seq):
    return self.umi2stub.get(seq, None)

  """Returns a list of all adapter objects"""
  def getAll(self):
    return self.stubs

  """Returns a list of all valid UMI sequences"""
  def getValidUmis(self):
    return self.umi2stub.keys()

  """Returns the length of the longest stub sequence (for trimming)"""
  def getMaxLength(self):
    return max([s.length for s in self.stubs])


class CodeTableEntry:
  """An entry in the CodeTable"""
  def __init__(self, msg, codeword, dist, ns, revcomp):
    self.msg = msg
    self.codeword_num = int(codeword)
    self.dist = int(dist)
    self.ns = int(ns)
    self.revcomp = revcomp

  def __str__(self):
    return "{0:02}+{1}".format(self.codeword_num, self.dist)

"""A code table mapping from messages to decoded code words"""
class CodeTable:

  """Loads the codetable from the given file (generated from the adapter set)"""
  def __init__(self, codetable_file, i2_is_reverse_complement = False):

    self.codewords = {} # dictionary of valid codewords
    self.codetable = {} # mapping from message to codeword

    # Index read 2 can be reverse complemented in the FastQ data (e.g. HiSeq4000) or not (e.g. MiSeq/NextSeq)
    # Depending on this setting, the codetable_i2 will contain the reverse complement or raw message strings, respectively.
    self.i2_is_reverse_complement = i2_is_reverse_complement
    self.codetable_i2 = {} # mapping from message to codeword for index read 2

    # read codetable
    with open(codetable_file) as f:
      header = f.readline().rstrip()
      if not header == CODETABLE_HEADER:
        sys.stderr.write("WARNING: Header in code table file does not match. Expected '{0}' but found '{1}'\n".format(CODETABLE_HEADER, header))
      cols = f.readline().rstrip().split()
      i = 1
      for line in f:
        fields = line.rstrip().split()
        msg = fields[0]
        codeword_num = int(fields[1])
        dist = int(fields[2])
        ns = int(fields[3])
        revcomp = fields[4]
        entry = CodeTableEntry(msg, codeword_num, dist, ns, revcomp)
        # message -> codeword index
        if msg in self.codetable:
          sys.stderr.write("Warning: Message {0} already in codetable\n".format(msg))
        self.codetable[msg] = entry

        if self.i2_is_reverse_complement:
          # revcomp(message) -> codeword index
          if revcomp in self.codetable_i2:
            sys.stderr.write("Warning: Message {0} already in reverse complement codetable\n".format(revcomp))
          self.codetable_i2[revcomp] = entry
        else:
          # ignore reverse complements and store raw message in i2 table
          self.codetable_i2[msg] = entry

        if dist == 0:
          # codeword index -> codeword
          if msg in self.codewords:
            sys.stderr.write("Warning: Codeword {0} already in codeword list\n".format(msg))
          self.codewords[codeword_num] = msg

  """Returns the decoded CodeTableEntry if message could be decoded or CODEWORD_NOT_FOUND"""
  def decode(self, message):
      return self.codetable.get(message,CODEWORD_NOT_FOUND)

  """Decode index read 2 which can be reverse complemented or not.

     Behaviour depends on self.i2_is_reverse_complement. Otherwise, it is the same as self.decode(), i.e.
     Returns the decoded CodeTableEntry if message could be decoded or CODEWORD_NOT_FOUND.
  """
  def decode_i2(self, message):
      return self.codetable_i2.get(message,CODEWORD_NOT_FOUND)

  """Returns the number of entries in the code table"""
  def getTableSize():
      return len(self.codetable)

  """Returns a codeword by number"""
  def getCodewordByNum(self, num):
     return self.codewords[num]

"""The sample to barcode associations and corresponding output files"""
class Samples2Barcodes():

  # -- constructor --
  """initializes the sample-barcode association table"""
  def __init__(self, s2b_filename, out_dir, use_gzip, use_index_reads, write_unidentified):

    # -- member variables --
    self.out_dir = out_dir
    self.index2barcode_sample = [] # list of barcode-sample associations
    self.barcode2index = {} # lookup table from barcode to index in list above (and in file handle lists)
    self.demux_files_r1 = [] # file handles for per-barcode output fastq R1 files
    self.demux_files_r2 = [] # file handles for per-barcode output fastq R2 files
    self.demux_files_i1 = [] # file handles for per-barcode output fastq I1 files
    self.demux_files_i2 = [] # file handles for per-barcode output fastq I2 files
    self.use_gzip = use_gzip # if True, compress output files with gzip

    self.unidentified_index = 0 # initial index in file handle list (ignored if write_unidentified=False)

    # initializes self.index2barcode_sample, self.barcode2index, self.unidentified_index
    self._read_sample2barcode(s2b_filename)

    # initializes self.demux_files_r{1,2}
    self._init_demux_files(out_dir, use_gzip, use_index_reads, write_unidentified)

  # -- public methods --

  """Returns the number of sample barcodes read from the sample2barcode file"""
  def get_number_of_barcodes(self):
    return len(self.index2barcode_sample)

  """Returns the barcodes from the sample2barcode file as a list"""
  def get_barcodes(self):
    return [ bs[0] for bs in self.index2barcode_sample ]

  """Returns the sample name associated with the given barcode"""
  def barcode2name(self, barcode):
    return self.index2barcode_sample[self.barcode2index[barcode]][1]

  """Returns true iff the parameter is one of the expected sample barcodes read from the sample2barcode file"""
  def is_sample_barcode(self, barcode):
    return barcode in self.barcode2index

  """Writes a read to the demultiplexed fastq file corresponding to the barcode.
    read_num: 1=R1, 2=R2, 3=I1, 4=I2
    This will fail if class was not instantiated with write_unidentified=True
  """
  def write_unidentified_read(self, read_data, read_num):
    if read_num == 1:
      out = self.demux_files_r1[self.unidentified_index]
      out.write(read_data.toString() + "\n")
    elif read_num == 2:
      out = self.demux_files_r2[self.unidentified_index]
      out.write(read_data.toString() + "\n")
    elif read_num == 3:
      out = self.demux_files_i1[self.unidentified_index]
      out.write(read_data.toString() + "\n")
    elif read_num == 4:
      out = self.demux_files_i2[self.unidentified_index]
      out.write(read_data.toString() + "\n")
    else:
      raise Exception("Unknown read number %s. Expected one of 1,2,3,4." % read_num)

  """Writes a read to the demultiplexed fastq file corresponding to the barcode.
    read_num: 1=R1, 2=R2, 3=I1, 4=I2
  """
  def write_demuxed_read(self, read_data, barcode, read_num):
    if read_num == 1:
      out = self.demux_files_r1[self.barcode2index[barcode]]
      out.write(read_data.toString() + "\n")
    elif read_num == 2:
      out = self.demux_files_r2[self.barcode2index[barcode]]
      out.write(read_data.toString() + "\n")
    elif read_num == 3:
      out = self.demux_files_i1[self.barcode2index[barcode]]
      out.write(read_data.toString() + "\n")
    elif read_num == 4:
      out = self.demux_files_i2[self.barcode2index[barcode]]
      out.write(read_data.toString() + "\n")
    else:
      raise Exception("Unknown read number %s. Expected one of 1,2,3,4." % read_num)

  """Closes the open file handles"""
  def close(self):
    for f in self.demux_files_r1:
      f.close()
    for f in self.demux_files_r2:
      f.close()
    for f in self.demux_files_i1:
      f.close()
    for f in self.demux_files_i2:
      f.close()

  # -- private methods --

  """Reads the sample2barcode file and returns a list of (barcode,sample) pairs and a dictionary of barcode to list index"""
  def _read_sample2barcode(self, s2b_filename):
    self.barcode2index = {}
    self.index2barcode_sample = []
    s2b = open(s2b_filename)
    #n = 0
    for n, line in enumerate(s2b):
      fields = line[:-1].split('\t')
      if len(fields) > 1:
        sample = fields[0]
        barcode = fields[1]
        if len(fields) < 4:
          sys.stderr.write("Warning: Expected at least 4 columns but found {0} in line {1} of '{2}'\n".format(len(fields), n+1, s2b_filename))
          sample_type = None
          adapter_type = None
        else:
          sample_type = fields[2]
          adapter_type = fields[3]
          # if adapter_type != ADAPTER_TYPE:
          #   sys.stderr.write("Warning: Expected adapter type '{0}' but found '{1}' in line {2}, column 4 of '{3}'\n".format(ADAPTER_TYPE, adapter_type, n+1, s2b_filename))
        m = re.match(r'([ACGT]{8})' + DUAL_INDEX_DELIMITER + r'([ACGT]{8})',barcode)
        if m:
          i1 = m.group(1)
          i2 = m.group(2)
        else:
          m = re.match(r'([ACGT]{8})', barcode)
          if m:
            # fallback case: expect same barcode in i1 and i2
            i1 = m.group(1)
            i2 = i1
          else:
            sys.stderr.write("Warning: Unexpected barcode format in column 2 of {0}. Ignoring line.\n".format(s2b_filename))
            continue
        key = i1 + DUAL_INDEX_DELIMITER + i2
        self.barcode2index[key] = n
        self.index2barcode_sample.append((barcode,sample))
        #n = n + 1
    self.unidentified_index = n+1 # set index in file handle list to num_samples + 1

  """Opens output files for the demultiplexed reads and stores file handles in self.demux_files_r{1,2}"""
  def _init_demux_files(self, out_dir, use_gzip, use_index_files, write_unidentified):
    # file handles for multiplexed samples
    for i in range(0, len(self.index2barcode_sample)):
      (barcode, sample) = self.index2barcode_sample[i]
      r1_filename = get_demux_filename(out_dir, r1f, sample, i+1, "R1", barcode)
      r2_filename = get_demux_filename(out_dir, r2f, sample, i+1, "R2", barcode)
      if use_index_files:
        i1_filename = get_demux_filename(out_dir, i1f,  sample, i+1, "I1", barcode)
        i2_filename = get_demux_filename(out_dir, i2f, sample, i+1, "I2", barcode)
      if use_gzip:
        self.demux_files_r1.append(gzip.open(r1_filename,'w'))
        self.demux_files_r2.append(gzip.open(r2_filename,'w'))
        if use_index_files:
          self.demux_files_i1.append(gzip.open(i1_filename,'w'))
          self.demux_files_i2.append(gzip.open(i2_filename,'w'))
      else:
        self.demux_files_r1.append(open(r1_filename,'w'))
        self.demux_files_r2.append(open(r2_filename,'w'))
        if use_index_files:
          self.demux_files_i1.append(open(i1_filename,'w'))
          self.demux_files_i2.append(open(i2_filename,'w'))
    # file handles for unidentified reads (append after samples)
    if write_unidentified:
      r1_filename = get_unidentified_filename(out_dir, "R1")
      r2_filename = get_unidentified_filename(out_dir, "R2")
      if use_index_files:
        i1_filename = get_unidentified_filename(out_dir, "I1")
        i2_filename = get_unidentified_filename(out_dir, "I2")
      if use_gzip:
        self.demux_files_r1.append(gzip.open(r1_filename,'w'))
        self.demux_files_r2.append(gzip.open(r2_filename,'w'))
        if use_index_files:
          self.demux_files_i1.append(gzip.open(i1_filename,'w'))
          self.demux_files_i2.append(gzip.open(i2_filename,'w'))
      else:
        self.demux_files_r1.append(open(r1_filename,'w'))
        self.demux_files_r2.append(open(r2_filename,'w'))
        if use_index_files:
          self.demux_files_i1.append(open(i1_filename,'w'))
          self.demux_files_i2.append(open(i2_filename,'w'))

# ------------- MAIN --------------
if __name__ == '__main__':

  # check environment variable
  if not envar_name in os.environ:
    sys.stderr.write("ERROR: Environment variable '{0}' not found. This script expects data files in ${0}/demux.\n".format(envar_name))
    sys.exit(1)

  # check data files
  check_infile(adapter_file)
  check_infile(codetable_samples_file)
  check_infile(codetable_umis_file)

  if len(sys.argv) < 5:
    sys.stdout.write("CAPP-Seq demultiplexing script for FLEX adapters\n")
    sys.stdout.write("Usage: {0} R1.fastq.gz R2.fastq.gz sample2barcode.txt out_dir [I1.fastq.gz I2.fastq.gz])\n".format(sys.argv[0]))
    sys.exit(1)

  r1f = sys.argv[1]
  r2f = sys.argv[2]
  s2bf = sys.argv[3]
  out_dir = sys.argv[4]
  if len(sys.argv) > 5:
    i1f = sys.argv[5]
    if len(sys.argv) < 7:
      error_and_exit("Both index files are required. Only one found.")
    i2f = sys.argv[6]
    use_index_files = True
  else:
    use_index_files = False

  # check file parameters
  check_infile(r1f)
  check_infile(r2f)
  check_infile(s2bf)
  check_outdir(out_dir)
  if use_index_files:
    check_infile(i1f)
    check_infile(i2f)

  # are input files gzip compressed?
  if len(r1f) >= 3 and r1f[-3:] == ".gz" and not force_uncompressed:
    use_gzip = True
    demux_file_tempblate += ".gz"
    unid_file_template += ".gz"
  else:
    use_gzip = False

  # read sample-barcode assignments
  if verbosity > 0:
    sys.stdout.write("Loading sample definitions...\n")
  s2b = Samples2Barcodes(s2bf, out_dir, use_gzip, use_index_files, write_unidentified)
  if verbosity > 0:
    sys.stdout.write("{0} sample barcodes found\n".format(s2b.get_number_of_barcodes(),s2bf))
  if s2b.get_number_of_barcodes() == 0:
    error_and_exit("No sample definitions found.")

  # barcode count output file
  if write_stats_file:
    stats_filename = get_stats_filename(out_dir, r1f)
    stats_file = open(stats_filename,'w')
    out = stats_file
  else:
    out = stdout

  # initialize input fastq files
  r1p = fq.FastQParser(r1f)
  r2p = fq.FastQParser(r2f)
  if use_index_files:
    i1p = fq.FastQParser(i1f)
    i2p = fq.FastQParser(i2f)
    # check first read for correct format
    test_i1 = fq.FastQParser(i1f)
    test_i2 = fq.FastQParser(i2f)
    test_r1 = test_i1.next()
    test_r2 = test_i2.next()
    h1 = test_r1.header
    h2 = test_r2.header
    s1 = test_r1.seq
    s2 = test_r2.seq
    test_i1.close()
    test_i2.close()
    if not h1.split(' ')[1][0] == "1":
      error_and_exit("{0} does not seem to be an I1 file".format(i1f))
    if not h2.split(' ')[1][0] == "2":
      error_and_exit("{0} does not seem to be an I2 file".format(i2f))
    if not re.match(r'[ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}[+][ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}', s1):
      error_and_exit("Unexpected read header format in file {0}".format(i1f))
    if not re.match(r'[ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}[+][ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}', s2):
      error_and_exit("Unexpected read header format in file {0}".format(i2f))

  else:
    # check first read for correct format
    test_r1 = fq.FastQParser(r1f)
    test_r2 = fq.FastQParser(r2f)
    r1 = test_r1.next()
    r2 = test_r2.next()
    h1 = r1.header
    h2 = r2.header
    test_r1.close()
    test_r2.close()
    if not h1.split(' ')[1][0] == "1":
      error_and_exit("{0} does not seem to be an R1 file".format(r1f))
    if not h2.split(' ')[1][0] == "2":
      error_and_exit("{0} does not seem to be an R2 file".format(r2f))
    if not re.search(r':[ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}[+][ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}$', h1):
      print(h1)
      error_and_exit("Unexpected read header format in file {0}".format(r1f))
    if not re.search(r':[ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}[+][ACGTN]{'+str(INDEX_BARCODE_LENGTH)+r'}$', h2):
      error_and_exit("Unexpected read header format in file {0}".format(r2f))

  # auto-detect whether index2 is reverse-complement or not
  if i2_reverse_complement_auto_detect:
    if verbosity > 0:
      sys.stdout.write("Auto-detecting I2 orientation...\n")
    codetable_samples = CodeTable(codetable_samples_file, True) # for auto-detection, load with reverse complement
    match_raw = 0
    match_revcomp = 0
    if use_index_files:
      auto_r2 = fq.FastQParser(i2f)
    else:
      auto_r2 = fq.FastQParser(r2f)
    for idx in range(i2_auto_detect_reads):      
      r2 = auto_r2.next()
      if use_index_files:
        bc2 = r2.header
      else:
        bc2 = r2.header[-INDEX_BARCODE_LENGTH:]
      if codetable_samples.decode(bc2) != CODEWORD_NOT_FOUND:
        match_raw += 1
      if codetable_samples.decode_i2(bc2)!= CODEWORD_NOT_FOUND:
        match_revcomp += 1
    auto_r2.close()
    if verbosity > 1:
      sys.stdout.write("Fwd: {0}, Rev: {1}\n".format(match_raw, match_revcomp))
    if match_raw + match_revcomp > i2_auto_detect_reads / 2:
      # The majority of reads should be decodable, otherwise something is wrong
      if match_raw > match_revcomp:
        if verbosity > 0:
          sys.stdout.write("Based on the first {0} reads, this seems to be MiSeq/NextSeq-style data. Using raw I2 sequences.\n".format(i2_auto_detect_reads))
        i2_reverse_complement = False
      else:
        if verbosity > 0:
          sys.stdout.write("Based on the first {0} reads, this seems to be HiSeq-style data. Using reverse complement I2 sequences.\n".format(i2_auto_detect_reads))
        i2_reverse_complement = True
    else:
      if verbosity > 0:
        sys.stdout.write("Based on the first {0} reads, could not determine whether I2 is reverse complement or not. Using default {1} I2 sequences.\n".format(i2_auto_detect_reads, "reverse complement" if i2_reverse_complement else "raw"))

  # load adapters and codetable
  if verbosity > 0:
    sys.stdout.write("Loading adapter definitions...\n")
  stubs = Stubs(adapter_file)
  if verbosity > 0:
    sys.stdout.write("Loading code tables...\n")
  codetable_samples = CodeTable(codetable_samples_file, i2_reverse_complement)
  codetable_umis = CodeTable(codetable_umis_file)

  # initialize statistics

  c = 0 # number of total fragments processed (for progress output)
  unid_frags = 0 # unidentified fragments
  phix_frags = 0 # suspected phix fragments
  error_corrected_frags = 0 # fragments rescued by error correction (either of I1/I2/UMI1/UMI2)
  # observed sample barcodes I2/I2
  i1_counts = {} # i1 sample barcode
  i2_counts = {} # i2 sample barcode
  swap_counts = {} # observed i1/i2 combinations after decoding, irrespective of UMI validity, for swapping assessment
  # observed UMIs per sample
  demux_frags = {} # number of properly demultiplexed fragments per sample (= reads in demuxed files)
  umi1_counts = {} # umi1 counts per sample
  umi2_counts = {} # umi2 counts per sample
  for bc in s2b.get_barcodes():
    demux_frags[bc] = 0
    umi1_counts[bc] = {}
    umi2_counts[bc] = {}
  correct_first_base_i1 = 0
  correct_first_base_i2 = 0
  correct_pm_i1 = 0
  correct_pm_i2 = 0
  if verbosity > 0:
    sys.stdout.write("Processing reads...\n")

  start_time = datetime.datetime.now()
  out.write(start_time.strftime("%Y-%m-%d %H:%M:%S") + '\n')
  out.write(VERSION + '\n')
  out.flush()

  # process reads (time critical inner loop)
  for read1 in r1p:
    read2 = r2p.next()

    c += 1 # fragment count
    found = False

    # get sample barcode from index read or header
    if use_index_files:
      iread1 = i1p.next()
      iread2 = i2p.next()
      barcode1 = iread1.seq
      barcode2 = iread2.seq
    else:
      # parse index barcodes from r1 header (assuming that i1+i2 is at the end of the header)
      barcodes = read1.header[-(2*INDEX_BARCODE_LENGTH+1):]
      barcode1 = barcodes[:INDEX_BARCODE_LENGTH]
      barcode2 = barcodes[INDEX_BARCODE_LENGTH+1:]

    # now we have the raw index barcode sequences available
    i1_counts[barcode1] = i1_counts.get(barcode1, 0) + 1
    i2_counts[barcode2] = i2_counts.get(barcode2, 0) + 1

    # check for typical phix signature (specific to our design)
    if barcode1 == PHIX_I1 or barcode2 == PHIX_I2:
      phix_frags += 1

    # decode index barcode
    dec1 = codetable_samples.decode(barcode1)
    dec2 = codetable_samples.decode_i2(barcode2)

    # decode UMI
    umi1 = read1.seq[1:7]
    umi2 = read2.seq[1:7]
    udec1 = codetable_umis.decode(umi1)
    udec2 = codetable_umis.decode(umi2)

    if dec1 != CODEWORD_NOT_FOUND and dec2 != CODEWORD_NOT_FOUND:
      new_barcode =  codetable_samples.getCodewordByNum(dec1.codeword_num) + DUAL_INDEX_DELIMITER + codetable_samples.getCodewordByNum(dec2.codeword_num)
      # new_barcode is now error corrected and i2 has been reverse complemented to match i1
      swap_counts[new_barcode] = swap_counts.get(new_barcode,0) + 1
      if s2b.is_sample_barcode(new_barcode):

        # sample has been identified

        # count UMIs per sample
        umi1_counts[new_barcode][umi1] = umi1_counts[new_barcode].get(umi1,0) + 1
        umi2_counts[new_barcode][umi2] = umi2_counts[new_barcode].get(umi2,0) + 1
        if udec1 != CODEWORD_NOT_FOUND and udec2 != CODEWORD_NOT_FOUND:
          stub1 = stubs.getByUmi(codetable_umis.getCodewordByNum(udec1.codeword_num))
          stub2 = stubs.getByUmi(codetable_umis.getCodewordByNum(udec2.codeword_num))
          if stub1 and stub2:
            # everything ok, demultiplex
            demux_frags[new_barcode] += 1
            found = True
            length1 = stub1.length
            length2 = stub2.length
            # first base error?
            if read1.seq[0] == stub1.full_insert[0]:
              correct_first_base_i1 += 1
            if read2.seq[0] == stub2.full_insert[0]:
              correct_first_base_i2 += 1
            # punctuation marks
            pm_errors1 = 0
            pm_errors2 = 0
            if read1.seq[length1-2] != stub1.punct[0]:
              pm_errors1 += 1
            if read1.seq[length1-1] != stub1.punct[1]:
              pm_errors1 += 1
            if read2.seq[length2-2] != stub2.punct[0]:
              pm_errors2 += 1
            if read2.seq[length2-1] != stub2.punct[1]:
              pm_errors2 += 1
            if pm_errors1 == 0:
              correct_pm_i1 += 1
            if pm_errors2 == 0:
              correct_pm_i2 += 1
            # errors corrected?
            if dec1.dist > 0 or dec2.dist > 0 or udec1.dist > 0 or udec2.dist > 0:
              error_corrected_frags += 1
            if demultiplex:
              # trim sequences
              read1.seq = read1.seq[length1:]
              read1.quals = read1.quals[length1:]
              read2.seq = read2.seq[length2:]
              read2.quals = read2.quals[length2:]
              # rewrite header
              add_to_header = fastq_header_demux.format(dec1.codeword_num,dec1.dist,dec1.ns,dec2.codeword_num,dec2.dist,dec2.ns,udec1.codeword_num,udec1.dist,udec1.ns,udec2.codeword_num,udec2.dist,udec2.ns,pm_errors1,pm_errors2)
              read1.header = read1.header.split(' ')[0] + ':' + add_to_header
              read2.header = read2.header.split(' ')[0] + ':' + add_to_header
              # write to fastq
              s2b.write_demuxed_read(read1, new_barcode, 1)
              s2b.write_demuxed_read(read2, new_barcode, 2)
              if use_index_files:
                s2b.write_demuxed_read(iread1, new_barcode, 3)
                s2b.write_demuxed_read(iread2, new_barcode, 4)

    else:
      new_barcode = barcode1 + DUAL_INDEX_DELIMITER + barcode2 # Note: barcode2 is still rev-comp

    if not found:
      unid_frags += 1
      if demultiplex and write_unidentified:
        # trim off maximum stub length so that reads can be mapped
        stub1_seq = read1.seq[:MAX_STUB_LENGTH]
        stub2_seq = read2.seq[:MAX_STUB_LENGTH]
        read1.seq = read1.seq[MAX_STUB_LENGTH:]
        read2.seq = read2.seq[MAX_STUB_LENGTH:]
        read1.quals = read1.quals[MAX_STUB_LENGTH:]
        read2.quals = read2.quals[MAX_STUB_LENGTH:]
        add_to_header = fastq_header_unid.format(barcode1,barcode2,stub1_seq,stub2_seq)
        read1.header = read1.header.split(' ')[0] + ':' + add_to_header
        read2.header = read2.header.split(' ')[0] + ':' + add_to_header
        # write to unidentified file
        s2b.write_unidentified_read(read1, 1)
        s2b.write_unidentified_read(read2, 2)
        if use_index_files:
          s2b.write_unidentified_read(iread1, 3)
          s2b.write_unidentified_read(iread2, 4)

    # Debug: output read info
    if verbosity >= 3:
        sys.stdout.write("{0}\t{1}\t{2}\t{3}\t->\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(barcode1, barcode2, umi1, umi2, s2b.barcode2name(new_barcode) if s2b.is_sample_barcode(new_barcode) else None, dec1, dec2, udec1, udec2, "PASS" if found else "FAIL"))

    # Progress indicator
    if verbosity >= 2:
      if c % 10000000 == 0:
        sys.stdout.write(c)

  # end of inner loop
  out.write(str(c) + '\n')
  sum_demux = sum(demux_frags.values())
  out.write(str(sum_demux) + '\n')
  out.write(str(unid_frags) + '\n')

  sorted_swap_counts = dict_items_by_val(swap_counts) # these are error corrected and i2 reverse complemented

  # print sample barcode stats
  sum_i1 = sum(i1_counts.values())
  sum_i2 = sum(i2_counts.values())
  assert sum_i1 == sum_i2
  sum_all = sum_i1

  # Demultiplexed fragments per sample
  sum_samples = 0
  invalid_umis = 0
  for bc in s2b.get_barcodes():
    bc_count = demux_frags[bc]
    invalid_umis += (swap_counts.get(bc,0)-bc_count)
    assert swap_counts.get(bc,0) == sum(umi1_counts.get(bc, {}).values())
    sum_samples = sum_samples + bc_count
    out.write("{0}\t{1}\t{2}\t".format(bc, bc_count,s2b.barcode2name(bc)))
  out.write('\n')
  assert (sum_all-sum_samples) == unid_frags, "{0} != {1}".format(sum_all-sum_samples, unid_frags)
  out.write(str(correct_first_base_i1) + '\t' + str(correct_first_base_i2) + '\n')
  out.write(str(correct_pm_i1) + '\t' + str(correct_pm_i2) + '\n')

  invalid_index = unid_frags - invalid_umis - phix_frags
  out.write(str(phix_frags) + '\n')
  out.write(str(invalid_index) + '\n')
  out.write(str(invalid_umis) + '\n')
  assert(sum_all == c)

  # check if non-sample barcodes have high counts
  for (bc, count) in sorted_swap_counts[0:s2b.get_number_of_barcodes()]:
    if not s2b.is_sample_barcode(bc):
      i1 = bc[:INDEX_BARCODE_LENGTH]
      i2r = reverse_complement(bc[INDEX_BARCODE_LENGTH+1:])
      out.write("{0}\t{1}\t{2}\t{3}\t".format(bc,count, i1, i2r))
  out.write('\n')

  # swapping
  # expected i1/i2 (from sample2barcode)
  expected_i1_species = {}
  expected_i2_species = {}
  for bc in s2b.get_barcodes():
    i1 = bc[:INDEX_BARCODE_LENGTH]
    i2 = bc[INDEX_BARCODE_LENGTH+1:]
    expected_i1_species[i1] = True
    expected_i2_species[i2] = True

  # counts
  expected_frags = 0
  swapped_frags = 0
  for i1 in expected_i1_species:
    out.write("{0}\t".format(i1))
  out.write("\n")
  for i2 in expected_i2_species:
    out.write("{0}\t".format(i2))
  out.write("\n")
  for i1 in expected_i1_species:
    for i2 in expected_i2_species:
      barcode = i1 + DUAL_INDEX_DELIMITER + i2
      count = swap_counts.get(barcode,0)
      if s2b.is_sample_barcode(barcode):
        expected_frags += count
      else:
        swapped_frags += count
      out.write("{0}\t".format(count))
    out.write(";")
  out.write('\n')

  total_frags = expected_frags + swapped_frags
  out.write(str(total_frags) + '\n')
  out.write(str(expected_frags) + '\n')
  out.write(str(swapped_frags) + '\n')

  # UMIs
  # umi1
  umi1_all = {}
  sum_u1 = 0
  for bc in s2b.get_barcodes():
    for umi, count in umi1_counts.get(bc, {}).items():
      umi1_all[umi] = umi1_all.get(umi, 0) + count
      sum_u1 += count
  sorted_umi1_all = dict_items_by_val(umi1_all)
  non_decodable1 = 0
  perfect_match1 = 0
  corrected1 = 0
  n_corrected1 = 0
  out.write(str(sum_u1) + '\n')
  for idx, (umi, count) in enumerate(sorted_umi1_all):
    udec = codetable_umis.decode(umi)
    msg = "U%02d+%d+%d\t" %(udec.codeword_num,udec.dist,udec.ns) if udec != CODEWORD_NOT_FOUND else "\t"
    if udec != CODEWORD_NOT_FOUND and stubs.getByUmi(codetable_umis.getCodewordByNum(udec.codeword_num)):
      msg += stubs.getByUmi(codetable_umis.getCodewordByNum(udec.codeword_num)).name
      if udec.dist == 0:
        perfect_match1 += count
      else:
        corrected1 += count
        if udec.ns > 0:
          n_corrected1 += count
    else:
      non_decodable1 += count
    out.write("{0}\t{1}\t{2}\t".format(umi, count, msg))
  out.write("\n")
  out.write(str(perfect_match1) + '\n')
  out.write(str(corrected1) + '\n')
  out.write(str(n_corrected1) + '\n')
  out.write(str(non_decodable1) + '\n')

  if non_decodable1 > sum_u1 / 2:
    out.write("\nWARNING: The majority of UMI sequences in read 1 could not be identified. This suggests problems with the data or UMI definitions.\n")

  # umi2
  umi2_all = {}
  sum_u2 = 0
  for bc in s2b.get_barcodes():
    for umi, count in umi2_counts.get(bc, {}).items():
      umi2_all[umi] = umi2_all.get(umi, 0) + count
      sum_u2 += count
  sorted_umi2_all = dict_items_by_val(umi2_all)
  out.write(str(sum_u2) + '\n')
  non_decodable2 = 0
  perfect_match2 = 0
  corrected2 = 0
  n_corrected2 = 0
  for idx, (umi, count) in enumerate(sorted_umi2_all):
    udec = codetable_umis.decode(umi)
    msg = "V%02d+%d+%d\t" %(udec.codeword_num,udec.dist,udec.ns) if udec != CODEWORD_NOT_FOUND else "\t"
    if udec != CODEWORD_NOT_FOUND and stubs.getByUmi(codetable_umis.getCodewordByNum(udec.codeword_num)):
      msg += stubs.getByUmi(codetable_umis.getCodewordByNum(udec.codeword_num)).name
      if udec.dist == 0:
        perfect_match2 += count
      else:
        corrected2 += count
        if udec.ns > 0:
          n_corrected2 += count
    else:
      non_decodable2 += count
    out.write("{0}\t{1}\t{2}\t".format(umi, count, msg))
  out.write("\n")
  out.write(str(perfect_match2) + '\n')
  out.write(str(corrected2) + '\n')
  out.write(str(n_corrected2) + '\n')
  out.write(str(non_decodable2) + '\n')
  assert sum_u1 == sum_u2

  out.write(str(sum_all) + '\n')
  sorted_i1_counts = dict_items_by_val(i1_counts) # these are raw counts without any filters
  non_decodable = 0
  other_barcode = 0
  perfect_match = 0
  corrected = 0
  n_corrected = 0
  phix = 0
  for idx, (seq, count) in enumerate(sorted_i1_counts):
    dec = codetable_samples.decode(seq)
    if dec == CODEWORD_NOT_FOUND:
      non_decodable += count
      if seq == PHIX_I1:
        msg = "PhiX\t"
        phix += count
      else:
        msg = "\t"
    else:
      msg = "I%02d+%d+%d\t" %(dec.codeword_num,dec.dist,dec.ns)
      codeword = codetable_samples.getCodewordByNum(dec.codeword_num)
      bc = DUAL_INDEX_DELIMITER.join([codeword]*2) # try if this is a barcode with same I1/I2
      if s2b.is_sample_barcode(bc):
        sample_name = s2b.barcode2name(bc)
        msg += sample_name
        if dec.dist == 0:
          perfect_match += count
        else:
          corrected += count
          if dec.ns > 0:
            n_corrected += count
      else:
        other_barcode += count
    out.write("{0}\t{1}\t{2}\t".format(seq, count, msg))
  out.write('\n')
  out.write(str(perfect_match) + '\n')
  out.write(str(corrected) + '\n')
  out.write(str(n_corrected) + '\n')
  out.write(str(other_barcode) + '\n')
  out.write(str(non_decodable) + '\n')
  out.write(str(phix) + '\n')

  sorted_i2_counts = dict_items_by_val(i2_counts) # TODO: should we reverse complement these?
  non_decodable = 0
  other_barcode = 0
  perfect_match = 0
  corrected = 0
  n_corrected = 0
  phix = 0
  for idx, (raw_seq, count) in enumerate(sorted_i2_counts):
    seq = reverse_complement(raw_seq)
    dec = codetable_samples.decode(seq)
    if dec == CODEWORD_NOT_FOUND:
      non_decodable += count
      if seq == reverse_complement(PHIX_I2):
        msg = "PhiX\t"
        phix += count
      else:
        msg = "\t"
    else:
      msg = "J%02d+%d+%d\t" %(dec.codeword_num,dec.dist,dec.ns)
      codeword = codetable_samples.getCodewordByNum(dec.codeword_num)
      bc = DUAL_INDEX_DELIMITER.join([codeword]*2) # try if this is a barcode with same I1/I2
      if s2b.is_sample_barcode(bc):
        sample_name = s2b.barcode2name(bc)
        msg += sample_name
        if dec.dist == 0:
          perfect_match += count
        else:
          corrected += count
          if dec.ns > 0:
            n_corrected += count
      else:
        other_barcode += count
    out.write("{0}\t{1}\t{2}\t".format(seq, count, msg))
  out.write("\n")
  out.write(str(perfect_match) + '\n')
  out.write(str(corrected) + '\n')
  out.write(str(n_corrected) + '\n')
  out.write(str(other_barcode) + '\n')
  out.write(str(non_decodable) + '\n')
  out.write(str(phix) + '\n')
  other = unid_frags-phix_frags-invalid_umis-swapped_frags
  out.write(str(error_corrected_frags) + '\n')
  out.write(str(other) + '\n')
  stop_time = datetime.datetime.now()
  out.write(stop_time.strftime("%Y-%m-%d %H:%M:%S"))

  # clean up
  if write_stats_file:
    stats_file.close()
  if demultiplex:
    s2b.close()

  if verbosity > 0:
    sys.stdout.write("Statistics written to '{0}'\n".format(stats_filename))
    sys.stdout.write("Done.\n")
