import sys
import os
import re
import pysam

input_file = sys.argv[1]
fa_index = sys.argv[2]
output_dir = sys.argv[3] 
name_root = os.path.basename(input_file)

try:
	os.makedirs(output_dir)
except Exception as e:
	pass

chr_lengths = {}
for line in open(fa_index, "r"):
	col = line.split("\t")
	chr_lengths[col[0]] = int(col[1])

d = {}
reg = re.compile("(\d+?)([A-Z])+?")

def read_removeclipping(read):
    if read.is_unmapped:
        return 1

    newcigar = []
    clip_5 = 0
    clip_3 = 0

    changed = False
    inseq = False
    for op, length in read.cigar:
        if op == 5:  # H
            changed = True
        elif op == 4:  # S
            changed = True
            if not inseq:
                clip_5 = length
            else:
                clip_3 = length
        else:
            inseq = True
            newcigar.append((op, length))

    if not changed:
        return 0

    read.cigar = newcigar
    orig_length = len(read.seq)

    s = read.seq
    q = read.qual

    if clip_3:
        read.seq = s[clip_5:-clip_3]
        if q:
            read.qual = q[clip_5:-clip_3]
    else:
        read.seq = s[clip_5:]
        if q:
            read.qual = q[clip_5:]

    newtags = []
    if clip_5:
        newtags.append(('ZA', clip_5))
    if clip_3:
        newtags.append(('ZB', clip_3))

    newtags.append(('ZC', float(clip_5 + clip_3) / orig_length))

    read.tags = read.tags + newtags

    return 2

def parse_cigar(cigar, start):
    matches = reg.findall(cigar)
   	#print cigar
    #print matches
    size = 0
    for el in matches:
       if el[1] == "M" or el[1] == "N" or el[1] == "D":
         start += int(el[0])
       elif el[1] == "S" or el[1] == "I":
         pass
       else:
         print("Error: Cigar string not in the expected format: " + cigar + "\n")
         return -1
    return start


def process_soft_clipping(read):
	to_clip = False
	res = re.findall(reg, read.cigarstring)
	 
	if res[-1][1] == "S":
		length = int(res[-1][0])
		chr_end = chr_lengths[read.reference_name]
		read_end = parse_cigar(read.cigarstring, read.reference_start)
		#print read
		#print chr_end, read_end, length
		if length <= 2 and read_end + length < chr_end: 
			res[-2] = (str(int(res[-2][0]) + length), res[-2][1])
			res = res[:-1]
		else:
			to_clip = True

	if res[0][1] == "S":
		length = int(res[0][0])
		if length <= 2 and (read.reference_start - length) >= 1:
			res[1] = (str(int(res[1][0]) + length), res[1][1])
			res = res[1:]
			read.reference_start = read.reference_start - length
		else:
			to_clip = True

	read.cigarstring = "".join(["".join(x) for x in res])

	if to_clip == True:
		read_removeclipping(read) 
		
valid_flags = {99 : True, 147 : True, 83 : True, 163 : True}

cnt = 0
samfile = pysam.AlignmentFile(input_file, "rb")

output_samfile = pysam.AlignmentFile(os.path.join(output_dir, name_root.replace(".bam", "_tmp.bam")), "wb", template = samfile)
prev_read = None
error_messages = 0
for read in samfile.fetch(until_eof = True):
	cnt += 1
	#if cnt % 100000 == 0:
	#	print("(1) Processed %s"%cnt)
		
	if not read.flag in valid_flags:
		#print "skipping read as it's not properly mapped", read.flag, read
		continue

	process_soft_clipping(read)
	
	if prev_read is None:
		prev_read = read
		continue
	else:
		if prev_read.qname != read.qname:
			error_messages += 1
			if error_messages < 100:
				print("[ERROR] Reads don't seem to be paired:" + prev_read.qname + " " + read.qname)
		else:
			
			read.setTag("MC", prev_read.cigarstring, value_type="Z")
			prev_read.setTag("MC", read.cigarstring, value_type="Z")
			
			try:
				read_tags = read.get_tag("ZA")
				prev_read.setTag("ZD", read_tags, value_type="i")
			except:
				pass
			try:
				read_tags = read.get_tag("ZB")
				prev_read.setTag("ZE", read_tags, value_type="i")
			except:
				pass
			try:
				read_tags = prev_read.get_tag("ZA")
				read.setTag("ZD", read_tags, value_type="i")
			except:
				pass
			try:
				read_tags = prev_read.get_tag("ZB")
				read.setTag("ZE", read_tags, value_type="i")
			except:
				pass
			
			read.next_reference_start = prev_read.reference_start
			prev_read.next_reference_start = read.reference_start
			output_samfile.write(prev_read)
			output_samfile.write(read)
			prev_read = None
			read = None
samfile.close()
output_samfile.close()

os.rename(os.path.join(output_dir, name_root.replace(".bam", "_tmp.bam")), os.path.join(output_dir, name_root))

