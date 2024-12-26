import sys
import os
import re
import pysam
import gzip

input_file = sys.argv[1]
output_dir = sys.argv[2]
#name_root = os.path.basename(input_file).split(".Aligned")[0]
name_root = os.path.basename(input_file).split(".bam")[0]
fq1 = input_file.split("__")[0]
fq1 = fq1.split(".sorted")[0]

fastp = fq1.replace("QC_prededupe", "fastp/good")
fastp1 = fastp + "_R1.fastq.gz"
fastp2 = fastp + "_R2.fastq.gz"

demux = fq1.replace("exp_analysis/QC_prededupe", "demultiplexed")
demux1 = demux + "_R1.fastq"
demux2 = demux + "_R2.fastq"

#os.system("zcat " + fastp1 + " | head ")

try:
	os.makedirs(output_dir)
except Exception as e:
	pass

def parse_cigar(cigar, start):
	#print cigar, start
	splices = re.findall("([0-9]+)([SMN]+)", cigar) 
	segments =[]
	end = start
	size = 0
	read_size = 0
	for splice in splices:
		if splice[1] == "M":
			size += int(splice[0])
			read_size += int(splice[0])
			segments.append((end, end + int(splice[0])))
		if splice[1] == "S":
			read_size += int(splice[0])
			continue
		end += int(splice[0])
	
	#print splices
	#print start, end, size, segments
	return {"start" : start, "end" : end, "size" : size, "read_size" : read_size, "segments" : segments}

samfile = pysam.AlignmentFile(input_file, "rb")

cnt = 0
prev_chr = -1
d = {}
overlapping = 0
non_overlapping = 0
inconsistent = 0
counts = {}

to_filt = {}
#f = open(os.path.join(output_dir, "%s.size"%name_root), "w")
for read in samfile.fetch(until_eof = True):
	cnt += 1
	if cnt % 1000000 == 0: 
		#print("Processed %s"%cnt)
		break
		
	if not (re.match("^([0-9]+[SMN])+$", read.cigarstring)):
		continue
		
	if not read.flag in [83, 163, 99, 147]:
		continue
		
	if read.reference_id != prev_chr:
		#print("Change from chr %s to chr %s"%(prev_chr, read.reference_id))
		d = {}
		prev_chr = read.reference_id 
	
	if not read.qname in d:
		d[read.qname] = read
	else:
		mate = d[read.qname]
		details = [parse_cigar(read.cigarstring, read.reference_start), parse_cigar(mate.cigarstring, mate.reference_start)]

		first = -1
		second = -1
		if details[0]["start"] <= details[1]["start"] <= details[0]["end"]:
			first = 0
			second = 1
			
		if details[1]["start"] <= details[0]["start"] <= details[1]["end"]:
			first = 1
			second = 0
		
		if first == -1:
			non_overlapping += 1
			#f.write("%s\t%s\n"%(read.qname, "-1"))
			continue
		
		overlapping += 1
		
		#if details[1]["start"] <= details[0]["start"] <= details[0]["end"] < details[1]["end"] or\
		#   details[0]["start"] <= details[1]["start"] <= details[1]["end"] < details[0]["end"] or\
		#   details[1]["start"] < details[0]["start"] <= details[0]["end"] <= details[1]["end"] or\
		#   details[0]["start"] < details[1]["start"] <= details[1]["end"] <= details[0]["end"]:
		#	inconsistent += 1
		#	print details
		#	#print "Another weird case!!!"
		#	#print details
		#	continue
		
		first_idx = 0
		second_idx = 0
		overlap = 0
		
		while details[first]["segments"][first_idx][1] < details[second]["segments"][second_idx][0]:
			overlap += details[first]["segments"][first_idx][1] - details[first]["segments"][first_idx][0]
			first_idx += 1
		
		ok = False
		#print details
		while first_idx < len(details[first]["segments"]) and second_idx < len(details[second]["segments"]):
			if details[first]["segments"][first_idx][0] <= details[second]["segments"][second_idx][0] <= details[first]["segments"][first_idx][1]:
				ok = True
			if details[second]["segments"][second_idx][0] <= details[first]["segments"][first_idx][0] <= details[second]["segments"][second_idx][1]:
				ok = True

			if not ok: 
				#print("breaking")
				break

			overlap += max(details[first]["segments"][first_idx][1], details[second]["segments"][second_idx][1]) - min(details[first]["segments"][first_idx][0], details[second]["segments"][second_idx][0])
			first_idx += 1
			second_idx += 1
		
		if first_idx < len(details[first]["segments"]):
			inconsistent += 1
			continue 

		while second_idx < len(details[second]["segments"]):
			overlap += details[second]["segments"][second_idx][1] - details[second]["segments"][second_idx][0]
			second_idx += 1
			
		#molecule_size = details[first]["size"] + details[second]["size"] - overlap 
		molecule_size = overlap

		if molecule_size in counts:
			counts[molecule_size] += 1
		else:
			counts[molecule_size] = 1
		del d[read.qname]
samfile.close()
#f.close()
f = open(os.path.join(output_dir, "%s.molecule_size.freq"%name_root), "w")
f.write("Molecule_size\tCount\n")
for k, v in counts.items():
	f.write("%s\t%s\n"%(k, v))
f.close()

f = open(os.path.join(output_dir, "%s.molecule_size.summary"%name_root), "w")
f.write("Overlapping\tNon_Overlapping\tInconsistent\n")
f.write("%s\t%s\t%s\n"%(overlapping, non_overlapping, inconsistent))
f.close()

exit(0)

def dump_fq(fp, reads, name = "", gz = True):
	g = open(os.path.basename(fp +  name), "w")
	if gz:
		f = gzip.open(fp, 'rb')
	else:
		f = open(fp, 'rb')
	printing = False
	start = 0
	cnt = 0
	for line in f:
		cnt += 1
		if cnt % 1000000  == 0:
			print(cnt)
		if not printing:
			s = line.strip().split(" ")
			if s[0] in reads:
				printing = True
			else:
				#print s
				#break
				continue
			#print(line.strip())
			g.write(line)
			start = 1
		else:
			if start <= 4:
				#print(line.strip())
				g.write(line)
				start += 1
				if start == 4:
					printing = False
	g.close()

dump_fq(fastp1, to_filt, ".fastp", True)
dump_fq(fastp2, to_filt, ".fastp", True)
dump_fq(demux1, to_filt, ".demux", False)
dump_fq(demux2, to_filt, ".demux", False)
