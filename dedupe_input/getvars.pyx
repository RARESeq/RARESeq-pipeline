import re

reg = re.compile("([0-9]+?)([NMDIS])")
def parse_cigar(cigar, start):
    matches = reg.findall(cigar)
    print cigar
    print matches
    size = 0
    for el in matches:
       if el[1] == "M" or el[1] == "N":
         start += int(el[0])
       elif el[1] == "S":
         pass
       else:
         print "Error: Cigar string not in the expected format: $cigar\n"
         return -1
    return start

