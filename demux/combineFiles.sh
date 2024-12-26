#!/bin/bash
echo "Combining ${2}"
cat ${1}/scratch/R*${2} > ${1}/demultiplexed/${2}
rm ${1}/scratch/R*${2}

#### Uncomment to zip
# echo "Zipping ${2}"
# gzip ${1}/demultiplexed/${2}
