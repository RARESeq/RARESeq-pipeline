#!/bin/bash
echo
echo "Applying fastp to $1 on ${15} threads"

if [ -z $DOCKERRUN ]
then
fpath=${pipeline_dir}
else
fpath=""
fi

#pigz -d $1 $2
f1=$1
f2=$2

f1_basename=$(basename $f1)
f2_basename=$(basename $f2)

mkdir -p ${14}

${fpath}fastp -i ${f1/.gz/} -I ${f2/.gz/} -o ${14}/${f1_basename/.gz/} -O ${14}/${f2_basename/.gz/} --trim_front1=${9} --trim_tail1=${10} --trim_front2=${11} --trim_tail2=${12} \
	-c --overlap_diff_limit=$5 --length_required=$6 --overlap_len_require=$4 --qualified_quality_phred=$7 --unqualified_percent_limit=$8 --html=$3.html --json=$3.json \
	--failed_out=${13} --thread=${15}
 
pigz -f -p ${15} ${14}/${f1_basename/.gz/}
pigz -f -p ${15} ${14}/${f2_basename/.gz/}
