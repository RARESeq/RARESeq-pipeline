input_file=$1
exons_bed=$2
output_dir=$3
mkdir -p $output_dir

file_name=$(basename $input_file)
tag=${file_name/.bam/}

samtools view $input_file  |\
awk -v OFS='\t' '{
		if($6 ~ /^([0-9]+[NM])+$/)
		{
			#print "OK", $6;
			split($6, M_splits, "M");
			size = 0;
			for(spl in M_splits)
			{
				if(length(M_splits[spl]) > 0)
				{
					#print M_splits[spl];
					if(index(M_splits[spl], "N"))
					{
						split(M_splits[spl], N_splits, "N");
						size += N_splits[2];
					}else{
						size += M_splits[spl];
					}
					#print size;
				}
			}
			if(size >= 100)
			{
				print $3, $4, $4, $1"__"$2"__"$6;
			}
		}
	}'|\
	bedtools intersect -a stdin -b $exons_bed -wo >$output_dir/$tag.bed

cat $output_dir/$tag.bed | awk 'BEGIN{OFS="\t";}
	{
	split($4, res, "__");
	if($1 == "chrM")
	{
		next;
	}
	if(index(res[3], "N"))
	{
		spl ++;
	}else{
		not_spl ++;
	}
	}
	END{
		print "Split", spl;
		print "NotSplit", not_spl;
	}' >$output_dir/$tag.contamination_stats.txt

rm -f $output_dir/$tag.bed
