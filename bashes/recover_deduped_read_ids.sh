input_file=$1
samtools view $input_file | awk '{
   res=$1
   n=split(res,resArr,":")
   #print n;
   s=resArr[1];
   for(i = 4; i <=n; i++)
   {	
   	s=s":"resArr[i];
   }
   print(s);
}' >$input_file.read_ids