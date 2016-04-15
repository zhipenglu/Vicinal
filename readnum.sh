#!/bin/bash
# batch execution of samtools view 
# Author: Zhipeng Lu
# Contact: zhipengluchina@gmail.com

if [ $# -lt 3 ]
then
   echo 
   echo "Not enough arguments."
   echo "Usage: $0 indexed_bam_file ncRNA_list output_file" 
   echo "ncRNA_list file format example: 'geneID geneName chr2L:865365-865493 note'"
   echo
else
   while read RNA_ID RNA_name RNA_pos Other
   do 
   RNA_record="$RNA_ID\t$RNA_name\t$RNA_pos\t$(samtools view $1 $RNA_pos | wc -l)\t$other"
   echo -e $RNA_record >> $3
   done < $2
fi

