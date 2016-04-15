# Vicinal
Detection of ncRNA ends using chimeric reads from RNA-seq data

by Zhipeng Lu, UNC at Chapel Hill, NC 27517, USA
2014-01-14

Copyright 2014, Zhipeng Lu, University of North Carolina at Chapel Hill
Contact information: zhipengluchina@gmail.com

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or 
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Citations:
Lu, Z. and Matera A.G. (2014). Vicinal: a method for the determination of ncRNA ends using chimeric reads from RNA-seq experiments. Nucleic Acid Research 42 (9), e79-e79
Lu, Z., Guan, X., Schmidt, C.A. and Matera, A.G. (2014) RIP-seq analysis of eukaryotic Sm proteins identifies three major categories of Sm-containing ribonucleoproteins. Genome Biology, 15, R7.

The Vicinal package contains three python script files: samsoftfilter.py, Vicinal_1.0.py and Vicinal_2.0.py, and one shell script readnum.sh. The samsoftfilter.py script filters Bowtie2 mapped reads in SAM files, and select partially mapped reads. The Vicinal scripts then maps the unmapped fragments to the vicinity of the mapped fragments. The python scripts require Python 2.7 and above. Some of the text processing steps used functions from samtools, bedtools and unix utilities. The readnum.sh extracts and counts the number of mapped chimeric reads in genomic intervals for ncRNAs (or other genomic features: e.g. genes) via batch execution of samtools.

The following is the complete pipeline for the analysis of RNA chimeric reads to determine the accurate ncRNA ends. Basic usage of each program in the pipeline is presented. However if you want to know the details of the programs other than Vicinal, please refer the respective manuals. Understanding the pipeline and using it requires basic knowledge of the command line utilities. Good luck with your analysis.

Note: Sometimes the --rad option has conflicts with the ref when it is close the end. So make sure there is enough distance at least the distance used by '--rad' from the end when analyzing genes in isolation

Overall flow of data: file.fastq --> file_local.sam --> file_soft.sam --> file_1.wig + file_2.wig + file_chim.sam



1. Bowtie1 mapping [Optional, for comparison with bowtie2 local mapping]
$ bowtie -v 2 <index> -q Lu001.fastq -S Lu001.sam 
$ bowtie -v 2 <index> -q Lu002.fastq -S Lu002.sam 
......
Combine all the sam files and make a bedgraph file [optional, for display]

2. Bowtie2 Mapping 
$ bowtie2 --sensitive-local --no-unal -x <index> -U <Lu001.fastq> -S Lu001_local.sam
......
Combine all the sam files and make a bedgraph file [optional, for comparison with bowtie2 local mapping]

3. Filtering (header lines are discarded in the filtering, and should be added back after combining several files)
$ python samsoftfilter.py [-s <softmin>] Lu001_local.sam Lu001_soft.sam
$ python samsoftfilter.py [-s <softmin>] Lu002_local.sam Lu002_soft.sam
......
Combining filtered reads
Extract header from sam file then concatenate the remaining, finally add back one copy of header. Note that it is faster to run Vicinal on each filtered sam file and then combine the results, if you can run them on computing clusters.

4. Mapping the second parts
$ python Vicinal_1.0.py --len <READLEN> [--rad RADIUS] refdir softclipsam outname
Note: Vicinal_1.0.py is more suitable for smaller genomes (like fly or nematode), while Vicinal_2.0.py is more suitable for bigger genomes (like mouse or human), because Vicinal_1.0.py initializes an empty dictionary to enable fast storage of the mapping data and it takes large amount of memory, approximately 200 times the size of the biggest chromosome (6G memory for a 30M fly chromosome, 80G memory for a 400M human chromosome)
(output files: chimeric.sam, chimeric1.wig, chimeric2.wig)

5. Converting wig files to bedgraph files using the Kent utilities (bedgraph files are smaller and are easy to work with using bedtools)
$ ./wigToBigWig in.wig chrom.sizes out.bw
$ ./bigWigToBedGraph in.bigWig out.bedGraph
(Zipping bedGraph files further saves space, and saves time when uploading to genome browser)

6. Combining bedgraph files [also used for reads with different lengths]
$ ./unionBedGraphs -i a.bg b.bg c.bg ... 
Then add up from the 4th to the last columns and add an appropriate header. Example:
$ awk '{sum=0; for (i=4; i <=NF; i++) {sum+=$i} print $1, $2, $3, sum}' < input.bg > output.bg

7. Extracting chimeric reads
$ samtools view -bS -o chimeric.bam chimeric.sam
$ samtools sort chimeric.bam chimeric.sorted
$ samtools index chimeric.sorted.bam
$ samtools view chimeric.sorted.bam chr2L:3,046,746-3,046,904 | cut -f10 | sort | uniq -c 
Then align the reads manually to identify the boundaries of the RNA

8. Generate a list of ncRNAs with numbers of chimeric reads
(need a reference list of ncRNAs and their genomic coordinates)
$ ./readnum.sh indexed_bam_file ncRNA_list output_file_name
Here is an example record for the ncRNA_list file:
CR32864-RA      7SLRNA:CR32864  chr2L:865365-865493	FBgn0000003

9. Secondary structure prediction. You can try Vienna RNA package webserver, UNAfold, and then draw secondary structure manually or using VARNA (Darty et al. 2009, http://varna.lri.fr/). VARNA allows manual manipulation of the overall structure and detailed labeling, and is very easy to use.


