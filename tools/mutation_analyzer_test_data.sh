#!/bin/bash

# 1. check whether samtools, Trimmomatic and bwa are installed 
# --> I already checked and inserted the right paths and names
# 2. download the following files and name them accordong to this:
# - choose a basename for all files e.g.: "ACH_old1_BAT_FS1_RL36"
# - the file with the DCS Mutations should be called: "FreeBayes_ACH_old1_BAT_FS1_RL36.vcf"
# - the BAM file of the DCS reads: "DCS_ACH_old1_BAT_FS1_RL36.bam"
# - the BAM file of the SSCS reads: "SSCS_ACH_old1_BAT_FS1_RL36.bam"
# - the file with the output of align families: "Aligned_Families_ACH_old1_BAT_FS1_RL36.tabular"
# 3. adapt the paths and names that are indicated by comments below and copy all variables to the terminal.
# Then copy all the statements from 4.-8. one after the other and execute them. 
# The order of the statements is important.

# change this path to the folder where your data is
path="/home/admin-isse/Documents/VariantAnalyzerGalaxy/tools/test-data"
# change this path so that it points to your reference fasta 
ref="reference.fasta"
# make sure that it is indexed or perform:
# bwa index $ref

# adapt basename to the files you want to analyze
basename="test"

# Threshold for displaying mutations in the final output file. 
# Only mutations occuring less than thresh times are displayed. 
# Default of 0 displays all.
thresh=0

# Threshold for Phred score.
# Only reads with Phred score higher than this threshold are used.
# Default is 20.
phred=20

outfile1=$path/"FreeBayes_"$basename".vcf"
file2=$path/"DCS_"$basename".bam"
file3=$path/"Aligned_Families_"$basename".tabular"

outfile=$path/"Interesting_Reads_"$basename".fastq"
pickle_file=$path/"tag_count_dict_"$basename".json"

file=$path"/Interesting_Reads_"$basename""

file2s=$path/"SSCS_"$basename".bam"
sscs_counts_pickle=$path/"SSCS_counts_"$basename".json"

outfile2=$path/"Variant_Analyzer_summary_"$basename".xlsx"
outfile3=$path/"Variant_Analyzer_summary_"$basename".csv"
outfile4=$path/"Variant_Analyzer_allele_frequencies_"$basename".xlsx"
outfile5=$path/"Variant_Analyzer_tiers_"$basename".xlsx"

# if you did not download the .bai files:
samtools index $file2
samtools index $file2s

# 5. Creates fastq file of reads of tags with mutation.
python2.7 mut2read.py --mutFile $outfile1 --bamFile $file2 --familiesFile $file3 --outputFastq $outfile --outputJson $pickle_file

# 6. Trim reads of fastq files.
# I already changed this path for you
java -jar /home/admin-isse/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 10 -trimlog $file.trim.log $outfile $file.trim.fastq LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36

# 7. Align trimmed reads to reference.
/home/admin-isse/bwa-0.7.17/bwa mem -t 10 $ref $file.trim.fastq > $file.trim.sam
samtools sort $file.trim.sam > $file.trim.bam
samtools index $file.trim.bam

# 8. Calculates statistics about number of ab/ba for mutation and wildtype for each mutation.
python2.7 mut2sscs.py --mutFile $outfile1 --bamFile $file2s --outputJson $sscs_counts_pickle

# 9. Looks for reads with mutation at known positions and calculates frequencies and stats.
python2.7 read2mut.py --mutFile $outfile1 --bamFile $file.trim.bam --inputJson $pickle_file --sscsJson $sscs_counts_pickle --outputFile $outfile2 --outputFile_csv $outfile3 --outputFile2 $outfile4 --outputFile3 $outfile5 --thresh $thresh --phred $phred --trim 10 --chimera_correction --softclipping_dist 15 --reads_threshold 1.0

