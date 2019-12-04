# Variant Analyzer

This tools-set allows deeper insights into duplex sequencing variant calls. The tools can be used within the [Galaxy platform](http://usegalaxy.org) under the section Du Novo, on your local Galaxy installation from the [toolshed](https://toolshed.g2.bx.psu.edu/view/iuc/variant_analyzer) or as stand-alone Python scrips. They can be used as stand-alone tools or as part of the [Du Novo Analysis Pipeline](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1039-4).

## Dependencies
All tools were built for Python 2.7, but can be easily adapted to newer versions.

## Usage
A detailed description of all tools can be found on [Galaxy](http://usegalaxy.org) with all parameters, input and output files.

### DCS mutations to tags/reads
Takes a tabular file with mutations, a BAM file of aligned DCS reads, and a 
tabular file with aligned families as input and prints all tags of reads that 
carry a mutation to a user specified output file and creates a fastq file of 
reads of tags with a mutation.

**Input** 

**Dataset 1 (--mutFile):** Tabular file with duplex consesus sequence (DCS) mutations as 
generated by the **Variant Annotator** tool.

**Dataset 2 (--bamFile):** BAM file of aligned DCS reads. This file can be obtained by the 
tool [Map with BWA-MEM](https://arxiv.org/abs/1303.3997).

**Dataset 3 (--familiesFile):** Tabular file with reads as produced by the 
**Du Novo: Align families** tool of the [Du Novo Analysis Pipeline](https://doi.org/10.1186/s13059-016-1039-4).

**Output**

The output is a json file containing dictonaries of the tags of reads containing mutations 
in the DCS and a fastq file of all reads of these tags.

`$ python2 mut2read.py --mutFile DCS_mutations.tabular --bamFile DCS.bam         --familiesFile aligned_families.tabular --outputFastq interesting_reads.fastq --outputJson tag_count_dict.json`

### DCS mutations to SSCS stats
Takes a tabular file with DCS mutations and a BAM file of aligned SSCS reads 
as input and writes statistics about tags of reads that carry a mutation in the 
SSCS at the same position a mutation is called in the DCS to a user specified output file..

**Input** 

**Dataset 1 (--mutFile):** Tabular file with duplex consesus sequence (DCS) mutations as 
generated by the **Variant Annotator** tool.

**Dataset 2 (--bamFile):** BAM file of aligned single stranded consensus sequence (SSCS) 
reads. This file can be obtained by the tool [Map with BWA-MEM](https://arxiv.org/abs/1303.3997).

**Output**

The output is a json file containing dictonaries with stats of tags that carry a mutation in the SSCS 
at the same position a mutation is called in the DCS.

`$ python2 mut2sscs.py --mutFile DCS_mutations.tabular --bamFile SSCS.bam --outputJson SSCS_counts.json`

### Call specific mutations in reads
Takes a tabular file with mutations, a BAM file of aligned raw reads, and JSON files 
created by the tools **DCS mutations to tags/reads** and **DCS mutations to SSCS stats** 
as input and calculates frequencies and stats for DCS mutations based on information 
from the raw reads.

**Input** 

**Dataset 1 (--mutFile):** Tabular file with duplex consesus sequence (DCS) mutations as 
generated by the **Variant Annotator** tool.

**Dataset 2 (--bamFile):** BAM file of aligned raw reads. This file can be obtained by the 
tool [Map with BWA-MEM](https://arxiv.org/abs/1303.3997).

**Dataset 3 (--inputJson):** JSON file generated by the **DCS mutations to tags/reads** tool 
containing dictonaries of the tags of reads containing mutations 
in the DCS.

**Dataset 4 (--sscsJson):** JSON file generated by the **DCS mutations to SSCS stats** tool 
stats of tags that carry a mutation in the SSCS at the same position a mutation 
is called in the DCS.

**Output**

The output is an XLSX file containing frequencies stats for DCS mutations based 
on information from the raw reads. In addition to that a tier based 
classification is provided based on the amout of support for a true variant call.

`$ python2 read2mut.py --mutFile DCS_mutations.tabular --bamFile interesting_reads.bam --inputJson tag_count_dict.json --sscsJson SSCS_counts.json --thresh 0 --phred 20 --trim 10 --outputFile variant_analyzer_output.xlsx`
