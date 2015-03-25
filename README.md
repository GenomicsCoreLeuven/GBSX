GBSX: a toolkit for experimental design and demultiplexing genotyping by sequencing experiments
====

## Overview

Genotyping by Sequencing is an emerging technology for cost effective variant 
discovery and genotyping. However, current analysis tools do not fulfill all 
experimental design and analysis needs.  
  
GBSX is a package of tools to first aid in experimental design, including choice 
of enzymes and barcode design. Secondly, it provides a first analysis step to 
demultiplex samples using in-line barcodes, providing fastq files that can 
easily be plugged into existing variant analysis pipelines.

## Download

The perl script for in silico digests and the compiled program for all other 
analyses can be found in the releases directory. The latest directory has the 
latest version. However, previous versions are still available. The complete 
source code can be found in the src directory. Example data and results for the 
tool can be found in the example directory.

##Licence

All parts of this tool is licenced under GPLv3.  
A copy of this licence is included under LICENSE.

## Contact

[Genomics Core](http://www.genomicscore.be "Genomics Core website")  
Center for Human Genetics  
UZ – KU Leuven  
Herestraat 49 PO box 602  
B-3000 Leuven, Belgium  

Mail: [genomicscore@uzleuven.be](mailto:genomicscore@uzleuven.be "")


## Citing GBSX

We ask that you cite this paper if you use GBSX in work that leads to 
publication.

>    Herten,K et al. (2015) GBSX: a toolkit for experimental design and 
demultiplexing genotyping by sequencing experiments  BMC Bioinformatics 2015, 16:73 doi:10.1186/s12859-015-0514-3

## Help

Genotyping By Sequencing demultipleXing toolkit (GBSX) is a toolkit with an 
inline barcode demultiplexer for usage in the analysis of single read or 
paired-end genotyping by sequence (GBS) data, a barcode generator, a barcode 
discovery tool, and a restriction enzyme predictor. GBSX can easily be 
incorperated as a preceding analysis step for already deployed SNP pipelines.

### Restriction Enzyme Predictor
 
mandatory parameters:  
*    `-d`    digest sequence  
*     `-l`     read length   
*     `-f`     file of reference fasta file location(s)  

optional parameters:  
*    `-e`    enzyme name to use (default: Enzyme)  
*     `-g`     genome name to use in bed file name (default: genome)  
*     `-n`     minimum size fragments to include (default: 100)  
*     `-m`     maximum size fragments to use (default: 1000)  
*     `-E`     second enzyme name to use (default: Enzyme2)  
*     `-D`     digest sequence for a second enzyme (default: not declared)  
*     `-R`     digest sequence for a third enzyme (default: not declared)  
 

### Barcode Generator

mandatory parameters:  
*    `-b`    the number of barcodes needed   
*    `-e`    the enzyme used for the experiment   
  
optional parameters:  
*    `-ef`    the enzyme file. This option adds new enzymes.   
                The file must be tab delimited: First column the enzyme name, 
second column the cutsites remains (comma separated) (example: enzyme ApeKI and 
restriction site G^CWGC: "ApeKI \tab CAGC,CTGC").   
*    `-nb`    the maximum number of bootstraps that must be executed. 
(optional, standard 10000). By the start of a new bootstrap a complete new 
design is made. The best scored design (most random barcodes and best scored 
bases distribution is kept as result)   
*    `-bt`    the number of barcode tries. (standard 20) If a random barcode 
does not fit into the current design try this number of times with a new random 
barcode before restarting the bootstrap.  
*    `-o`    the output directory (standard current working directory)   
*    `-us`    try to find the ultime match: the best barcode combination with 
the best bases distribution (standard false)  true: continue even when the right
 number of barcodes is found.  
*    `-bf`    a file with all barcodes that are used as basic set (this file is 
one of the possible output files)  
*    `-nf`    a file with all barcodes that may not be used in the design. If 
this file contains barcodes that are also found in the basic set file, these 
barcodes will be replaced in the design by new random barcodes.  


### Demultiplexer

This program demultiplexes fastq or fastq.gz files obtained from sequencing with 
inline barcodes.  
Like used in GBS, RAD, ... protocols.

These parameters are mandatory:   
*    `-f1`    the name and path of the fastq or fastq.gz file to demultiplex  
*    `-i`    the name and path of the info file. This is a tab delimeted file 
without headings, with three columns: sample, sequence of the barcode, name of 
the enzyme  
  
These parameters are optional:   
*    `-f2`    the name of the second fastq or fastq.gz file (only with 
paired-end sequencing)  
*    `-o`    the name of the output directory (standard the directory of the 
call)  
*    `-lf`    use long file names (standard false) filename is standard the 
sample name, long file names is sample name _ barcode _ enzyme	   
*    `-rad`    if the data is rad data or not (-rad true for RAD data, -rad 
false for GBS data) standard false (GBS)  
*    `-gzip`    the input and output are/must be gziped (.gz) (standard false: 
input and output are .fastq, if true this is .fastq.gz)  
*    `-m`    the allowed mismatches in the barcodes + enzymes (standard this 
value is 1)  
*    `-mb`    the allowed mismatches in the barcodes (overrides the option -m)  
*    `-me`    the allowed mismatches in the enzymes (overrides the option -m)  
*    `-minsl`    the minimum allowed length for the sequences (standard 0, 
rejected sequences are found in the stats for each sample in the rejected.count 
column. The sequences are found untrimmed in the undetermined file.)  
*    `-n`    keep sequences where N occurs as a "nucleotide" (standard true)  
*    `-ca`    the common adaptor used in the sequencing (standard (only first 
piece) AGATCGGAAGAGCG) currently only used for adaptor ligase see -al and when 
-rad is true) (minimum length is 10)  
*    `-s`    the posible distance of the start. This is the distance count from 
the start of the read to the first basepair of the barcode or enzyme 
(standard 0, maximum 20)  
*    `-cc`    Checks the complete read for the enzyme (if false, stops at the 
first possible enzyme cutsite) (use values true or false, standard is true). If 
used, the sequence after the enzyme site is compared to the adaptors, if the 
first basepairs of the sequence are compaired to the first basepairs of the 
adaptor  
*    `-kc`    Keep the enzyme cut-site remains (standard true) (example: enzyme 
ApeKI and restriction site G^CWGC: "ApeKI \tab CAGC,CTGC")  
*    `-ea`    Add enzymes from the given file (keeps the standard enzymes, and 
add the new) (enzyme file: no header, enzyme name tab cutsites (multiple 
cutsites are comma separeted)) (only use once, not use -er) (example: enzyme 
ApeKI and restriction site G^CWGC: "ApeKI \tab CAGC,CTGC")  
*    `-er`    Replace enzymes from the given file (do not keep the standard 
enzymes) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites 
are comma separeted)) (only use once, not use -ea)  
*    `-al`    check for adaptor ligase: no (for no check) or a positive integer 
(starts at 0), for the number of mismatches (only checks 10 basepairs of 
the adaptor), standard 1  
*    `-scb`    Use self correcting barcodes (barcodes created by the 
barcodeGenerator) (standard false)  
*    `-malg`    the used algorithm to find mismatches and indels, possible 
algorithms:   
    *    `hammings (Standard)`    Checks for mismatches (no indels)  
    *    `knuth`    Faster than hammings, but can miss some locations  
    *    `indelmis`    Checks for mismatches and indels, the barcode/enzyme/
adaptor with the least errors (mismatches or indels) is taken  
    *    `misindel`    Checks for mismatches and indels, the mismatches are 
supperior to the indels (faster than indelmis, but errors can be higher)  
  
*    `-q`    the kind of quality scores used in the fastq file (including how 
phred scores are encoded):   
    *    `Illumina1.8 (Standard)`  
    *    `Illumina1.5`  
    *    `Illumina1.3`  
    *    `Sanger`  
    *    `Solid`  
  
  
  
Possible Standard Enzymes for the info file: (NA is no enzyme)  
* `ApeKI`  
* `PstI`  
* `EcoT22I`  
* `PasI`  
* `HpaII`  
* `MspI`  
* `PstI-EcoT22I`  
* `PstI-MspI`  
* `PstI-TaqI`  
* `SbfI-MspI`  
* `AsiSI-MspI`  
* `BssHII-MspI`  
* `FseI-MspI`  
* `SalI-MspI`  
* `ApoI`  
* `BamHI`  
* `MseI`  
* `Sau3AI`  
* `RBSTA`  
* `RBSCG`  
* `NspI`  
* `NAN`  
  
  
  
### Barcode Discovery

This program searches for possible barcodes and barcode enzyme combinations.   
Designed for the discovery of sequencing errors, or unused barcodes when a large
proportion of the demultiplex is undetermined.  

Mandatory parameters:  
* `-f1`    the name of the input file (mandatory)   
Optional parameters:  
* `-min`    the minimum length of the barcode (standard 6)   
* `-max`    the maximum length of the barcode (standard 10)   
* `-gzip`    use gzip files as input and output (standard false)   
* `-o`    the output directory (standard the directory of execution)   
* `-ea`    Add enzymes from the given file (keeps the standard enzymes, and add 
the new) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites 
are comma separeted)) (example: enzyme ApeKI and restriction site G^CWGC: 
"ApeKI \tab CAGC,CTGC") (only use once, not use -er)  
* `-er`    Replace enzymes from the given file (do not keep the standard 
enzymes) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites 
are comma separeted)) (example: enzyme ApeKI and restriction site G^CWGC: 
"ApeKI \tab CAGC,CTGC") (only use once, not use -ea)  
* `-barmin`    The minimum occurance of a barcode before it is shown in the 
results (standard: 200)   
* `-barmax`    The maximum occurance of barcodes shown in the output (increasing
 this number will increase ram usage, but gives a slightly better result) 
(standard: 100)  
* `-barmis`    The percentage of mismatches that may occure between barcodes 
(integer between 1 and 10) (standard: 10)  

## Tutorial

See the Tutorial file and the example folder.



## Change Logs

v1.0
* The original version

v1.0.1
* While demultiplexing, the number of demultiplexed reads are shown for 
  every 100000 reads