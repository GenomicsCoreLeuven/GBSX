---
layout: default
title: Demultiplexer
permalink: /documentation/demultiplexer
---

## Demultiplexer

This program demultiplexes fastq or fastq.gz files obtained from sequencing with 
inline barcodes.  
Like used in GBS, RAD, ... protocols.

These parameters are mandatory:   
*    `-f1`    the name and path of the fastq or fastq.gz file to demultiplex  
*    `-i`    the name and path of the info file. This is a tab delimeted file 
without headings, with three (or more) columns: sample, sequence of the barcode, 
name of the enzyme, name of the second enzyme (optional, can be an empty 
string), the second barcode (optional, can be an empty string),mismatches 
for the barcode (optional)
 
  
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
*    `-kc`    Keep the enzyme cut-site remains (standard true) (example: enzyme 
ApeKI and restriction site G^CWGC: "ApeKI \tab CAGC,CTGC")  
*    `-ea`    Add enzymes from the given file (keeps the standard enzymes, and 
add the new) (enzyme file: no header, enzyme name tab cutsites (multiple 
cutsites are comma separeted)) (only use once, not use -er) (example: enzyme 
ApeKI and restriction site G^CWGC: "ApeKI \tab CAGC,CTGC")  
*    `-er`    Replace enzymes from the given file (do not keep the standard 
enzymes) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites 
are comma separeted)) (only use once, not use -ea) 
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
* `AvaII`  
* `NAN`  
  
  
