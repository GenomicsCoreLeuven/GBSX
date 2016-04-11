---
layout: default
title: Barcode Discovery
permalink: /documentation/barcode_discovery
---  
  
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

