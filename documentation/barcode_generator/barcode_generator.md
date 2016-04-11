---
layout: default
title: Barcode Generator
permalink: /documentation/barcode_generator
---

## Barcode Generator

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

