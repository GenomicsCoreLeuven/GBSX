---
layout: default
title: Tutorial
permalink: /tutorial/
---


Tutorial
========

For the tutorial we use 2 small human chromosomes (hg19): chr21 and chr22.
These chromosomes where downloaded from UCSC:\\
* [Human chr21 of USCS](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz "Human chr21")  
* [Human chr22 of USCS](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz "Human chr22")  

To run the tutorial yourself please download these two fasta files and edit a location file (example format as _chromosomeLinks.txt_) indicating where these files are located on your machine.

The restriction enzyme used is ApekI. The cutsite is G^CWGC. More information about this enzyme can be found on [the neb website](https://www.neb.com/products/r0643-apeki "the neb website about ApekI").

### 1) GBSX Digest
For the usage of this perl module, you must have BioPerl installed (installation instructions on [http://www.bioperl.org/wiki/Installing_BioPerl](http://www.bioperl.org/wiki/Installing_BioPerl, "http://www.bioperl.org/wiki/Installing_BioPerl")).
For running the script we have to configure a txt-file with on each line the full path to a fasta file to use in the digest (see the file chromosomeLinks.txt as an example, DO edit this file to indicate the directories on your machine).
This script can easily been run by:  
```
perl GBSX_digest_v1.0.pl -d G^CWGC -l 100 -f chromosomeLinks.txt
```  
where -d is the digest site, -l is the sequencing length, -f is the file edited for your fasta file links.

This will generate 2 files (the output can be found under the 1_digest folder):  
* _genome.Enzyme.100nt.digest.bed_  
* _genome.Enzyme.100nt.digest_results_  
The bed file contains all start-stop locations of the sequenced portion (based on a read length of 100bp) of 
all fragments that are above or equal to 100bp and below or equal to 1000bp in length.
The results file give a brief summary of the fragments, including the total number of fragments and the 
number of fragments between 100-1000bp.  For fragments between 100-1000bp, the results file also indicates 
the number of fragments per chromosome, the distribution of distances between fragments, and the distribution of fragment lengths.

When using this digest-tool for research, it is possible to determine the SNPs that will be sequenced. 
This generated bed file can be used with bedtools, dbSNP, and can be uploaded in ucsc as a custom track. 
This can be illustrated by downloading all known [snp positions from dbSNP or ucsc](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/snp138Common.txt.gz "UCSC hg19 snp138 common SNPs"). (filtered version can be found in _A_chromosomes_)
From this file a bed file is created for the chromosomes 21 and 22.  
```
    zcat snp138Common.txt.gz | awk '{OFS="\t"; print $2,$3,$4,$5}' | awk '{if ($1 == 22 || $1 == 21) {print $0}}' > snp138.chr21-22.bed
```  
Bedtools can be used to get the snp positions:  
```
    /home/GenomicsCore/bedtools/intersectBed -a 1_digest/genome.Enzyme.100nt.digest.bed -b A_chromosomes/snp138.chr21-22.bed > snps.bed
```  
In total 90271 snps are found in these GBS fragments.

### 2) GBSX Generate barcodes
For the creation of unique barcodes, the BarcodeGenerator is used. 
For this tutorial 5 barcodes will be created, using the ApeKI enzyme:  
```
    java -jar ../../releases/latest/GBSX_v1.0.jar --BarcodeGenerator -b 5 -e ApeKI
```  
The result of this tool are 2 files:  
* _barcode_list.txt_        where all barcode sequences are listed  
* _barcode_summary.txt_     where there is a summary of the creation, as well as a base occurance matrix.  

The outcome of this tool will always be different. This is because the generated barcodes
are completely random. The best practice is to manually check the position matrix to see if every bases occurs equally on every position. 

### 3) GBSX Demultiplex
Because of the random barcodes, use for this part of the tutorial the barcodes in the _2_barcodes_ folder.
For the demultiplexing we adjust the original _barcode_list.txt_ file. Sample names are added in the first column, the second column will be 
the generated barcodes, the last column will be the restriction enzyme, here ApeKI.  

In _3_simulated_data_, a possible dataset is found to demultiplex. This dataset is based on fragments between 50bp and 200bp from chr21 and chr22 (derived from the digest script).
The dataset has a depth of 10, with possible snps and errors.  

For the demultiplexing, the standard parameters should be fine (starting in the example folder):  
```
    java -jar ../releases/latest/GBSX_v1.0.jar --Demultiplexer -f1 B_simulated_data/simulation.R1.fastq.gz -f2 B_simulated_data/simulation.R2.fastq.gz -i 3_demultiplex/barcode_list.txt -gzip true -o 3_demultiplex/
```  


### 4) GBSX Barcode predictor
For the illustration of the barcode predictor, the output of the simulator is analyzed. With the standard parameters this results in the same stats as for the demultiplexing (for no mismatches).  
```
    java -jar ../releases/latest/GBSX_v1.0.jar --BarcodeDiscovery -f1 B_simulated_data/simulation.R1.fastq.gz -gzip true -max 16 -o 4_predictor
```  
when the results are filterd on the used enzyme ApeKI, the used barcodes are the most used.
Also when there is searched for barcodes without the use of an enzyme, the used barcodes pop up (with a C at the end, the first base of the restriction site).

The result of this tool are multiple files starting with counts. These files have the most counted sequences with as length defined in the file name.
A more important file is the possibleBarcodes.cnt file. This file contains the sequences that are most likely to be barcodes, the number of occurances, and 
(if asked) the possible enzyme(s).


