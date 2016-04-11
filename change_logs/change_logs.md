---
layout: default
title: Change Logs
permalink: /change_logs/
---

## Change Logs

v1.0
* The original version  

v1.0.1
* While demultiplexing, the number of demultiplexed reads are shown for 
  every 100000 reads

v1.1
* Possible to simulate and demultiplex dual barcode experiments (in paired end modus only)
* Updated barcode recognition for paired end modus in the demultiplexer: when a read can be assigned to multiple samples, 
the read is considered as unvalid (previous was first sample)

v1.1.1
* Updated output and stats for dual barcode experiments

v1.1.2
* Updated barcode recognition for single read modus in the demultiplexer: when a read can be assigned to multiple samples, 
the read is considered as unvalid (previous was first sample)

v1.1.3
* On request added the enzyme AvaII

v1.1.4
* Update adaptor ligase finding algorithm
* Removed unneeded, confusing parameters -cc and -al
* Removed unused code

v1.1.5
 * Update digest (removed possible input file parsing error)
 * Updated Single Read Demultiplexing

v1.2
 * Deleted Demultiplexer option -m
 * Code Clean-up
