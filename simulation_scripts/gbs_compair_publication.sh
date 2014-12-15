#!/bin/bash
workingdir="/home/koen/tmp/gbs_test";
gbsdir="/home/koen/software/GBSX-1.0/releases/latest";
stacksdir="/home/koen/software/stacks-1.21";
sabredir="/home/koen/software/sabre-master";
javadir="/usr/bin";
cd $workingdir;

paired="true";
readlength="100";
#depth="3";
simulation_number=10;

counter=0;
while [ $counter -lt $simulation_number ]
do
    cd $workingdir;
    counter=`expr $counter + 1`;
    mkdir test_$counter;
    cd test_$counter;
    barcodeNr=`shuf -i 5-20 -n 1`;
    depth=`shuf -i 2-10 -n 1`;
    echo "[COMPAIR][RUN]number=$counter;barcodes=$barcodeNr;depth=$depth;";
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --BarcodeGenerator -b $barcodeNr -e ApeKI;

echo "[INFO]Barcodes=$barcodeNr";
echo "[INFO]depth=$depth";

    awk '{OUT="barcode.stacks.apeki." length($0) ".txt";}; OUT{print >OUT}' barcode_list.txt;
    awk '{print $1 "\t" $1 "\tApeKI"}' barcode_list.txt > barcode.gbsx.apeki.txt;
	mv barcode_list.txt barcode_list.apeki.txt;
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --BarcodeGenerator -b $barcodeNr -e NA;
    awk '{OUT="barcode.stacks.na." length($0) ".txt";}; OUT{print >OUT}' barcode_list.txt;
    awk '{print $1 "\t" $1 "\tNA"}' barcode_list.txt > barcode.gbsx.NA.txt;
    awk '{print $1 "\t" $1 ".R1.fastq\t" $1 ".R2.fastq"}' barcode_list.txt > barcode.sabre.pe.txt;
    awk '{print $1 "\t" $1 ".R1.fastq"}' barcode_list.txt > barcode.sabre.se.txt;
	mv barcode_list.txt barcode_list.na.txt;

    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --GBSsimulator -f $workingdir/chr21.out.fa -p $paired -l $readlength -rpl $depth -b barcode_list.apeki.txt;
	mv output.R1.fastq output.apeki.R1.fastq;
	mv output.R2.fastq output.apeki.R2.fastq;
	mv output.error.txt output.apeki.error.txt;
	
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --GBSsimulator -f $workingdir/chr21.out.fa -p $paired -l $readlength -rpl $depth -b barcode_list.na.txt;
	mv output.R1.fastq output.na.R1.fastq;
	mv output.R2.fastq output.na.R2.fastq;
	mv output.error.txt output.na.error.txt;
    totalCounts=`cat output.apeki.R1.fastq | wc -l`;
    totalReads=`expr $totalCounts / 4`;
    readCounts=`expr $totalReads / $barcodeNr`;
    longestBarcodeLength=`awk 'BEGIN{longest=0;}{if (length($1) > longest){longest=length($1);};}END{printf longest;}' barcode_list.apeki.txt`;
    trimlengths=`expr $readlength - $longestBarcodeLength`;

###############################
####
####
####  PAIRED END
####
####
###############################


###############################
####
####
####  PAIRED END
####
####
####  GBS
####
###############################

##### GBSX GBS paired end simulation 
	mkdir GBSX_GBS_pe;

    echo -n "[COMPAIR] GBSX_GBS start: "
    date;
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --Demultiplexer -f1 output.apeki.R1.fastq -f2 output.apeki.R2.fastq -i barcode.gbsx.apeki.txt -o GBSX_GBS_pe;
    echo -n "[COMPAIR] GBSX_GBS end: "
    date;

    for i in `ls GBSX_GBS_pe/*.R1.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/GBSX\_GBS\_pe\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> GBSX_GBS_pe.sum.tsv;
            echo $readCounts >> GBSX_GBS_pe.sum.tsv;
        fi
    done

	rm -r GBSX_GBS_pe;


###############################
####
####
####  PAIRED END
####
####
####  RAD
####
###############################

##### GBSX RAD paired end simulation
	mkdir GBSX_RAD_pe;

    echo -n "[COMPAIR] GBSX_RAD start: "
    date;
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --Demultiplexer -f1 output.apeki.R1.fastq -f2 output.apeki.R2.fastq -i barcode.gbsx.apeki.txt -o GBSX_RAD_pe -rad true;
    echo -n "[COMPAIR] GBSX_RAD end: "
    date;

    for i in `ls GBSX_RAD_pe/*.R1.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/GBSX\_RAD\_pe\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> GBSX_RAD_pe.sum.tsv;
            echo $readCounts >> GBSX_RAD_pe.sum.tsv;
        fi
    done

	rm -r GBSX_RAD_pe;

#### Stacks RAD paired end simulation
	mkdir stacks_radtags_pe;

    echo -n "[COMPAIR] stacks_radtags start: "
    date;
    for i in `ls barcode.stacks.apeki.*.txt`;
    do
        $stacksdir/process_radtags -1 output.apeki.R1.fastq -2 output.apeki.R2.fastq -o stacks_radtags_pe/ -e apeKI -b $i --adapter_1 AGATCGGAAGAGCG --adapter_2 AGATCGGAAGAGCG --adapter_mm 1 --inline_null -s 0;
    done
    echo -n "[COMPAIR] stacks_radtags end: "
    date;

    for i in `ls stacks_radtags_pe/*.1.fq`;
    do
        if [[ $i != *rem* ]]
        then
            barcode=`echo $i | sed 's/stacks\_radtags\_pe\///g' | sed 's/\.1\.fq//g' | sed 's/sample\_//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> stacks_radtags_pe.sum.tsv;
            echo $readCounts >> stacks_radtags_pe.sum.tsv;
        fi
    done

	rm -r stacks_radtags_pe;


###############################
####
####
####  PAIRED END
####
####
####  Inline Barcodes
####
###############################

##### GBSX inline barcode paired end simulation
	mkdir GBSX_NA_pe;

    echo -n "[COMPAIR] GBSX_NA start: "
    date;
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --Demultiplexer -f1 output.na.R1.fastq -f2 output.na.R2.fastq -i barcode.gbsx.NA.txt -o GBSX_NA_pe -rad true;
    echo -n "[COMPAIR] GBSX_NA end: "
    date;

    for i in `ls GBSX_NA_pe/*.R1.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/GBSX\_NA\_pe\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> GBSX_NA_pe.sum.tsv;
            echo $readCounts >> GBSX_NA_pe.sum.tsv;
        fi
    done

	rm -r GBSX_NA_pe;


##### Stacks inline barcode paired end simulation
	mkdir stacks_shortreads_pe;

    echo -n "[COMPAIR] stacks_shortreads start: "
    date;
    for i in `ls barcode.stacks.na.*.txt`;
    do
        $stacksdir/process_shortreads -1 output.na.R1.fastq -2 output.na.R2.fastq -o stacks_shortreads_pe/ -b $i --adapter_1 AGATCGGAAGAGCG --adapter_2 AGATCGGAAGAGCG --adapter_mm 1 --inline_null -s 0;
    done
    echo -n "[COMPAIR] stacks_shortreads end: "
    date;

    for i in `ls stacks_shortreads_pe/*.1.fq`;
    do
        if [[ $i != *rem* ]]
        then
            barcode=`echo $i | sed 's/stacks\_shortreads\_pe\///g' | sed 's/\.1\.fq//g' | sed 's/sample\_//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> stacks_shortreads_pe.sum.tsv;
            echo $readCounts >> stacks_shortreads_pe.sum.tsv;
        fi
    done
	
	rm -r stacks_shortreads_pe;

##### Sabre inline barcode paired end simulation
	mkdir sabre_pe;

    echo -n "[COMPAIR] sabre start: "
    date;
	cd sabre_pe;
    $sabredir/sabre pe -f ../output.na.R1.fastq -r ../output.na.R2.fastq -b ../barcode.sabre.pe.txt -m 1 -u undetermined.R1.fastq -w undetermined.R2.fastq;
	cd ..;
    echo -n "[COMPAIR] sabre end: "
    date;

    for i in `ls sabre_pe/*.R1.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/sabre\_pe\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> sabre_pe.sum.tsv;
            echo $readCounts >> sabre_pe.sum.tsv;
        fi
    done

	rm -r sabre_pe;


###############################
####
####
####  SIGNLE READS
####
####
####  GBS
####
###############################
    
##### GBSX GBS signle read simulation
	mkdir GBSX_GBS_sr;

    echo -n "[COMPAIR] GBSX_GBS start: "
    date;
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --Demultiplexer -f1 output.apeki.R1.fastq -i barcode.gbsx.apeki.txt -o GBSX_GBS_sr;
    echo -n "[COMPAIR] GBSX_GBS end: "
    date;

    for i in `ls GBSX_GBS_sr/*.R1.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/GBSX\_GBS\_sr\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> GBSX_GBS_sr.sum.tsv;
            echo $readCounts >> GBSX_GBS_sr.sum.tsv;
        fi
    done

	rm -r GBSX_GBS_sr;

###############################
####
####
####  SIGNLE READS
####
####
####  RAD
####
###############################


##### GBSX RAD signle read simulation
	mkdir GBSX_RAD_sr;

    echo -n "[COMPAIR] GBSX_RAD start: "
    date;
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --Demultiplexer -f1 output.apeki.R1.fastq -i barcode.gbsx.apeki.txt -o GBSX_RAD_sr -rad true;
    echo -n "[COMPAIR] GBSX_RAD end: "
    date;

    for i in `ls GBSX_RAD_sr/*.R1.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/GBSX\_RAD\_sr\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> GBSX_RAD_sr.sum.tsv;
            echo $readCounts >> GBSX_RAD_sr.sum.tsv;
        fi
    done

	rm -r GBSX_RAD_sr;

##### Stacks RAD signle read simulation
	mkdir stacks_radtags_sr;

    echo -n "[COMPAIR] stacks_radtags start: "
    date;
    for i in `ls barcode.stacks.apeki.*.txt`;
    do
        $stacksdir/process_radtags -f output.apeki.R1.fastq -o stacks_radtags_sr/ -e apeKI -b $i --adapter_1 AGATCGGAAGAGCG --adapter_2 AGATCGGAAGAGCG --adapter_mm 1 --inline_null -s 0;
    done
    echo -n "[COMPAIR] stacks_radtags end: "
    date;

    for i in `ls stacks_radtags_sr/*.fq`;
    do
        if [[ $i != *rem* ]]
        then
            barcode=`echo $i | sed 's/stacks\_radtags\_sr\///g' | sed 's/\.fq//g' | sed 's/sample\_//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> stacks_radtags_sr.sum.tsv;
            echo $readCounts >> stacks_radtags_sr.sum.tsv;
        fi
    done

	rm -r stacks_radtags_sr;



###############################
####
####
####  SIGNLE READS
####
####
####  Inline Barcodes
####
###############################


##### GBSX inline barcode signle read simulation
	mkdir GBSX_NA_sr;

    echo -n "[COMPAIR] GBSX_NA start: "
    date;
    $javadir/java -jar $gbsdir/GBSX_v1.0.jar --Demultiplexer -f1 output.na.R1.fastq -i barcode.gbsx.NA.txt -o GBSX_NA_sr -rad true;
    echo -n "[COMPAIR] GBSX_NA end: "
    date;

    for i in `ls GBSX_NA_sr/*.R1.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/GBSX\_NA\_sr\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> GBSX_NA_sr.sum.tsv;
            echo $readCounts >> GBSX_NA_sr.sum.tsv;
        fi
    done

	rm -r GBSX_NA_sr;


##### Stacks inline barcode signle read simulation
	mkdir stacks_shortreads_sr;

    echo -n "[COMPAIR] stacks_shortreads start: "
    date;
    for i in `ls barcode.stacks.na.*.txt`;
    do
        $stacksdir/process_shortreads -f output.na.R1.fastq -o stacks_shortreads_sr/ -b $i --adapter_1 AGATCGGAAGAGCG --adapter_2 AGATCGGAAGAGCG --adapter_mm 1 --inline_null -s 0;
    done
    echo -n "[COMPAIR] stacks_shortreads end: "
    date;

    for i in `ls stacks_shortreads_sr/*.fq`;
    do
        if [[ $i != *rem* ]]
        then
            barcode=`echo $i | sed 's/stacks\_shortreads\_sr\///g' | sed 's/\.fq//g' | sed 's/sample\_//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> stacks_shortreads_sr.sum.tsv;
            echo $readCounts >> stacks_shortreads_sr.sum.tsv;
        fi
    done

	rm -r stacks_shortreads_sr;

##### sabre inline barcode signle read simulation
	mkdir sabre_sr;

    echo -n "[COMPAIR] sabre start: "
    date;
	cd sabre_sr;
    $sabredir/sabre se -f ../output.na.R1.fastq -b ../barcode.sabre.se.txt -m 1 -u undetermined.R1.fastq;
	cd ..;
    echo -n "[COMPAIR] GBSX_NA end: "
    date;

    for i in `ls sabre_sr/*.fastq`;
    do
        if [[ $i != *undetermined* ]]
        then
            barcode=`echo $i | sed 's/sabre\_sr\///g' | sed 's/\.R1\.fastq//g'`;
            perl $workingdir/compair_demultiplex.pl $workingdir/chr21.out.fa $i $trimlengths $barcode >> sabre_sr.sum.tsv;
            echo $readCounts >> sabre_sr.sum.tsv;
        fi
    done

	rm -r sabre_sr;



###### Remove all simulated files


done


cd $workingdir;

counter=0;
while [ $counter -lt $simulation_number ]
do
	counter=`expr $counter + 1`;
	echo $counter;
	cat $workingdir/test_$counter/GBSX_GBS_pe.sum.tsv >> GBSX_GBS_pe.sum.tsv;
	cat $workingdir/test_$counter/GBSX_RAD_pe.sum.tsv >> GBSX_RAD_pe.sum.tsv;
	cat $workingdir/test_$counter/GBSX_NA_pe.sum.tsv >> GBSX_NA_pe.sum.tsv;
	cat $workingdir/test_$counter/stacks_radtags_pe.sum.tsv >> stacks_radtags_pe.sum.tsv;
	cat $workingdir/test_$counter/stacks_shortreads_pe.sum.tsv >> stacks_shortreads_pe.sum.tsv;
	cat $workingdir/test_$counter/sabre_pe.sum.tsv >> sabre_pe.sum.tsv;

	cat $workingdir/test_$counter/GBSX_GBS_sr.sum.tsv >> GBSX_GBS_sr.sum.tsv;
	cat $workingdir/test_$counter/GBSX_RAD_sr.sum.tsv >> GBSX_RAD_sr.sum.tsv;
	cat $workingdir/test_$counter/GBSX_NA_sr.sum.tsv >> GBSX_NA_sr.sum.tsv;
	cat $workingdir/test_$counter/stacks_radtags_sr.sum.tsv >> stacks_radtags_sr.sum.tsv;
	cat $workingdir/test_$counter/stacks_shortreads_sr.sum.tsv >> stacks_shortreads_sr.sum.tsv;
	cat $workingdir/test_$counter/sabre_sr.sum.tsv >> sabre_sr.sum.tsv;
done

awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' GBSX_GBS_pe.sum.tsv > GBSX_GBS_pe.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' GBSX_GBS_pe.sum.tsv > GBSX_GBS_pe.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' GBSX_GBS_pe.sum.tsv > GBSX_GBS_pe.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' GBSX_RAD_pe.sum.tsv > GBSX_RAD_pe.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' GBSX_RAD_pe.sum.tsv > GBSX_RAD_pe.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' GBSX_RAD_pe.sum.tsv > GBSX_RAD_pe.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' GBSX_NA_pe.sum.tsv > GBSX_NA_pe.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' GBSX_NA_pe.sum.tsv > GBSX_NA_pe.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' GBSX_NA_pe.sum.tsv > GBSX_NA_pe.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' stacks_radtags_pe.sum.tsv > stacks_radtags_pe.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' stacks_radtags_pe.sum.tsv > stacks_radtags_pe.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' stacks_radtags_pe.sum.tsv > stacks_radtags_pe.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' stacks_shortreads_pe.sum.tsv > stacks_shortreads_pe.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' stacks_shortreads_pe.sum.tsv > stacks_shortreads_pe.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' stacks_shortreads_pe.sum.tsv > stacks_shortreads_pe.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' sabre_pe.sum.tsv > sabre_pe.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' sabre_pe.sum.tsv > sabre_pe.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' sabre_pe.sum.tsv > sabre_pe.fault_dem.tsv;

join -j 1 GBSX_GBS_pe.cor_all.tsv GBSX_RAD_pe.cor_all.tsv > GBSX_GBS_RAD_pe.cor_all.tsv;
join -j 1 GBSX_GBS_RAD_pe.cor_all.tsv GBSX_NA_pe.cor_all.tsv > GBSX_pe.cor_all.tsv;
join -j 1 stacks_radtags_pe.cor_all.tsv stacks_shortreads_pe.cor_all.tsv > stacks_pe.cor_all.tsv;
join -j 1 GBSX_pe.cor_all.tsv stacks_pe.cor_all.tsv > GBSX_stacks_pe.cor_all.tsv;
join -j 1 GBSX_stacks_pe.cor_all.tsv sabre_pe.cor_all.tsv > all_pe.cor_all.tsv;
rm GBSX_GBS_RAD_pe.cor_all.tsv; rm GBSX_pe.cor_all.tsv; rm stacks_pe.cor_all.tsv GBSX_stacks_pe.cor_all.tsv;

join -j 1 GBSX_GBS_pe.dem_all.tsv GBSX_RAD_pe.dem_all.tsv > GBSX_GBS_RAD_pe.dem_all.tsv;
join -j 1 GBSX_GBS_RAD_pe.dem_all.tsv GBSX_NA_pe.dem_all.tsv > GBSX_pe.dem_all.tsv;
join -j 1 stacks_radtags_pe.dem_all.tsv stacks_shortreads_pe.dem_all.tsv > stacks_pe.dem_all.tsv;
join -j 1 GBSX_pe.dem_all.tsv stacks_pe.dem_all.tsv > GBSX_stacks_pe.dem_all.tsv;
join -j 1 GBSX_stacks_pe.dem_all.tsv sabre_pe.dem_all.tsv > all_pe.dem_all.tsv;
rm  GBSX_GBS_RAD_pe.dem_all.tsv; rm GBSX_pe.dem_all.tsv; rm stacks_pe.dem_all.tsv GBSX_stacks_pe.dem_all.tsv;

join -j 1 GBSX_GBS_pe.fault_dem.tsv GBSX_RAD_pe.fault_dem.tsv > GBSX_GBS_RAD_pe.fault_dem.tsv;
join -j 1 GBSX_GBS_RAD_pe.fault_dem.tsv GBSX_NA_pe.fault_dem.tsv > GBSX_pe.fault_dem.tsv;
join -j 1 stacks_radtags_pe.fault_dem.tsv stacks_shortreads_pe.fault_dem.tsv > stacks_pe.fault_dem.tsv;
join -j 1 GBSX_pe.fault_dem.tsv stacks_pe.fault_dem.tsv > GBSX_stacks_pe.fault_dem.tsv;
join -j 1 GBSX_stacks_pe.fault_dem.tsv sabre_pe.fault_dem.tsv > all_pe.fault_dem.tsv;
rm GBSX_GBS_RAD_pe.fault_dem.tsv; rm stacks_pe.fault_dem.tsv GBSX_stacks_pe.fault_dem.tsv;



awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' GBSX_GBS_sr.sum.tsv > GBSX_GBS_sr.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' GBSX_GBS_sr.sum.tsv > GBSX_GBS_sr.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' GBSX_GBS_sr.sum.tsv > GBSX_GBS_sr.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' GBSX_RAD_sr.sum.tsv > GBSX_RAD_sr.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' GBSX_RAD_sr.sum.tsv > GBSX_RAD_sr.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' GBSX_RAD_sr.sum.tsv > GBSX_RAD_sr.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' GBSX_NA_sr.sum.tsv > GBSX_NA_sr.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' GBSX_NA_sr.sum.tsv > GBSX_NA_sr.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' GBSX_NA_sr.sum.tsv > GBSX_NA_sr.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' stacks_radtags_sr.sum.tsv > stacks_radtags_sr.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' stacks_radtags_sr.sum.tsv > stacks_radtags_sr.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' stacks_radtags_sr.sum.tsv > stacks_radtags_sr.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' stacks_shortreads_sr.sum.tsv > stacks_shortreads_sr.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' stacks_shortreads_sr.sum.tsv > stacks_shortreads_sr.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' stacks_shortreads_sr.sum.tsv > stacks_shortreads_sr.fault_dem.tsv;
awk 'BEGIN{counter=1;} {if ($2 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($2/$5);} counter++;}' sabre_sr.sum.tsv > sabre_sr.cor_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $5 == 0){print counter"\t"0;}else{print counter"\t"($3/$5);} counter++;}' sabre_sr.sum.tsv > sabre_sr.dem_all.tsv;
awk 'BEGIN{counter=1;} {if ($3 == 0 || $4 == 0){print counter"\t"0;}else{print counter"\t"($4/$3);} counter++;}' sabre_sr.sum.tsv > sabre_sr.fault_dem.tsv;

join -j 1 GBSX_GBS_sr.cor_all.tsv GBSX_RAD_sr.cor_all.tsv > GBSX_GBS_RAD_sr.cor_all.tsv;
join -j 1 GBSX_GBS_RAD_sr.cor_all.tsv GBSX_NA_sr.cor_all.tsv > GBSX_sr.cor_all.tsv;
join -j 1 stacks_radtags_sr.cor_all.tsv stacks_shortreads_sr.cor_all.tsv > stacks_sr.cor_all.tsv;
join -j 1 GBSX_sr.cor_all.tsv stacks_sr.cor_all.tsv > GBSX_stacks_sr.cor_all.tsv;
join -j 1 GBSX_stacks_sr.cor_all.tsv sabre_sr.cor_all.tsv > all_sr.cor_all.tsv;
rm GBSX_GBS_RAD_sr.cor_all.tsv; rm GBSX_sr.cor_all.tsv; rm stacks_sr.cor_all.tsv GBSX_stacks_sr.cor_all.tsv;

join -j 1 GBSX_GBS_sr.dem_all.tsv GBSX_RAD_sr.dem_all.tsv > GBSX_GBS_RAD_sr.dem_all.tsv;
join -j 1 GBSX_GBS_RAD_sr.dem_all.tsv GBSX_NA_sr.dem_all.tsv > GBSX_sr.dem_all.tsv;
join -j 1 stacks_radtags_sr.dem_all.tsv stacks_shortreads_sr.dem_all.tsv > stacks_sr.dem_all.tsv;
join -j 1 GBSX_sr.dem_all.tsv stacks_sr.dem_all.tsv > GBSX_stacks_sr.dem_all.tsv;
join -j 1 GBSX_stacks_sr.dem_all.tsv sabre_sr.dem_all.tsv > all_sr.dem_all.tsv;
rm  GBSX_GBS_RAD_sr.dem_all.tsv; rm GBSX_sr.dem_all.tsv; rm stacks_sr.dem_all.tsv GBSX_stacks_sr.dem_all.tsv;

join -j 1 GBSX_GBS_sr.fault_dem.tsv GBSX_RAD_sr.fault_dem.tsv > GBSX_GBS_RAD_sr.fault_dem.tsv;
join -j 1 GBSX_GBS_RAD_sr.fault_dem.tsv GBSX_NA_sr.fault_dem.tsv > GBSX_sr.fault_dem.tsv;
join -j 1 stacks_radtags_sr.fault_dem.tsv stacks_shortreads_sr.fault_dem.tsv > stacks_sr.fault_dem.tsv;
join -j 1 GBSX_sr.fault_dem.tsv stacks_sr.fault_dem.tsv > GBSX_stacks_sr.fault_dem.tsv;
join -j 1 GBSX_stacks_sr.fault_dem.tsv sabre_sr.fault_dem.tsv > all_sr.fault_dem.tsv;
rm GBSX_GBS_RAD_sr.fault_dem.tsv; rm stacks_sr.fault_dem.tsv GBSX_stacks_sr.fault_dem.tsv;

Rscript create_plots.R;
