/*
 * This is GBSX v1.0. A toolkit for experimental design and demultiplexing genotyping by sequencing experiments. \n *  \n * Copyright 2014 KU Leuven
 * 
 * 
 * This file is part of GBSX.
 *
 * GBSX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GBSX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GBSX.  If not, see <http://www.gnu.org/licenses/>.
 */
package be.uzleuven.gc.logistics.GBSX.demultiplexer.application;

import be.uzleuven.gc.logistics.GBSX.GBSX;
import be.uzleuven.gc.logistics.GBSX.barcodeCorrector.application.BarcodeCorrectingAlgorithm;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.GBSdemultiplex;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.application.distanceAlgorithms.FindingDistanceAlgorithm;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.application.distanceAlgorithms.MismatchIndelDistance;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.ErrorInLogException;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.InvalidReadEnum;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.InvalidReadException;
import be.uzleuven.gc.logistics.GBSX.utils.exceptions.StopExcecutionException;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.infrastructure.EnzymeFileParser;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqBufferedReader;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqBufferedWriter;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.infrastructure.InfoFileParser;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.infrastructure.fileInteractors.LoggerFile;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.BasePair;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqParts;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.DemultiplexArguments;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.DemultiplexParameters;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeComparator;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeEnum;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqPairBufferedReader;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqPairBufferedWriter;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqRead;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author Koen Herten for the KU Leuven
 */
public final class FastqDemultiplex {
    
    private final int MAXIMUM_DISTANCE_BETWEEN_START_AND_BARCODE = 20;
    private final int PROCESS_REPORTING_INTERVAL = 100000;
    
    private ArrayList<Sample> sampleList;
    private int longestBarcodeLength;
    private final LoggerFile loggerFile;
    /**
     * correction log is only used in the correction of the pair end demultiplex
     */
    private final CorrectionLog correctionLog;
    
    private final DemultiplexParameters parameters;
    private final FindingDistanceAlgorithm findingDistanceAlgorithm;
    
    /**
     * this makes the bases for the demultiplexing. The given array contains all the needed arguments.
     * <br> If not, there will be thrown exceptions like IllegalArgumentException or RunTimeException. This will end the program and the JVM on a correct way
     * @param args String[] | all the arguments needed for the excecution. 
     */
    public FastqDemultiplex(String[] args) throws StopExcecutionException{
        this.parameters = new DemultiplexParameters();
        if (args.length == 0) {
            throw new IllegalArgumentException("\n\nNo arguments given.\n\n");
        }
        if (args[0].equals("version") || args[0].equals("-version") || args[0].equals("-v")){
            System.out.println(GBSdemultiplex.VERSION);
            System.out.println("This is " + GBSX.VERSION + ". A toolkit for experimental design and demultiplexing genotyping by \n" +
"sequencing experiments.");
            throw new StopExcecutionException();
        }
        
        if (args[0].equals("help") || args[0].equals("-help") || args[0].equals("-h")){
            System.out.println("");
            System.out.println("");
            System.out.println("This is the help of the demultiplexer.");
            System.out.println();
            System.out.println("KU Leuven");
            System.out.println("Licenced under " + GBSdemultiplex.LICENCE);
            System.out.println(GBSdemultiplex.VERSION);
            System.out.println();
            System.out.println("This program demultiplexes fastq or fastq.gz files direved of Sequencing with inline barcodes.");
            System.out.println("Like used in GBS, RAD, ... protocols.");
            System.out.println();
            System.out.println(this.parameters.getParametersHelp());
            System.out.println();
            System.out.println("Possible Standard Enzymes for the info file: (NAN is no enzyme)");
            for (EnzymeEnum enzyme : EnzymeEnum.values()){
                System.out.println("\t" + enzyme.getName());
            }
            System.out.println();
            System.out.println("Developed by the KU Leuven 2014");
            System.out.println("Licenced under " + GBSdemultiplex.LICENCE);
            System.out.println("For licence information use -licence");
            System.out.println();
            throw new StopExcecutionException();
        }
        
        if (args[0].equals("licence") || args[0].equals("-licence") || args[0].equals("-l")){
            System.out.println("This file is part of GBSX.\n" +
                "\n" +
                "GBSX is free software: you can redistribute it and/or modify\n" +
                "it under the terms of the GNU General Public License as published by\n" +
                "the Free Software Foundation, either version 3 of the License, or\n" +
                "(at your option) any later version.\n" +
                "\n" +
                "GBSX is distributed in the hope that it will be useful,\n" +
                "but WITHOUT ANY WARRANTY; without even the implied warranty of\n" +
                "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n" +
                "GNU General Public License for more details.\n" +
                "\n" +
                "You should have received a copy of the GNU General Public License\n" +
                "along with GBSX.  If not, see <http://www.gnu.org/licenses/>.");
            throw new StopExcecutionException();
        }
        
        //parse parameters
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];

            // An option.
            if (arg.startsWith("--") || arg.startsWith("-")) {
                //check if the option is invalid
                if (i + 1 >= args.length || args[i + 1].startsWith("-")) {
                    //invalid option
                    throw new RuntimeException("Value required for option "
                            + arg);
                }else{
                    //get the argument
                    DemultiplexArguments argument = DemultiplexArguments.INVALID_ARGUMENT.getArgument(arg);
                    this.parameters.setParameter(argument, args[++i]);
                }
            }
        }
        if (! this.parameters.areRequiredParametersSet()){
            throw new RuntimeException(this.parameters.getErrorRequiredParametersSet());
        }
        File outDir = new File(this.parameters.getOutputDirectory());
        if (! outDir.exists()){
            outDir.mkdirs();
        }
        //parse possible enzyme files
        if (this.parameters.mustAddEnzymes() || this.parameters.mustReplaceEnzymes()){
            String fileName;
            if (this.parameters.mustAddEnzymes()){
                fileName = this.parameters.getAddEnzymeFile();
            }else{
                fileName = this.parameters.getReplaceEnzymeFile();
            }
            try {
                EnzymeFileParser enzymeFileParser = new EnzymeFileParser();
                ArrayList<Enzyme> enzymeList = new ArrayList();
                enzymeList.addAll(enzymeFileParser.parseEnzymeFile(fileName));
                if (enzymeList.isEmpty()){
                    throw new RuntimeException("No valid enzymes found in the enzyme file");
                }
                if (this.parameters.mustAddEnzymes()){
                    this.parameters.getEnzymeCollection().addEnzymes(enzymeList);
                }else{
                    this.parameters.getEnzymeCollection().replaceEnzymes(enzymeList);
                }
            } catch (FileNotFoundException ex) {
                Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Enzyme file.", ex);
            } catch (IOException ex) {
                Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Enzyme file", ex);
            }
        }
        //parse the info file
        try {
            InfoFileParser infoFileParser = new InfoFileParser(this.parameters.getEnzymeCollection());
            this.sampleList = infoFileParser.parseInfoFile(this.parameters.getInfoFile());
            this.longestBarcodeLength = 0;
            for (Sample sample : this.sampleList){
                if (sample.getBarcode().length() > this.longestBarcodeLength){
                    this.longestBarcodeLength = sample.getBarcode().length();
                }
            }
            if (this.sampleList.isEmpty()){
                throw new RuntimeException("No valid info found in the info file");
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
            throw new RuntimeException("Couldn't open the info file.", ex);
        } catch (IOException ex) {
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
            throw new RuntimeException("Couldn't open the info file", ex);
        }
        try {
            //Open log file
            this.loggerFile = new LoggerFile(this.parameters);
        } catch (IOException ex) {
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
            throw new RuntimeException("Couldn't create a log file.", ex);
        } catch (ErrorInLogException ex){
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
            throw new RuntimeException("Couldn't create a log file.", ex);
        }
        this.findingDistanceAlgorithm = new FindingDistanceAlgorithm(this.parameters.getFindingsAlgorithm());
        this.correctionLog = new CorrectionLog(this.sampleList, this.loggerFile);
        //check if all are simple, or all are double
        int doubleBarcodeCount = 0;
        for (Sample sample : this.sampleList){
            if (sample.has2barcodes()){
                doubleBarcodeCount++;
            }
        }
        if (doubleBarcodeCount != 0 && doubleBarcodeCount < this.sampleList.size()){
            //if a mix of simple and double, throw an error
            throw new RuntimeException("Mix of single and double barcodes found.");
        }
        //double barcodes only possible with paired end sequencing
        if (doubleBarcodeCount == this.sampleList.size()){
            this.parameters.setDoubleBarcodes();
            if(this.parameters.getFastqFile2() == null){
                //single end demultiplex
            throw new RuntimeException("Use Double Barcodes, but no paired end fastq file found.");
            }
        }
    }

    
    /**
     * Start the processing of the demultiplex files
     * <br> can only if the constructor is executed right.
     * <br> the log file is updated and written.
     * @see FastqDemultiplex#processPairEndFastqFiles() 
     * @see FastqDemultiplex#processSingleReadFastqFiles() 
     * @see FastqDemultiplex#writeToLog(java.lang.String) 
     */
    public void processFastqFiles(){
        if (this.parameters.getFastqFile2() == null){
            this.processSingleReadFastqFiles();
        }else{
            this.processPairEndFastqFiles();
        }
        try {
            this.loggerFile.closeLog();
        } catch (ErrorInLogException ex) {
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    /**
     * the complete processing of pair end fastq files.
     * <br> Makes for every sample found in the info file 2 new fastq.gz file (R1 and R2)
     * <br> creates and keeps stats for the whole processing
     * <br> reads every read from both fastq file (a read is the description, the sequence and the quality => 4 lines in the file)
     * <br> parses the whole read to determinate the barcode, enzyme, complete sequence and complete quality
     * <br> parses the reverse read
     * <br> saves the parsed read and parsed reverse read to the corresponding sample, or to the undetermined file
     * <br> updates the stats
     * <br> at the end it closes all used files (the read fastq.gz, the new fastq.gz and the stats)
     * 
     * @see FastqBufferedReader
     * @see FastqBufferedWriter
     * @see DemultiplexStats
     * @see FastqDemultiplex#parseFastqRead(java.util.Map, java.util.Map) 
     */
    private void processPairEndFastqFiles(){
        //open the files and the reader
        File fastqFile1 = new File(this.parameters.getFastqFile1());
        File fastqFile2 = new File(this.parameters.getFastqFile2());
        try {
            //create buffered readers for the files
            FastqPairBufferedReader fastqReader = new FastqPairBufferedReader(fastqFile1, fastqFile2, this.parameters.mustBeZipped());
            
            //create the outputfiles
            HashMap<Sample, FastqPairBufferedWriter> sampleFiles = new HashMap();
            for (Sample sample : this.sampleList){
                String sampleFileR1name = this.parameters.getOutputDirectory() + System.getProperty("file.separator") + sample.getSampleID() + ".R1" + this.parameters.getFileExtension();
                String sampleFileR2name = this.parameters.getOutputDirectory() + System.getProperty("file.separator") + sample.getSampleID() + ".R2" + this.parameters.getFileExtension();
                if (this.parameters.useLongFileNames()){
                    sampleFileR1name = this.parameters.getOutputDirectory() + System.getProperty("file.separator") + sample.getSampleID() + "_" + sample.getBarcode() + "_" + sample.getEnzyme().getName().toLowerCase() + ".R1" + this.parameters.getFileExtension();
                    sampleFileR2name = this.parameters.getOutputDirectory() + System.getProperty("file.separator") + sample.getSampleID() + "_" + sample.getBarcode() + "_" + sample.getEnzyme().getName().toLowerCase() + ".R2" + this.parameters.getFileExtension();
                }
                sampleFiles.put(sample, new FastqPairBufferedWriter(new File(sampleFileR1name), new File(sampleFileR2name), this.parameters.mustBeZipped()));
            }
            FastqPairBufferedWriter undeterminedFastqFile = new FastqPairBufferedWriter(new File(this.parameters.getOutputDirectory() + System.getProperty("file.separator") + "undetermined" + ".R1" + this.parameters.getFileExtension()), 
                    new File(this.parameters.getOutputDirectory() + System.getProperty("file.separator") + "undetermined" + ".R2" + this.parameters.getFileExtension()), this.parameters.mustBeZipped());
            
            //create the stats file
            DemultiplexStats statsFile = new DemultiplexStats(this.sampleList, this.loggerFile, this.parameters.getFastqQualityScore());
            
            
            //init vars needed in the loop to go furter
            FastqRead fastq1 = null;
            FastqRead fastq2 = null;
            
            HashMap<String, FastqRead> readMap;
            int process_count = 0;
            while (((readMap = fastqReader.next()) != null) && ((fastq1 = readMap.get("r1")) != null) && ((fastq2 = readMap.get("r2")) != null)){
                process_count++;
                if (process_count % this.PROCESS_REPORTING_INTERVAL == 0){
                    System.out.println(process_count + " reads demultiplexed");
                }
                //read the next fastq line
                try {
                    //try to parse the fastq
                    ProcessedFragment newReads = this.parseFastqRead_new(fastq1, fastq2);
                    //if a read is empty: correct it
                    if (newReads.getRead1().getSequence().equals("")){
                        newReads = new ProcessedFragment(newReads.getSample(), new FastqRead(newReads.getRead1().getDescription(), "N", "#"), newReads.getRead2(), newReads.getMismatch());
                    }
                    if (newReads.getRead2().getSequence().equals("")){
                        newReads = new ProcessedFragment(newReads.getSample(), newReads.getRead1(), new FastqRead(newReads.getRead2().getDescription(), "N", "#"), newReads.getMismatch());
                    }
                    //check if a sequence must be rejected (to short)
                    if (newReads.getRead1().getSequence().length() < this.parameters.getMinimumSequenceLength()
                            || (! this.parameters.keepSequencesWithN() && newReads.getRead1().getSequence().contains("N"))){
                        statsFile.addRejectedRead(newReads.getSample());
                        undeterminedFastqFile.write(fastq1, fastq2);
                        statsFile.addUndeterminedStat();
                    }else{
                        //copies the read to the corresponding sample
                        sampleFiles.get(newReads.getSample()).write(newReads.getRead1(), newReads.getRead2());
                        //updates the stats
                        String totalQuality = newReads.getRead1().getQuality() + newReads.getRead2().getQuality();
                        statsFile.addStat(newReads.getSample(), newReads.getMismatch(), totalQuality);
                    }

                } catch (InvalidReadException ex) {
                    if (ex.getInvalidRead() != InvalidReadEnum.READ1 && ex.getInvalidRead() != InvalidReadEnum.READ2){
                        this.writeToLog("Unknown error in " + fastq1.getDescription() + " and " + fastq2.getDescription());
                    }
                    //write unknown fastq to the undetermined file
                    try{
                        undeterminedFastqFile.write(fastq1, fastq2);
                        statsFile.addUndeterminedStat();
                    } catch (IOException ioex) {
                        this.writeToLog("ERROR in writing the fastq files for " + fastq1.getDescription());
                        Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ioex);
                    }
                }
            }
            
            System.out.println(process_count + " reads demultiplexed");
            
            fastqReader.close();
            for (FastqPairBufferedWriter writer : sampleFiles.values()){
                writer.close();
            }
            undeterminedFastqFile.close();
            
            int mostMismatches = this.parameters.getAllowedMismatchesBarcode() * 2;
            for (Sample sample : this.sampleList){
                if (sample.getBarcodeMismatches() * 2 > mostMismatches){
                    mostMismatches = sample.getBarcodeMismatches();
                }
            }
            statsFile.saveStats(mostMismatches, this.parameters.getOutputDirectory());
            this.correctionLog.writeToFile(this.parameters.getOutputDirectory());
        } catch (IOException ex) {
            this.writeToLog("ERROR in reading the fastq files.");
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    /**
     * the complete processing of single read fastq files.
     * <br> Makes for every sample found in the info file a new fastq.gz file
     * <br> creates and keeps stats for the whole processing
     * <br> reads every read in the fastq file (a read is the description, the sequence and the quality => 4 lines in the file)
     * <br> parses the whole read to determinate the barcode, enzyme, complete sequence and complete quality
     * <br> saves the parsed read to the corresponding sample, or to the undetermined file
     * <br> updates the stats
     * <br> at the end it closes all used files (the read fastq.gz, the new fastq.gz and the stats)
     * 
     * @see FastqBufferedReader
     * @see FastqBufferedWriter
     * @see DemultiplexStats
     * @see FastqDemultiplex#parseFastqRead(java.util.Map) 
     */
    private void processSingleReadFastqFiles(){
        //open the files and the reader
        File fastqFile1 = new File(this.parameters.getFastqFile1());
        try {
            //create buffered readers for the files
            FastqBufferedReader fastq1Reader = new FastqBufferedReader(fastqFile1, this.parameters.mustBeZipped());
            
            //create the outputfiles
            HashMap<Sample, FastqBufferedWriter> sampleFiles = new HashMap();
            for (Sample sample : this.sampleList){
                String sampleFileName = this.parameters.getOutputDirectory() + System.getProperty("file.separator") + sample.getSampleID() + ".R1" + this.parameters.getFileExtension();
                if (this.parameters.useLongFileNames()){
                    sampleFileName = this.parameters.getOutputDirectory() + System.getProperty("file.separator") + sample.getSampleID() + "_" + sample.getBarcode() + "_" + sample.getEnzyme().getName().toLowerCase() + ".R1" + this.parameters.getFileExtension();
                }
                sampleFiles.put(sample, new FastqBufferedWriter(new File(sampleFileName), this.parameters.mustBeZipped()));
            }
            FastqBufferedWriter undeterminedFastqFile = new FastqBufferedWriter(new File(this.parameters.getOutputDirectory() + System.getProperty("file.separator") + "undetermined" + this.parameters.getFileExtension()), this.parameters.mustBeZipped());
            
            //create the stats file
            DemultiplexStats statsFile = new DemultiplexStats(this.sampleList, this.loggerFile, this.parameters.getFastqQualityScore());
            
            //init vars needed in the loop to go furter
            FastqRead fastq1;
            int process_count = 0;
            while ((fastq1 = fastq1Reader.next()) != null) {
                process_count++;
                if (process_count % this.PROCESS_REPORTING_INTERVAL == 0){
                    System.out.println(process_count + " reads demultiplexed");
                }
                //read the next fastq line
                try {
                    //try to parse the fastq
                    ProcessedFragment newReads = this.parseFastqRead(fastq1);
                    //check if a sequence is empty: correct it
                    if (newReads.getRead1().getSequence().equals("")){
                        newReads = new ProcessedFragment(newReads.getSample(), new FastqRead(newReads.getRead1().getDescription(), "N", "#"), newReads.getMismatch());
                    }
                    //check if a sequence must be rejected (to short)
                    if ((newReads.getRead1().getSequence().length() < this.parameters.getMinimumSequenceLength())
                            || (! this.parameters.keepSequencesWithN() && newReads.getRead1().getSequence().contains("N"))){
                        statsFile.addRejectedRead(newReads.getSample());
                        undeterminedFastqFile.write(fastq1);
                        statsFile.addUndeterminedStat();
                    }else{
                        //copies the read to the corresponding sample
                        sampleFiles.get(newReads.getSample()).write(newReads.getRead1());
                        //updates the stats
                        statsFile.addStat(newReads.getSample(), newReads.getMismatch(), newReads.getRead1().getQuality());
                    }

                } catch (InvalidReadException ex) {
                    //error int the read
                    if (ex.getInvalidRead() != InvalidReadEnum.READ1){
                        this.writeToLog("Unknown error in " + fastq1.getDescription());
                    }
                    //write unknown fastq to the undetermined file
                    try{
                        undeterminedFastqFile.write(fastq1);
                        statsFile.addUndeterminedStat();
                    } catch (IOException ioex) {
                        this.writeToLog("ERROR in writing the fastq files for " + fastq1.getDescription());
                        Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ioex);
                    }
                }  
            }
            
            System.out.println(process_count + " reads demultiplexed");
            
            
            //all reads are parsed => close all files
            fastq1Reader.close();
            for (FastqBufferedWriter writer : sampleFiles.values()){
                writer.close();
            }
            undeterminedFastqFile.close();
            int mostMismatches = this.parameters.getAllowedMismatchesBarcode();
            for (Sample sample : this.sampleList){
                if (sample.getBarcodeMismatches() > mostMismatches){
                    mostMismatches = sample.getBarcodeMismatches();
                }
            }
            statsFile.saveStats(mostMismatches, this.parameters.getOutputDirectory());
        } catch (IOException ex) {
            this.writeToLog("ERROR in reading the fastq files: " + ex.getMessage());
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    /**
     * 
     * @param log String | the log to be written to the logfile
     * @see LoggerFile#addToLog(java.lang.String) 
     */
    public void writeToLog(String log){
        try {
            this.loggerFile.addToLog(log + "\n");
        } catch (ErrorInLogException ex) {
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    /**
     * this method gets two reads from a pair-end fastq file. 
     * <br> Read1 is from the first read, Read2 is from the reverse read
     * <br> The reads are made of the 4 lines for every read: the information line, the sequence line, the + line and the quality line.
     * <br> The first read is searched for the barcode + enzyme site (looked to all known samples)
     * <br> if no barcode + enzyme (perfect match) is found, a InvalidReadException is thrown
     * <br> else these are trimed from the read (both from the sequence as the quality)
     * <br> then the read is searched for a cutsite of the same enzyme, if any that piece + the rest is removed (both from the sequence as the quality)
     * <br> The second read is searched for the enzyme site
     * <br> if no enzyme site (perfect match) is found, a InvalidReadException is thrown
     * <br> else these are trimed from the read (both sequence as quality)
     * <br> then the read is searched for a secund cutsite (this will be the complement of the enzyme site + barcode)
     * <br> if found these are removed + rest after (both from the sequence as the quality)
     * <br> a new array of Strings is created with 0 as the modified read of read1, and 1 as the modified read of read2
     * @param read1 Map<FastqParts, String> | the 4 lines for the first read: the information line, the sequence line, the + line and the quality line (as Fastq from BioJava)
     * @param read2 Map<FastqParts, String> | the 4 lines for the second read: the information line, the sequence line, the + line and the quality line (as Fastq from BioJava)
     * @return ProcessedRead | All information about the processed reads
     * @throws InvalidReadException if the first or second read doesn't has the right index (barcode + enzyme site for read1, enzyme site for read2)
     * @see FastqDemultiplex#findBestBarcode(java.lang.String) 
     * @see FastqDemultiplex#findRead1EnzymeLocation(java.lang.String, be.uzleuven.gc.logistics.gbsDemultiplex.model.Sample) 
     * @see ProcessedFragment
     */
    private ProcessedFragment parseFastqRead(FastqRead read1, FastqRead read2) throws InvalidReadException{
        //search for the optimal barcode
        //SampleBarcodeCombination sampleBarcodeCombination = this.findBestBarcode(read1.getSequence());
        SampleBarcodeCombination sampleBarcodeCombination = this.findGBSBarcode(read1.getSequence(), this.parameters.getStartDistance());
        if (sampleBarcodeCombination == null){
            //no optimal barcode found in the first read
            throw new InvalidReadException(read1, read2, InvalidReadEnum.READ1);
        }
        Sample sample = sampleBarcodeCombination.getSample();
        String barcodeEnzyme = sampleBarcodeCombination.getSample().getBarcode();
        String enzymeCutsite = sampleBarcodeCombination.getEnzymeCutsite();
        int barcodeEnzymeLength = sampleBarcodeCombination.getLengthFoundBarcode();
        if (! this.parameters.keepCutSites()){
            //the cutsite mustn't be kept
            barcodeEnzyme += sampleBarcodeCombination.getEnzymeCutsite();
            barcodeEnzymeLength += sampleBarcodeCombination.getLengthFoundEnzyme();
        }
        //remove the barcode and the enzyme site
        int read1BarcodeLocation = sampleBarcodeCombination.getLocation();
        String read1modifiedSequence = read1.getSequence().substring(read1BarcodeLocation + barcodeEnzymeLength, read1.getSequence().length() - (this.longestBarcodeLength - sample.getBarcode().length()));
        String read1modifiedQuality = read1.getQuality().substring(read1BarcodeLocation + barcodeEnzymeLength, read1.getSequence().length() - (this.longestBarcodeLength - sample.getBarcode().length()));
        
        if (this.parameters.keepCutSites()){
            //add the cutsite to the barcodeEnzyme when not already added (to have the correct complement)
            barcodeEnzyme += sampleBarcodeCombination.getEnzymeCutsite();
        }
        
        //find the next enzyme site (if there is any)
        int read1EndLocation = -1;
        int[] read1EndLocationLength = {-1, 0};
//        EnzymeComparator enzymeComparator = new EnzymeComparator();
//        if (enzymeComparator.compare(sample.getEnzyme(), this.parameters.getNeutralEnzyme()) != 0){
            read1EndLocationLength = this.findRead1EnzymeLocation(read1modifiedSequence, sample);
            read1EndLocation = read1EndLocationLength[0];
//        }
        String read1optimalSequence = read1modifiedSequence;
        String read1optimalQuality = read1modifiedQuality;
        if (read1EndLocation != -1){
            //compliment barcode found
            if (this.parameters.keepCutSites() && ! this.parameters.isRadData()){
                //cutsites must be kept
                read1EndLocation += read1EndLocationLength[1];
            }
            read1optimalSequence = read1modifiedSequence.substring(0, read1EndLocation);
            read1optimalQuality = read1modifiedQuality.substring(0, read1EndLocation);
        }
        
        
        //parsing of read2
        //find the first enzyme location
        
        //just cut the first basepairs because these will be the enzyme site
        int read2firstEnzymeLocation = 0;
        String read2modifiedSequence = read2.getSequence().substring(read2firstEnzymeLocation);
        String read2modifiedQuality = read2.getQuality().substring(read2firstEnzymeLocation);
        if (this.parameters.isRadData()){
            //RAD data
//            if (this.findingDistanceAlgorithm.isEquivalent(read2modifiedSequence.substring(0, sample.getBarcode().length()), sample.getBarcode(), this.parameters.getAllowedMismatchesBarcode(sample))){
//                //possible barcode found
//                read2modifiedSequence = read2modifiedSequence.substring(sample.getBarcode().length());
//                read2modifiedQuality = read2modifiedQuality.substring(sample.getBarcode().length());
//                if (! this.parameters.keepCutSites()){
//                    //remove the enzyme site
//                    String foundEnzyme = "";
//                    for (String enzymeSite : sample.getEnzyme().getInitialCutSiteRemnant()){
//                        if (this.findingDistanceAlgorithm.isEquivalent(read2modifiedSequence.substring(0, enzymeSite.length()), enzymeSite, this.parameters.getAllowedMismatchesEnzyme())){
//                            foundEnzyme = enzymeSite;
//                        }
//                    }
//                    read2modifiedSequence = read2modifiedSequence.substring(foundEnzyme.length());
//                    read2modifiedQuality = read2modifiedQuality.substring(foundEnzyme.length());
//                }
//            }
        }else{
            //GBS data
            if (! this.parameters.keepCutSites()){
                //remove the enzyme site
                String foundEnzyme = "";
                for (String enzymeSite : sample.getEnzyme().getInitialCutSiteRemnant()){
                    if (this.findingDistanceAlgorithm.isEquivalent(read2modifiedSequence.substring(0, enzymeSite.length()), enzymeSite, this.parameters.getAllowedMismatchesEnzyme())){
                        foundEnzyme = enzymeSite;
                    }
                }
                read2modifiedSequence = read2modifiedSequence.substring(foundEnzyme.length());
                read2modifiedQuality = read2modifiedQuality.substring(foundEnzyme.length());
            }
        }
        
        //find the next enzyme site (if there is any)
        //String complementBarcodeEnzyme = BasePair.getComplementSequence(barcodeEnzyme);
        int[] read2secondEnzymeLocationLength = this.findRead2EnzymeLocation(read2modifiedSequence, sample, enzymeCutsite);
        String read2optimalSequence = read2modifiedSequence;
        String read2optimalQuality = read2modifiedQuality;
        if (read2secondEnzymeLocationLength[0] != -1){
            int read2secondEnzymeLocation = read2secondEnzymeLocationLength[0];
            //enzyme site found
            if (this.parameters.keepCutSites()){
                //keep the enzyme sites
                read2secondEnzymeLocation += read2secondEnzymeLocationLength[1];
            }
            read2optimalSequence = read2modifiedSequence.substring(0, read2secondEnzymeLocation);
            read2optimalQuality = read2modifiedQuality.substring(0, read2secondEnzymeLocation);
        }
        
        
        //size selection
        boolean trimedR1 = false;
        boolean trimedR2 = false;
        int lengthR1 = read1optimalSequence.length();
        int lengthR2 = read2optimalSequence.length();
        String sequenceError = "TRIM\t" + "trimok" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
        
        
        if (read1optimalSequence.length() != read1.getSequence().length() - this.longestBarcodeLength){
            trimedR1 = true;
        }
        if (read2optimalSequence.length() != read2.getSequence().length()){
            trimedR2 = true;
        }
        
        if (! trimedR1 && ! trimedR2){
            //both original => ok
            this.correctionLog.addCorrecterTrimOk(sample);
        }else if (! trimedR1 && read2optimalSequence.length() >= read1optimalSequence.length()){
            //R1 is original, R2 is trimed, but same length or longer => ok
            this.correctionLog.addCorrecterTrimOk(sample);
        }else if (read1optimalSequence.length() == read2optimalSequence.length()){
            //R1 and R2 are same length => OK
            this.correctionLog.addCorrecterTrimOk(sample);
        }else{
            int compareLength = this.longestBarcodeLength;
            if (this.parameters.keepCutSites()){
                compareLength += sample.getPossibleEnzymeCutSiteLength();
            }
            //R1 and R2 have different sizes, find lowest and check
            if (read1optimalSequence.length() < read2optimalSequence.length()){
                //R1 is shortest
                if (read1optimalSequence.length() > compareLength){
                    String read1end = read1optimalSequence.substring(read1optimalSequence.length() - compareLength, read1optimalSequence.length());
                    String expectedStartR2 = BasePair.getComplementSequence(read1end);
                    if (this.parameters.keepCutSites()){
                        if (this.findingDistanceAlgorithm.isEquivalent(read2optimalSequence.substring(0, sample.getPossibleEnzymeCutSiteLength()), expectedStartR2.substring(0, sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesEnzyme())
                                && this.findingDistanceAlgorithm.isEquivalent(read2optimalSequence.substring(sample.getPossibleEnzymeCutSiteLength(), compareLength), expectedStartR2.substring(sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesBarcode(sample))){
                            //is equivalent:(read2(0-cutsite), complement read1(0-cutsite) with allowed mismatches enzyme)
                            //and is equivalent:(read2(cutsite-comparelength), complement read1(cutsite-comparelength) with allowed mismatches barcode)
                            //is same => trim R2
                            read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
                            read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
                            trimedR2 = true;
                            //sequenceError = "CORRECTED R2\n";
                            this.correctionLog.addCorrecterR2Corrected(sample);
                        }else{
                            if (read1optimalSequence.length() + sample.getBarcode().length() + this.parameters.getAdaptorCompareSize() >= read1.getSequence().length()){
                                //read1 is only checked on cutsite, not on adaptor, not corrected so wrong
                                int minus = this.longestBarcodeLength - sample.getBarcode().length();
                                read1optimalSequence = read1.getSequence().substring(sample.getBarcode().length(), read1.getSequence().length() - minus);
                                read1optimalQuality = read1.getQuality().substring(sample.getBarcode().length(), read1.getSequence().length() - minus);
                                trimedR1 = false;
                                this.correctionLog.addCorrecterR1Corrected(sample);
                            }else{
                                this.correctionLog.addCorrecterR2NotCorrected(sample);
                                //sequenceError = "NOTCORRECT R2\n";
                            }
                        }
                    }else{
                        //not keep cutsites
                        if (this.findingDistanceAlgorithm.isEquivalent(read2optimalSequence.substring(0, compareLength), expectedStartR2, this.parameters.getAllowedMismatchesBarcode(sample))){
                            //is same => trim R2
                            read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
                            read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
                            trimedR2 = true;
                            this.correctionLog.addCorrecterR2Corrected(sample);
                            //sequenceError = "CORRECTED R2\n";
                        }else{
                            if (read1optimalSequence.length() + sample.getBarcode().length() + this.parameters.getAdaptorCompareSize() + sample.getPossibleEnzymeCutSiteLength() + sample.getPossibleEnzymeCutSiteLength() >= read1.getSequence().length()){
                                //read1 is only checked on cutsite, not on adaptor, not corrected so wrong
                                int minus = this.longestBarcodeLength - sample.getBarcode().length();
                                read1optimalSequence = read1.getSequence().substring(sample.getBarcode().length() + sample.getPossibleEnzymeCutSiteLength(), read1.getSequence().length() - minus);
                                read1optimalQuality = read1.getQuality().substring(sample.getBarcode().length() + sample.getPossibleEnzymeCutSiteLength(), read1.getSequence().length() - minus);
                                trimedR1 = false;
                                this.correctionLog.addCorrecterR1Corrected(sample);
                            }else{
                                this.correctionLog.addCorrecterR2NotCorrected(sample);
                                //sequenceError = "NOTCORRECT R2\n";
                            }
                        }
                    }
                }
            }else if (read1optimalSequence.length() > read2optimalSequence.length()){
                //R2 is shortest
                if (read2optimalSequence.length() > compareLength){
                    String read2end = read2optimalSequence.substring(read2optimalSequence.length() - compareLength, read2optimalSequence.length());
                    String expectedStartR1 = BasePair.getComplementSequence(read2end);
                    if (this.parameters.keepCutSites()){
                        if (this.findingDistanceAlgorithm.isEquivalent(read1optimalSequence.substring(0, sample.getPossibleEnzymeCutSiteLength()), expectedStartR1.substring(0, sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesEnzyme())
                                && this.findingDistanceAlgorithm.isEquivalent(read1optimalSequence.substring(sample.getPossibleEnzymeCutSiteLength(), compareLength), expectedStartR1.substring(sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesBarcode(sample))){
                            read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                            read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                            trimedR1 = true;
                            this.correctionLog.addCorrecterR1Corrected(sample);
                            //sequenceError = "CORRECTED R1\n";
                        }else{
                            this.correctionLog.addCorrecterR1NotCorrected(sample);
                            //sequenceError = "NOTCORRECT R1\n";
                        }
                    }else{
                        //not keep cutsites
                        if (this.findingDistanceAlgorithm.isEquivalent(read1optimalSequence.substring(0, compareLength), expectedStartR1, this.parameters.getAllowedMismatchesBarcode(sample))){
                            read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                            read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                            trimedR1 = true;
                            this.correctionLog.addCorrecterR1Corrected(sample);
                            //sequenceError = "CORRECTED R1\n";
                        }else{
                            this.correctionLog.addCorrecterR1NotCorrected(sample);
                            //sequenceError = "NOTCORRECT R1\n";
                        }
                    }
                }
            }else{
                //same size (may not occure here)
                this.correctionLog.addCorrecterTrimOk(sample);
            }
        }
        
        
        
        
        if (trimedR1 && ! trimedR2){
            //R1 was trimed, but R2 not
            int mismatch = -1;
            if (read1optimalSequence.length() < (read1.getSequence().length() - this.longestBarcodeLength - sample.getBarcode().length() + 1)){
                mismatch = (MismatchIndelDistance.calculateEquivalentDistance(read2optimalSequence.substring(read1optimalSequence.length() + 1, (read1optimalSequence.length() + 1 + sample.getBarcode().length())), sample.getComplementBarcode(), 1)[0]);
            }
            read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
            read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
            sequenceError = "TRIM\t" + "tR1nR2" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2) + "\t" + mismatch;
            this.correctionLog.addTrimTrimR1NotR2Fail(sample);
        }else if (! trimedR1 && trimedR2){
            //R2 was trimed, not R1, so check sizes to diside to trim
            if (read1optimalSequence.length() > read2optimalSequence.length()){
                //read 1 is longer, so trim read 1
                int mismatch = -1;
                if (read2optimalSequence.length() < read1.getSequence().length() - this.parameters.getAdaptorCompareSize() - this.longestBarcodeLength){
                    mismatch = (MismatchIndelDistance.calculateEquivalentDistance(read1optimalSequence.substring(read2optimalSequence.length(), (read2optimalSequence.length() + this.parameters.getAdaptorCompareSize())), this.parameters.getCommonAdaptor(), 1)[0]);;
                }
                read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                sequenceError = "TRIM\t" + "nR1tR2not" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2) + "\t" + mismatch;
                this.correctionLog.addTrimNotR1TrimR2butFail(sample);
            }else if (read1optimalSequence.length() < read2optimalSequence.length()){
                //read 2 is longer, so is ok
                sequenceError = "TRIM\t" + "nR1tR2okl" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
                this.correctionLog.addTrimNotR1TrimR2butOk(sample);
            }else{
                //same length so ok
                sequenceError = "TRIM\t" + "nR1tR2oks" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
                this.correctionLog.addTrimNotR1TrimR2butOk(sample);
            }
        }else if (trimedR1 && trimedR2){
            //both are trimed so both have to be the same size
            if (read1optimalSequence.length() > read2optimalSequence.length()){
                //read 1 is longer, so trim read 1
                read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                sequenceError = "TRIM\t" + "tR1tR2notR1" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2);
                this.correctionLog.addTrimTrimR1TrimR2longR1(sample);
            }else if (read1optimalSequence.length() < read2optimalSequence.length()){
                //read 2 is longer, so trim read 2
                read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
                read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
                sequenceError = "TRIM\t" + "tR1tR2notR2" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2);
                this.correctionLog.addTrimTrimR1TrimR2longR2(sample);
            }else{
                //both reads have same length so ok
                sequenceError = "TRIM\t" + "tR1tR2ok" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
                this.correctionLog.addTrimTrimR1TrimR2ok(sample);
            }
        }else if (! trimedR1 && ! trimedR2){
            //perfect sequence, not trimmed
            this.correctionLog.addTrimNotR1NotR2(sample);
        }
        
        
        //save results
        FastqRead resultRead1 = new FastqRead(read1.getDescription(), read1optimalSequence, read1optimalQuality);
        FastqRead resultRead2 = new FastqRead(read2.getDescription(), read2optimalSequence, read2optimalQuality);
        return new ProcessedFragment(sample, resultRead1, resultRead2, sampleBarcodeCombination.getMismatches(), sequenceError);
    }
    
    /**
     * this method gets two reads from a pair-end fastq file. 
     * <br> Read1 is from the first read, Read2 is from the reverse read
     * <br> The reads are made of the 4 lines for every read: the information line, the sequence line, the + line and the quality line.
     * <br> The first read is searched for the barcode + enzyme site (looked to all known samples)
     * <br> if no barcode + enzyme (perfect match) is found, a InvalidReadException is thrown
     * <br> else these are trimed from the read (both from the sequence as the quality)
     * <br> then the read is searched for a cutsite of the same enzyme, if any that piece + the rest is removed (both from the sequence as the quality)
     * <br> The second read is searched for the enzyme site
     * <br> if no enzyme site (perfect match) is found, a InvalidReadException is thrown
     * <br> else these are trimed from the read (both sequence as quality)
     * <br> then the read is searched for a secund cutsite (this will be the complement of the enzyme site + barcode)
     * <br> if found these are removed + rest after (both from the sequence as the quality)
     * <br> a new array of Strings is created with 0 as the modified read of read1, and 1 as the modified read of read2
     * @param read1 Map<FastqParts, String> | the 4 lines for the first read: the information line, the sequence line, the + line and the quality line (as Fastq from BioJava)
     * @param read2 Map<FastqParts, String> | the 4 lines for the second read: the information line, the sequence line, the + line and the quality line (as Fastq from BioJava)
     * @return ProcessedRead | All information about the processed reads
     * @throws InvalidReadException if the first or second read doesn't has the right index (barcode + enzyme site for read1, enzyme site for read2)
     * @see FastqDemultiplex#findBestBarcode(java.lang.String) 
     * @see FastqDemultiplex#findRead1EnzymeLocation(java.lang.String, be.uzleuven.gc.logistics.gbsDemultiplex.model.Sample) 
     * @see ProcessedFragment
     */
    private ProcessedFragment parseFastqRead_new(FastqRead read1, FastqRead read2) throws InvalidReadException{
        //search for the optimal barcode
        SampleBarcodeCombination[] comb = this.findGBSBarcode2(read1.getSequence(), this.parameters.getStartDistance(), read2.getSequence());
        if (comb == null){
            //no optimal barcode found in the first read or both reads
            throw new InvalidReadException(read1, read2, InvalidReadEnum.READ1);
        }
        SampleBarcodeCombination sampleBarcodeCombination1 = comb[0];
        SampleBarcodeCombination sampleBarcodeCombination2 = null;
        if (comb.length == 2){
            sampleBarcodeCombination2 = comb[1];
        }
        Sample sample = sampleBarcodeCombination1.getSample();
        String barcodeEnzyme = sampleBarcodeCombination1.getSample().getBarcode();
        String enzymeCutsite = sampleBarcodeCombination1.getEnzymeCutsite();
        int barcodeEnzymeLength = sampleBarcodeCombination1.getLengthFoundBarcode();
        if (! this.parameters.keepCutSites()){
            //the cutsite mustn't be kept
            barcodeEnzyme += sampleBarcodeCombination1.getEnzymeCutsite();
            barcodeEnzymeLength += sampleBarcodeCombination1.getLengthFoundEnzyme();
        }
        //remove the barcode and the enzyme site
        int read1BarcodeLocation = sampleBarcodeCombination1.getLocation();
        String read1modifiedSequence = read1.getSequence().substring(read1BarcodeLocation + barcodeEnzymeLength, read1.getSequence().length() - (this.longestBarcodeLength - sample.getBarcode().length()));
        String read1modifiedQuality = read1.getQuality().substring(read1BarcodeLocation + barcodeEnzymeLength, read1.getSequence().length() - (this.longestBarcodeLength - sample.getBarcode().length()));
        
        if (this.parameters.keepCutSites()){
            //add the cutsite to the barcodeEnzyme when not already added (to have the correct complement)
            barcodeEnzyme += sampleBarcodeCombination1.getEnzymeCutsite();
        }
        
        //find the next enzyme site (if there is any)
        int read1EndLocation = -1;
        int[] read1EndLocationLength = {-1, 0};
//        EnzymeComparator enzymeComparator = new EnzymeComparator();
//        if (enzymeComparator.compare(sample.getEnzyme(), this.parameters.getNeutralEnzyme()) != 0){
            read1EndLocationLength = this.findRead1EnzymeLocation(read1modifiedSequence, sample);
            read1EndLocation = read1EndLocationLength[0];
//        }
        String read1optimalSequence = read1modifiedSequence;
        String read1optimalQuality = read1modifiedQuality;
        if (read1EndLocation != -1){
            //compliment barcode found
            if (this.parameters.keepCutSites() && ! this.parameters.isRadData()){
                //cutsites must be kept
                read1EndLocation += read1EndLocationLength[1];
            }
            read1optimalSequence = read1modifiedSequence.substring(0, read1EndLocation);
            read1optimalQuality = read1modifiedQuality.substring(0, read1EndLocation);
        }
        
        
        //parsing of read2
        //find the first enzyme location
        
        //just cut the first basepairs because these will be the enzyme site
        int read2firstEnzymeLocation = 0;
        if (sampleBarcodeCombination2 != null){
            read2firstEnzymeLocation = sampleBarcodeCombination2.getLengthFoundBarcode();
        }
        String read2modifiedSequence = read2.getSequence().substring(read2firstEnzymeLocation);
        String read2modifiedQuality = read2.getQuality().substring(read2firstEnzymeLocation);
        if (this.parameters.isRadData()){
            //RAD data
//            if (this.findingDistanceAlgorithm.isEquivalent(read2modifiedSequence.substring(0, sample.getBarcode().length()), sample.getBarcode(), this.parameters.getAllowedMismatchesBarcode(sample))){
//                //possible barcode found
//                read2modifiedSequence = read2modifiedSequence.substring(sample.getBarcode().length());
//                read2modifiedQuality = read2modifiedQuality.substring(sample.getBarcode().length());
//                if (! this.parameters.keepCutSites()){
//                    //remove the enzyme site
//                    String foundEnzyme = "";
//                    for (String enzymeSite : sample.getEnzyme().getInitialCutSiteRemnant()){
//                        if (this.findingDistanceAlgorithm.isEquivalent(read2modifiedSequence.substring(0, enzymeSite.length()), enzymeSite, this.parameters.getAllowedMismatchesEnzyme())){
//                            foundEnzyme = enzymeSite;
//                        }
//                    }
//                    read2modifiedSequence = read2modifiedSequence.substring(foundEnzyme.length());
//                    read2modifiedQuality = read2modifiedQuality.substring(foundEnzyme.length());
//                }
//            }
        }else{
            //GBS data
            if (! this.parameters.keepCutSites()){
                //remove the enzyme site
                String foundEnzyme = "";
                for (String enzymeSite : sample.getEnzyme().getInitialCutSiteRemnant()){
                    if (this.findingDistanceAlgorithm.isEquivalent(read2modifiedSequence.substring(0, enzymeSite.length()), enzymeSite, this.parameters.getAllowedMismatchesEnzyme())){
                        foundEnzyme = enzymeSite;
                    }
                }
                read2modifiedSequence = read2modifiedSequence.substring(foundEnzyme.length());
                read2modifiedQuality = read2modifiedQuality.substring(foundEnzyme.length());
            }
        }
        
        //find the next enzyme site (if there is any)
        //String complementBarcodeEnzyme = BasePair.getComplementSequence(barcodeEnzyme);
        int[] read2secondEnzymeLocationLength = this.findRead2EnzymeLocation(read2modifiedSequence, sample, enzymeCutsite);
        String read2optimalSequence = read2modifiedSequence;
        String read2optimalQuality = read2modifiedQuality;
        if (read2secondEnzymeLocationLength[0] != -1){
            int read2secondEnzymeLocation = read2secondEnzymeLocationLength[0];
            //enzyme site found
            if (this.parameters.keepCutSites()){
                //keep the enzyme sites
                read2secondEnzymeLocation += read2secondEnzymeLocationLength[1];
            }
            read2optimalSequence = read2modifiedSequence.substring(0, read2secondEnzymeLocation);
            read2optimalQuality = read2modifiedQuality.substring(0, read2secondEnzymeLocation);
        }
        
        
        //size selection
        boolean trimedR1 = false;
        boolean trimedR2 = false;
        int lengthR1 = read1optimalSequence.length();
        int lengthR2 = read2optimalSequence.length();
        String sequenceError = "TRIM\t" + "trimok" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
        
        
        if (read1optimalSequence.length() != read1.getSequence().length() - this.longestBarcodeLength){
            trimedR1 = true;
        }
        if (read2optimalSequence.length() != read2.getSequence().length()){
            trimedR2 = true;
        }
        
        if (! trimedR1 && ! trimedR2){
            //both original => ok
            this.correctionLog.addCorrecterTrimOk(sample);
        }else if (! trimedR1 && read2optimalSequence.length() >= read1optimalSequence.length()){
            //R1 is original, R2 is trimed, but same length or longer => ok
            this.correctionLog.addCorrecterTrimOk(sample);
        }else if (read1optimalSequence.length() == read2optimalSequence.length()){
            //R1 and R2 are same length => OK
            this.correctionLog.addCorrecterTrimOk(sample);
        }else{
            int compareLength = this.longestBarcodeLength;
            if (this.parameters.keepCutSites()){
                compareLength += sample.getPossibleEnzymeCutSiteLength();
            }
            //R1 and R2 have different sizes, find lowest and check
            if (read1optimalSequence.length() < read2optimalSequence.length()){
                //R1 is shortest
                if (read1optimalSequence.length() > compareLength){
                    String read1end = read1optimalSequence.substring(read1optimalSequence.length() - compareLength, read1optimalSequence.length());
                    String expectedStartR2 = BasePair.getComplementSequence(read1end);
                    if (this.parameters.keepCutSites()){
                        if (this.findingDistanceAlgorithm.isEquivalent(read2optimalSequence.substring(0, sample.getPossibleEnzymeCutSiteLength()), expectedStartR2.substring(0, sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesEnzyme())
                                && this.findingDistanceAlgorithm.isEquivalent(read2optimalSequence.substring(sample.getPossibleEnzymeCutSiteLength(), compareLength), expectedStartR2.substring(sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesBarcode(sample))){
                            //is equivalent:(read2(0-cutsite), complement read1(0-cutsite) with allowed mismatches enzyme)
                            //and is equivalent:(read2(cutsite-comparelength), complement read1(cutsite-comparelength) with allowed mismatches barcode)
                            //is same => trim R2
                            read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
                            read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
                            trimedR2 = true;
                            //sequenceError = "CORRECTED R2\n";
                            this.correctionLog.addCorrecterR2Corrected(sample);
                        }else{
                            if (read1optimalSequence.length() + sample.getBarcode().length() + this.parameters.getAdaptorCompareSize() >= read1.getSequence().length()){
                                //read1 is only checked on cutsite, not on adaptor, not corrected so wrong
                                int minus = this.longestBarcodeLength - sample.getBarcode().length();
                                read1optimalSequence = read1.getSequence().substring(sample.getBarcode().length(), read1.getSequence().length() - minus);
                                read1optimalQuality = read1.getQuality().substring(sample.getBarcode().length(), read1.getSequence().length() - minus);
                                trimedR1 = false;
                                this.correctionLog.addCorrecterR1Corrected(sample);
                            }else{
                                this.correctionLog.addCorrecterR2NotCorrected(sample);
                                //sequenceError = "NOTCORRECT R2\n";
                            }
                        }
                    }else{
                        //not keep cutsites
                        if (this.findingDistanceAlgorithm.isEquivalent(read2optimalSequence.substring(0, compareLength), expectedStartR2, this.parameters.getAllowedMismatchesBarcode(sample))){
                            //is same => trim R2
                            read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
                            read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
                            trimedR2 = true;
                            this.correctionLog.addCorrecterR2Corrected(sample);
                            //sequenceError = "CORRECTED R2\n";
                        }else{
                            if (read1optimalSequence.length() + sample.getBarcode().length() + this.parameters.getAdaptorCompareSize() + sample.getPossibleEnzymeCutSiteLength() + sample.getPossibleEnzymeCutSiteLength() >= read1.getSequence().length()){
                                //read1 is only checked on cutsite, not on adaptor, not corrected so wrong
                                int minus = this.longestBarcodeLength - sample.getBarcode().length();
                                read1optimalSequence = read1.getSequence().substring(sample.getBarcode().length() + sample.getPossibleEnzymeCutSiteLength(), read1.getSequence().length() - minus);
                                read1optimalQuality = read1.getQuality().substring(sample.getBarcode().length() + sample.getPossibleEnzymeCutSiteLength(), read1.getSequence().length() - minus);
                                trimedR1 = false;
                                this.correctionLog.addCorrecterR1Corrected(sample);
                            }else{
                                this.correctionLog.addCorrecterR2NotCorrected(sample);
                                //sequenceError = "NOTCORRECT R2\n";
                            }
                        }
                    }
                }
            }else if (read1optimalSequence.length() > read2optimalSequence.length()){
                //R2 is shortest
                if (read2optimalSequence.length() > compareLength){
                    String read2end = read2optimalSequence.substring(read2optimalSequence.length() - compareLength, read2optimalSequence.length());
                    String expectedStartR1 = BasePair.getComplementSequence(read2end);
                    if (this.parameters.keepCutSites()){
                        if (this.findingDistanceAlgorithm.isEquivalent(read1optimalSequence.substring(0, sample.getPossibleEnzymeCutSiteLength()), expectedStartR1.substring(0, sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesEnzyme())
                                && this.findingDistanceAlgorithm.isEquivalent(read1optimalSequence.substring(sample.getPossibleEnzymeCutSiteLength(), compareLength), expectedStartR1.substring(sample.getPossibleEnzymeCutSiteLength()), this.parameters.getAllowedMismatchesBarcode(sample))){
                            read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                            read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                            trimedR1 = true;
                            this.correctionLog.addCorrecterR1Corrected(sample);
                            //sequenceError = "CORRECTED R1\n";
                        }else{
                            this.correctionLog.addCorrecterR1NotCorrected(sample);
                            //sequenceError = "NOTCORRECT R1\n";
                        }
                    }else{
                        //not keep cutsites
                        if (this.findingDistanceAlgorithm.isEquivalent(read1optimalSequence.substring(0, compareLength), expectedStartR1, this.parameters.getAllowedMismatchesBarcode(sample))){
                            read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                            read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                            trimedR1 = true;
                            this.correctionLog.addCorrecterR1Corrected(sample);
                            //sequenceError = "CORRECTED R1\n";
                        }else{
                            this.correctionLog.addCorrecterR1NotCorrected(sample);
                            //sequenceError = "NOTCORRECT R1\n";
                        }
                    }
                }
            }else{
                //same size (may not occure here)
                this.correctionLog.addCorrecterTrimOk(sample);
            }
        }
        
        
        
        
        if (trimedR1 && ! trimedR2){
            //R1 was trimed, but R2 not
            int mismatch = -1;
            if (read1optimalSequence.length() < (read1.getSequence().length() - this.longestBarcodeLength - sample.getBarcode().length() + 1)){
                mismatch = (MismatchIndelDistance.calculateEquivalentDistance(read2optimalSequence.substring(read1optimalSequence.length() + 1, (read1optimalSequence.length() + 1 + sample.getBarcode().length())), sample.getComplementBarcode(), 1)[0]);
            }
            read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
            read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
            sequenceError = "TRIM\t" + "tR1nR2" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2) + "\t" + mismatch;
            this.correctionLog.addTrimTrimR1NotR2Fail(sample);
        }else if (! trimedR1 && trimedR2){
            //R2 was trimed, not R1, so check sizes to diside to trim
            if (read1optimalSequence.length() > read2optimalSequence.length()){
                //read 1 is longer, so trim read 1
                int mismatch = -1;
                if (read2optimalSequence.length() < read1.getSequence().length() - this.parameters.getAdaptorCompareSize() - this.longestBarcodeLength){
                    mismatch = (MismatchIndelDistance.calculateEquivalentDistance(read1optimalSequence.substring(read2optimalSequence.length(), (read2optimalSequence.length() + this.parameters.getAdaptorCompareSize())), this.parameters.getCommonAdaptor(), 1)[0]);;
                }
                read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                sequenceError = "TRIM\t" + "nR1tR2not" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2) + "\t" + mismatch;
                this.correctionLog.addTrimNotR1TrimR2butFail(sample);
            }else if (read1optimalSequence.length() < read2optimalSequence.length()){
                //read 2 is longer, so is ok
                sequenceError = "TRIM\t" + "nR1tR2okl" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
                this.correctionLog.addTrimNotR1TrimR2butOk(sample);
            }else{
                //same length so ok
                sequenceError = "TRIM\t" + "nR1tR2oks" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
                this.correctionLog.addTrimNotR1TrimR2butOk(sample);
            }
        }else if (trimedR1 && trimedR2){
            //both are trimed so both have to be the same size
            if (read1optimalSequence.length() > read2optimalSequence.length()){
                //read 1 is longer, so trim read 1
                read1optimalSequence = read1optimalSequence.substring(0, read2optimalSequence.length());
                read1optimalQuality = read1optimalQuality.substring(0, read2optimalSequence.length());
                sequenceError = "TRIM\t" + "tR1tR2notR1" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2);
                this.correctionLog.addTrimTrimR1TrimR2longR1(sample);
            }else if (read1optimalSequence.length() < read2optimalSequence.length()){
                //read 2 is longer, so trim read 2
                read2optimalSequence = read2optimalSequence.substring(0, read1optimalSequence.length());
                read2optimalQuality = read2optimalQuality.substring(0, read1optimalSequence.length());
                sequenceError = "TRIM\t" + "tR1tR2notR2" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length() + "\t" + (lengthR1 - lengthR2);
                this.correctionLog.addTrimTrimR1TrimR2longR2(sample);
            }else{
                //both reads have same length so ok
                sequenceError = "TRIM\t" + "tR1tR2ok" + "\t" + read1.getDescription() + "\t" + read2.getDescription() + "\tor1:" + lengthR1 + "\tor2:" + lengthR2 + "\tread1:" + read1optimalSequence.length() + "\tread2:" + read2optimalSequence.length();
                this.correctionLog.addTrimTrimR1TrimR2ok(sample);
            }
        }else if (! trimedR1 && ! trimedR2){
            //perfect sequence, not trimmed
            this.correctionLog.addTrimNotR1NotR2(sample);
        }
        
        
        //save results
        FastqRead resultRead1 = new FastqRead(read1.getDescription(), read1optimalSequence, read1optimalQuality);
        FastqRead resultRead2 = new FastqRead(read2.getDescription(), read2optimalSequence, read2optimalQuality);
        int mismatches = sampleBarcodeCombination1.getMismatches();
        if (sampleBarcodeCombination2 != null){
            mismatches += sampleBarcodeCombination2.getMismatches();
        }
        return new ProcessedFragment(sample, resultRead1, resultRead2, mismatches, sequenceError);
    }
     
    /**
     * this method gets two reads from a single read fastq file. 
     * <br> Read1 is from the first read
     * <br> The reads are made of the 4 lines for every read: the information line, the sequence line, the + line and the quality line.
     * <br> The first read is searched for the barcode + enzyme site (looked to all known samples)
     * <br> if no barcode + enzyme (perfect match) is found, a InvalidReadException is thrown
     * <br> else these are trimed from the read (both from the sequence as the quality)
     * <br> then the read is searched for a cutsite of the same enzyme, if any that piece + the rest is removed (both from the sequence as the quality)
     * @param read1 Map<FastqParts, String> | the 4 lines for the first read: the information line, the sequence line, the + line and the quality line (as Fastq from BioJava)
     * @return ProcessedRead | all information of the processed read
     * @throws InvalidReadException if the  doesn't has the right index (barcode + enzyme site)
     * @see FastqDemultiplex#findBestBarcode(java.lang.String) 
     * @see FastqDemultiplex#findRead1EnzymeLocation(java.lang.String, be.uzleuven.gc.logistics.gbsDemultiplex.model.Sample) 
     * @see ProcessedFragment
     */
    private ProcessedFragment parseFastqRead(FastqRead read1) throws InvalidReadException{
        //search for the optimal barcode
        //SampleBarcodeCombination sampleBarcodeCombination = this.findBestBarcode(read1.getSequence());
        SampleBarcodeCombination sampleBarcodeCombination = this.findGBSBarcode(read1.getSequence(), this.parameters.getStartDistance());
        if (sampleBarcodeCombination == null){
            //no optimal barcode found in the first read
            throw new InvalidReadException(read1, InvalidReadEnum.READ1);
        }
        Sample sample = sampleBarcodeCombination.getSample();
        int barcodeEnzymeLength = sampleBarcodeCombination.getLengthFoundBarcode();
        if (! this.parameters.keepCutSites()){
            //if the cutsites mustn't be kept
            barcodeEnzymeLength += sampleBarcodeCombination.getLengthFoundEnzyme();
        }
        //remove the barcode and the enzyme site
        int read1BarcodeLocation = sampleBarcodeCombination.getLocation();
        String read1modifiedSequence = read1.getSequence().substring(read1BarcodeLocation + barcodeEnzymeLength, read1.getSequence().length() - (this.longestBarcodeLength - sample.getBarcode().length()));
        String read1modifiedQuality = read1.getQuality().substring(read1BarcodeLocation + barcodeEnzymeLength, read1.getSequence().length() - (this.longestBarcodeLength - sample.getBarcode().length()));
        
        //find the next enzyme site (if there is any)
        int read1EndLocation = -1;
        int[] read1EndLocationLength = {-1, 0};
//        EnzymeComparator enzymeComparator = new EnzymeComparator();
//        if (enzymeComparator.compare(sample.getEnzyme(), this.parameters.getNeutralEnzyme()) != 0){
            read1EndLocationLength = this.findRead1EnzymeLocation(read1modifiedSequence, sample);
            read1EndLocation = read1EndLocationLength[0];
//        }
        String read1optimalSequence = read1modifiedSequence;
        String read1optimalQuality = read1modifiedQuality;
        if (read1EndLocation != -1){
            //compliment barcode found
            if (this.parameters.keepCutSites() && ! this.parameters.isRadData()){
                //if the cutsites must be kept
                read1EndLocation += read1EndLocationLength[1];
            }
            read1optimalSequence = read1modifiedSequence.substring(0, read1EndLocation);
            read1optimalQuality = read1modifiedQuality.substring(0, read1EndLocation);
        }
        //save results
        HashMap<FastqParts, String> result1 = new HashMap();
        result1.put(FastqParts.DESCRIPTION, read1.getDescription());
        result1.put(FastqParts.SEQUENCE, read1optimalSequence);
        result1.put(FastqParts.QUALITY, read1optimalQuality);
        
        FastqRead fastqRead = new FastqRead(read1.getDescription(), read1optimalSequence, read1optimalQuality);
        
        return new ProcessedFragment(sample, fastqRead, sampleBarcodeCombination.getMismatches());
    }
    
    
    
    /**
     * search for a enzyme cutsite or complement in the given sequence (cutsite found in given sample)
     * <br> returns -1 if no cutsite is found
     * <br> returns a location (int) if a cutsite is found
     * @param sequence String | the sequence where is searched in
     * @param sample Sample | the sample from which compliment cutsites is used
     * @return int | -1 if no cutsite is found, else the location of the cutsite
     */
    private int findFirstEnzymeLocation(String sequence, Sample sample){
        //init location
        int bestLocation = -1;
        //try every cutsite of the known enzyme
        for (String enzymeCutSite : sample.getEnzyme().getInitialCutSiteRemnant()){
            //try to find the location of the enzyme (-1 is returned if the barcode isn't found)
            int location = sequence.indexOf(enzymeCutSite);
            if (location != -1){
                //location found
                if (bestLocation == -1 || location < bestLocation){
                    //location is better then previous one
                    bestLocation = location;
                }
            }
        }
        return bestLocation;
    }
    
    /**
     * Goes over the sequence and search the if the piece of the sequence is equivalent to a possible enzyme cutsite (possible enzyme mismatches)
     * <br> if the found piece is equivalent and the complete check option is false, the location is returned
     * <br> if the complete check option is true:
     * <br> if the length of the basepairs after the found site is smaller then 5, the location is returned
     * <br> if the length of the basepairs is 5 or longer, the HammingsDistance is used to see if its equivalent to the start of the adaptor (use possible enzyme mismatches)
     * <br> if it is equivalent, the location is returned, else the sequence is continued
     * <br> if no cutsite is found (any more) -1 is returned
     * @param sequence String | the whole sequence to search the possible enzyme cut site
     * @param sample Sample | the sample of the sequence (contains the enzyme)
     * @return int[location, length] | the location of the cutsite, or -1 if no location is found, length of the cutsite (indels)
     * @see FindingDistanceAlgorithm#isEquivalent(java.lang.String, java.lang.String, int) 
     * @see FastqDemultiplex#getAllowedMismatchesEnzyme() 
     */
    private int[] findRead1EnzymeLocation(String sequence, Sample sample){
        //get all possible cutsites
        HashSet<String> cutsites = new HashSet<String>(sample.getEnzyme2().getComplementCutSiteRemnant());

        int[] bestIndex = {-1, 0};
        EnzymeComparator enzymeComparator = new EnzymeComparator();
        if (this.parameters.isRadData() || enzymeComparator.compare(sample.getEnzyme2(), this.parameters.getNeutralEnzyme()) == 0){
            //rad data, so check for adaptor (with out enzyme cut site)
            int adaptorMismatches = this.parameters.getAdaptorLigaseMismatches();
            if (adaptorMismatches == -1) adaptorMismatches = 0;
            int[] locationLength = {-1, 0};
            for(int mis = 0; mis <= adaptorMismatches && locationLength[0] != -1; mis++){
                locationLength = this.findingDistanceAlgorithm.indexOf(sequence, this.parameters.getCommonAdaptor().substring(0, this.parameters.getAdaptorCompareSize()), mis);
            }
            if (locationLength[0] == -1){
                int mis = adaptorMismatches;
                for (int adaptorSize = this.parameters.getAdaptorCompareSize(); adaptorSize >= 4 && locationLength[0] == -1; adaptorSize--){
                    if (adaptorSize < this.parameters.getAdaptorCompareSize() / 2){
                        mis = 0;
                    }
                    locationLength = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(sequence.length() - adaptorSize), this.parameters.getCommonAdaptor().substring(0, adaptorSize), mis);
                    if (locationLength[0] != -1){
                        locationLength[0] = locationLength[0] + (sequence.length() - adaptorSize);
                    }
                }
            }
            if (bestIndex[0] == -1 || bestIndex[0] > locationLength[0]){
                bestIndex = locationLength;
            }
        }else{
            //GBS data, so check for enzyme cut site and adaptor
            for (String enzyme : cutsites){
                String endSequence = "";
                int endSeqMismatches = 0;
                endSequence = this.parameters.getCommonAdaptor().substring(0, this.parameters.getAdaptorCompareSize());
                if (this.parameters.useDoubleBarcodes()){
                    endSequence = BasePair.getComplementSequence(sample.getBarcodeSecond()) + this.parameters.getCommonAdaptor();
                    endSequence = endSequence.substring(0, this.parameters.getAdaptorCompareSize());
                }
                endSeqMismatches = this.parameters.getAdaptorLigaseMismatches();
                if (endSeqMismatches < 0) endSeqMismatches = 0;
                int[] locationLength = this.findingDistanceAlgorithm.indexOf(sequence, enzyme, endSequence, this.parameters.getAllowedMismatchesEnzyme(), endSeqMismatches);
                
                if (locationLength[0] == -1){
                    //if the sequence is search and no enzyme + adaptor is found, look at the end of the sequence for the enzyme (no adaptor)
                    for (int i = 0; i < this.parameters.getAdaptorCompareSize(); i++){
                        int[] newLocationLength = this.findingDistanceAlgorithm.indexOf(sequence.substring(sequence.length() - this.parameters.getAdaptorCompareSize() - enzyme.length() + i), enzyme, endSequence.substring(0, this.parameters.getAdaptorCompareSize() - i), this.parameters.getAllowedMismatchesEnzyme(), 0);
                        if (newLocationLength[0] != -1){
                            newLocationLength[0] = sequence.length() - this.parameters.getAdaptorCompareSize() - enzyme.length() + i + newLocationLength[0];
                            locationLength = newLocationLength;
                            i = this.parameters.getAdaptorCompareSize();
                        }
                    }

                }
                if (locationLength[0] == -1){
                    //if the sequence is search and no enzyme + adaptor is found, look at the end of the sequence for the enzyme
                    for (int i = 0; i < enzyme.length(); i++){
                        int[] newLocationLength = this.findingDistanceAlgorithm.indexOf(sequence.substring(sequence.length() - enzyme.length() + i), enzyme.substring(0, enzyme.length() - i), 0);
                        if (newLocationLength[0] != -1){
                            newLocationLength[0] = sequence.length() - enzyme.length() + i + newLocationLength[0];
                            locationLength = newLocationLength;
                            i = enzyme.length();
                        }
                    }

                }
                if ((locationLength[0] < bestIndex[0] && locationLength[0] != -1) || bestIndex[0] == -1){
                    bestIndex = locationLength;
                    bestIndex[1] = enzyme.length();
                    if (bestIndex[0] + enzyme.length() >= sequence.length()){
                        bestIndex[0] = sequence.length();
                        bestIndex[1] = 0;
                    }
                }
            }
        }
        return bestIndex;
    }
    
    
    /**
     * if the complete check option is used:
     * <br> then the complement of the found barcode and enzyme combination is searched in the given sequence (used mismatches for barcode enzyme combination)
     * <br> if this is found the location is returned
     * <br> if the complete check option isn't used:
     * <br> the first possible enzyme cutsite is searched, if found the location is returned
     * <br> if no locations are found -1 is returned
     * @param sequence String | the sequence to search in
     * @param sample Sample | the sample with the possible enzyme site
     * @param foundEnzyme String | the found enzyme (NOT the complement) barcode is in sample
     * @return int[location, length] | the location of the first enzyme (don't do a complete check) or the location of the complement barcode enzyme combination (do a complete check) or -1 if no location is found
     */
    private int[] findRead2EnzymeLocation(String sequence, Sample sample, String foundEnzyme){
        if (this.parameters.completeCheck()){
            //uses the complete check option: look for the complement of the found enzyme => look for the reverse barcode => and look for the adaptor there after
            String complementFoundEnzyme = BasePair.getComplementSequence(foundEnzyme);
            int barcodeEnzymeLength = complementFoundEnzyme.length() + sample.getBarcode().length();
            int extraAdaptorSearch = this.longestBarcodeLength;
            if (extraAdaptorSearch > this.parameters.getCommonAdaptor().length()){
                extraAdaptorSearch = this.parameters.getCommonAdaptor().length();
            }
//            int extraAdaptorSearch = 4;
            int[] place = {-1, 0};
            boolean searchMore = true;
            int posloc = -1;
            int mismatches = 1 + this.parameters.getAllowedMismatchesBarcode(sample) + this.parameters.getAllowedMismatchesEnzyme() + this.parameters.getAdaptorLigaseMismatches();
            while (searchMore){
                //find index of complement enzyme + complement barcode
                int[] posplace = this.findingDistanceAlgorithm.indexOf(sequence.substring(posloc + 1), complementFoundEnzyme, sample.getComplementBarcode(), this.parameters.getAllowedMismatchesEnzyme(), this.parameters.getAllowedMismatchesBarcode(sample));
                
                if (posplace[0] == -1){
                    //no index found
                    searchMore = false;
                }else if (posloc + 1 + posplace[0] < sequence.length() - complementFoundEnzyme.length() - sample.getComplementBarcode().length()
//                        && posloc + 1 + posplace[0] + complementFoundEnzyme.length() + sample.getComplementBarcode().length() + extraAdaptorSearch < sequence.length()
                        ){
                    //if new found possible location is smaller then the maximum location
                    //and the the new found location + enzyme + barcode + extraAdaptor search is smaller than the length of the sequence
                    extraAdaptorSearch = this.parameters.getAdaptorCompareSize();
                    if (extraAdaptorSearch >= sequence.length() - (posloc + 1 + posplace[0] + complementFoundEnzyme.length() + sample.getComplementBarcode().length())){
                        extraAdaptorSearch = sequence.length() - (posloc + 1 + posplace[0] + complementFoundEnzyme.length() + sample.getComplementBarcode().length()) -1;
                    }
                    if (this.findingDistanceAlgorithm.isEquivalent(sequence.substring(posloc + 1 + posplace[0] + complementFoundEnzyme.length() + sample.getComplementBarcode().length(), posloc + 1 + posplace[0] + complementFoundEnzyme.length() + sample.getComplementBarcode().length() + extraAdaptorSearch), this.parameters.getCommonAdaptor().substring(0, extraAdaptorSearch), this.parameters.getAdaptorLigaseMismatches())){
                        //common adaptor found
                        int[] newMismatch = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(posloc + 1 + posplace[0], posloc + 1 + posplace[0] + complementFoundEnzyme.length() + sample.getComplementBarcode().length() + extraAdaptorSearch), complementFoundEnzyme + sample.getComplementBarcode() + this.parameters.getCommonAdaptor().substring(0, extraAdaptorSearch), (this.parameters.getAllowedMismatchesEnzyme() + this.parameters.getAllowedMismatchesBarcode(sample) + this.parameters.getAdaptorLigaseMismatches()));
                        posloc += 1 + posplace[0];
                        if (newMismatch[0] == -1){
                            //normaly can not occure
                        }else{
                            if (newMismatch[0] < mismatches){
                                //less mismatches found then previous try => correcter
                                place[0] = posloc;
                                mismatches = newMismatch[0];
                            }else{
                                //more mismatches found => probebly less correct
                            }
                        }
                        posloc += 1;
                        if (mismatches == 0){
                            searchMore = false;
                        }
                    }else{
                        posloc += 1 + posplace[0];
                    }
                }else if (posplace[0] != -1){
                    searchMore = false;
                }
                
            }
            place[1] = complementFoundEnzyme.length();
            if (place[0] == -1){
                //find with no barcode
//                int[] newplace = this.findingDistanceAlgorithm.indexOf(sequence.substring(sequence.length() - complementFoundEnzyme.length() - sample.getBarcode().length() - extraAdaptorSearch), complementFoundEnzyme, this.parameters.getAllowedMismatchesEnzyme());
                int[] newplace = this.findingDistanceAlgorithm.indexOf(sequence.substring(sequence.length() - complementFoundEnzyme.length() - sample.getBarcode().length()), complementFoundEnzyme, this.parameters.getAllowedMismatchesEnzyme());
                
                if (newplace[0] != -1){
//                    place[0] = (sequence.length() - complementFoundEnzyme.length() - sample.getBarcode().length() - extraAdaptorSearch + newplace[0]);
//                    if(this.findingDistanceAlgorithm.isEquivalent(sequence.substring(sequence.length() - sample.getBarcode().length() + newplace[0]), sample.getComplementBarcode().substring(0, sample.getBarcode().length() + complementFoundEnzyme.length() - newplace[0]), this.parameters.getAllowedMismatchesBarcode(sample))){
//                        place[0] = (sequence.length() - complementFoundEnzyme.length() - sample.getBarcode().length() + newplace[0]);
//                    }else{
//                        newplace[0] = -1;
//                    }
                    place[0] = (sequence.length() - complementFoundEnzyme.length() - sample.getBarcode().length() + newplace[0]);
                }
                else{
//                if (newplace[0] == -1){
                    //find with part enzyme
                    for (int ei=0; ei < complementFoundEnzyme.length(); ei++){
                        newplace = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(sequence.length() - complementFoundEnzyme.length() + ei), complementFoundEnzyme.substring(0, complementFoundEnzyme.length() - ei), 0);
                        if (newplace[0] != -1){
                            place[0] = (sequence.length() - complementFoundEnzyme.length() + ei);
                            place[1] = complementFoundEnzyme.length() - ei;
                            ei = complementFoundEnzyme.length();
                        }
                    }
                }
            }
            return place;
        }else{
            EnzymeComparator enzymeComparator = new EnzymeComparator();
            if (this.parameters.isRadData() || enzymeComparator.compare(sample.getEnzyme(), this.parameters.getNeutralEnzyme()) == 0){
                return new int[]{-1,0};
            }
            //uses complete digest option: look only for an enzyme cutsite
            //get all possible cutsites
            HashSet<String> cutsites = new HashSet<String>(sample.getEnzyme().getComplementCutSiteRemnant());
            //check every enzyme
            int [] locationLength = new int[] {-1, 0};
            for (String enzyme : cutsites){
                int[] newlocationLength;
                if ((newlocationLength = this.findingDistanceAlgorithm.indexOf(sequence, enzyme, this.parameters.getAllowedMismatchesEnzyme()))[0] != -1){
                    if (locationLength[0] == -1 || (newlocationLength[0] < locationLength[0] && newlocationLength[0] != -1)){
                        locationLength = newlocationLength;
                    }
                }
            }
            return locationLength;
        }
    }
    
    
    /**
     * search for the best barcode + enzyme combination in the given sequence
     * <br> uses all known samples, and tries to find the barcode-enzyme combination in front of the sequence
     * <br> returns null if no combination is found
     * <br> returns a SampleBarcodeCombination of the best found sample and the barcode-enzyme combination
     * @param sequence String | the sequence to search
     * @return null if no best combination is found, else SampleBarcodeCombination of the best combination
     */
    private SampleBarcodeCombination findBestBarcode(String sequence){
        //sample, barcode-enzyme and location init
        Sample bestSample = null;
        String bestEnzymeCutSite = "";
        int bestLocation = -1;
        int barcodelength = 0;
        int enzymelength = 0;
        //try every sample
        for (Sample sample : this.sampleList){
            //try every barcode
            String barcode = sample.getBarcode();
            for (String enzymecutsite : sample.getEnzyme().getInitialCutSiteRemnant()){
                //try to find the location of the barcode (-1 is returned if the barcode isn't found)
                int location = sequence.indexOf(barcode + enzymecutsite);
                if (location != -1){
                    //barcode found
                    if (bestLocation == -1 || location < bestLocation){
                        //location is beter then previous one
                        bestSample = sample;
                        bestEnzymeCutSite = enzymecutsite;
                        bestLocation = location;
                        barcodelength = barcode.length();
                        enzymelength = enzymecutsite.length();
                    }
                }
            }
        }
        if (bestSample == null){
            //no combination found
            return null;
        }
        //return the combination
        return new SampleBarcodeCombination(bestSample, bestEnzymeCutSite, bestLocation, 0, barcodelength, enzymelength);
        
    }
    
    /**
     * looks at the sequences first basepairs to determenate if there is a barcode + enzyme combination possible
     * <br> the hammings distance is used to allow mismatches.
     * @param sequence String | the sequence to look in
     * @param startDistance int | the max distance from the first basepair of the sequence to the first basepair of the barcode (0 is the start only), startDistance must be smaller then the MAXIMUM_DISTANCE_BETWEEN_START_AND_BARCODE
     * @return SampleBarcodeCombination | a class that is the combination of the found sample, and the found barcodeEnzyme combination + the location (currently 0)
     * @see FindingDistanceAlgorithm#isEquivalent(java.lang.String, java.lang.String, int) 
     */
    private SampleBarcodeCombination findGBSBarcode(String sequence, int startDistance){
        //for every distance
        for (int distance = 0; distance <= startDistance && distance <= this.MAXIMUM_DISTANCE_BETWEEN_START_AND_BARCODE; distance++){
            //try every sample
            for (Sample sample : this.sampleList){
                //try every barcode
                int[] barcodeLocationLength = new int[]{-1};
                if (this.parameters.useSelfCorrectingBarcodes()){
                    BarcodeCorrectingAlgorithm barcodeCorrectingAlgorithm = new BarcodeCorrectingAlgorithm();
                    String newBarcode = barcodeCorrectingAlgorithm.correctBarcode(sequence.substring(distance, distance + sample.getBarcode().length()));
                    if (newBarcode.equals(sample.getBarcode())){
                        //distance, length
                        barcodeLocationLength = new int[] {0, sample.getBarcode().length()};
                    }
                }else{
                    barcodeLocationLength = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(distance), sample.getBarcode(), this.parameters.getAllowedMismatchesBarcode(sample));
                }
                if (barcodeLocationLength[0] != -1){
                    //try every enzyme site
                    for (String enzymeCutSite : sample.getEnzyme().getInitialCutSiteRemnant()){
                        int[] cutsiteLocationLength;
                        if ((cutsiteLocationLength = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(distance + barcodeLocationLength[1]), enzymeCutSite, this.parameters.getAllowedMismatchesEnzyme()))[0] != -1){
                            //check on adaptor ligase
                            if (this.parameters.getAdaptorLigaseMismatches() != -1){
                                String adaptor = this.parameters.getCommonAdaptor();
                                if (this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(distance + barcodeLocationLength[1] + cutsiteLocationLength[1]), adaptor, this.parameters.getAdaptorLigaseMismatches())[0] != -1){
                                    return null;
                                }
                            }
                            return new SampleBarcodeCombination(sample, enzymeCutSite, distance, barcodeLocationLength[0], barcodeLocationLength[1], cutsiteLocationLength[1]);
                        }
                    }
                }
            }
        }
        return null;
    }
    
    /**
     * looks at the sequences first basepairs to determenate if there is a barcode + enzyme combination possible
     * <br> the hammings distance is used to allow mismatches.
     * @param sequence String | the sequence to look in
     * @param startDistance int | the max distance from the first basepair of the sequence to the first basepair of the barcode (0 is the start only), startDistance must be smaller then the MAXIMUM_DISTANCE_BETWEEN_START_AND_BARCODE
     * @return SampleBarcodeCombination | a class that is the combination of the found sample, and the found barcodeEnzyme combination + the location (currently 0)
     * @see FindingDistanceAlgorithm#isEquivalent(java.lang.String, java.lang.String, int) 
     */
    private SampleBarcodeCombination[] findGBSBarcode2(String sequence, int startDistance, String sequence2){
        //for every distance
        for (int distance = 0; distance <= startDistance && distance <= this.MAXIMUM_DISTANCE_BETWEEN_START_AND_BARCODE; distance++){
            HashSet<SampleBarcodeCombination> foundSampleSet = new HashSet<SampleBarcodeCombination>();
            //try every sample
            for (Sample sample : this.sampleList){
                //try every barcode
                int[] barcodeLocationLength = new int[]{-1};
                if (this.parameters.useSelfCorrectingBarcodes()){
                    BarcodeCorrectingAlgorithm barcodeCorrectingAlgorithm = new BarcodeCorrectingAlgorithm();
                    String newBarcode = barcodeCorrectingAlgorithm.correctBarcode(sequence.substring(distance, distance + sample.getBarcode().length()));
                    if (newBarcode.equals(sample.getBarcode())){
                        //distance, length
                        barcodeLocationLength = new int[] {0, sample.getBarcode().length()};
                    }
                }else{
                    barcodeLocationLength = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(distance), sample.getBarcode(), this.parameters.getAllowedMismatchesBarcode(sample));
                }
                if (barcodeLocationLength[0] != -1){
                    //try every enzyme site
                    boolean cutsiteFound = false;
                    String exactEnzymeCutSite = "";
                    int[] cutsiteLocationLength = new int[1];
                    cutsiteLocationLength[0] = 0;
                    for (String enzymeCutSite : sample.getEnzyme().getInitialCutSiteRemnant()){
                        if ((cutsiteLocationLength = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(distance + barcodeLocationLength[1]), enzymeCutSite, this.parameters.getAllowedMismatchesEnzyme()))[0] != -1){
                            //check on adaptor ligase
                            if (this.parameters.getAdaptorLigaseMismatches() != -1){
                                String adaptor = this.parameters.getCommonAdaptor();
                                if (this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence.substring(distance + barcodeLocationLength[1] + cutsiteLocationLength[1]), adaptor, this.parameters.getAdaptorLigaseMismatches())[0] != -1){
                                    return null;
                                    //adaptor ligase
                                }
                            }
                            exactEnzymeCutSite = enzymeCutSite;
                            cutsiteFound = true;
                        }
                    }
                    if (cutsiteFound){
                        foundSampleSet.add(new SampleBarcodeCombination(sample, exactEnzymeCutSite, distance, barcodeLocationLength[0], barcodeLocationLength[1], cutsiteLocationLength[1]));
                        
                    }
                }
            }
            //found possible sampleBarcode combinations
            if (foundSampleSet.isEmpty()){
                return null;
            }else if(! this.parameters.useDoubleBarcodes()){
                //only 1 barcode
                if (foundSampleSet.size() == 1){
                    //found exact 1 possibility == report
                    SampleBarcodeCombination[] combinationArray = new SampleBarcodeCombination[1];
                    combinationArray[0] = foundSampleSet.iterator().next();
                    return combinationArray;
                }else{
                    //found multiple possibilities == error
                    return null;
                }
            }else {
                //use 2 barcodes
                HashMap<Sample, SampleBarcodeCombination> testSamplesMap = new HashMap<Sample, SampleBarcodeCombination>();
                for (SampleBarcodeCombination sampleBarcodeCombination : foundSampleSet){
                    testSamplesMap.put(sampleBarcodeCombination.getSample(), sampleBarcodeCombination);
                }
                HashMap<Sample, SampleBarcodeCombination> resultMap = new HashMap<Sample, SampleBarcodeCombination>();
                for (Sample sample : testSamplesMap.keySet()){
                //try every barcode
                    int[] barcodeLocationLength = new int[]{-1};
                    if (this.parameters.useSelfCorrectingBarcodes()){
                        BarcodeCorrectingAlgorithm barcodeCorrectingAlgorithm = new BarcodeCorrectingAlgorithm();
                        String newBarcode = barcodeCorrectingAlgorithm.correctBarcode(sequence2.substring(distance, distance + sample.getBarcodeSecond().length()));
                        if (newBarcode.equals(sample.getBarcodeSecond())){
                            //distance, length
                            barcodeLocationLength = new int[] {0, sample.getBarcodeSecond().length()};
                        }
                    }else{
                        barcodeLocationLength = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence2.substring(distance), sample.getBarcodeSecond(), this.parameters.getAllowedMismatchesBarcode(sample));
                    }
                    if (barcodeLocationLength[0] != -1){
                        //try every enzyme site
                        boolean cutsiteFound = false;
                        String exactEnzymeCutSite = "";
                        int[] cutsiteLocationLength = new int[1];
                        cutsiteLocationLength[0] = 0;
                        for (String enzymeCutSite : sample.getEnzyme2().getInitialCutSiteRemnant()){
                            if ((cutsiteLocationLength = this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence2.substring(distance + barcodeLocationLength[1]), enzymeCutSite, this.parameters.getAllowedMismatchesEnzyme()))[0] != -1){
                                //check on adaptor ligase
                                if (this.parameters.getAdaptorLigaseMismatches() != -1){
                                    String adaptor = this.parameters.getCommonAdaptor();
                                    if (this.findingDistanceAlgorithm.calculateEquivalentDistance(sequence2.substring(distance + barcodeLocationLength[1] + cutsiteLocationLength[1]), adaptor, this.parameters.getAdaptorLigaseMismatches())[0] != -1){
                                        return null;
                                    }
                                }
                                exactEnzymeCutSite = enzymeCutSite;
                                cutsiteFound = true;
                            }
                        }
                        if(cutsiteFound){
                                resultMap.put(sample, new SampleBarcodeCombination(sample, exactEnzymeCutSite, distance, barcodeLocationLength[0], barcodeLocationLength[1], cutsiteLocationLength[1]));
                        }
                    }
                }
                //found result?
                if (resultMap.isEmpty() || resultMap.size() > 1){
                    //no result, or multiple results
                    return null;
                }else{
                    //found exact 1 possibility == report
                    SampleBarcodeCombination[] combinationArray = new SampleBarcodeCombination[2];
                    for (Sample sample : resultMap.keySet()){
                        combinationArray[0] = resultMap.get(sample);
                        combinationArray[1] = resultMap.get(sample);
                    }
                    return combinationArray;
                }
            }
        }
        return null;
    }
    
    /**
     * a combination of a sample, a barcode and enzyme, the mismatches (between sequence and barcode/enzyme) and the start location of the barcode in the sequence
     */
    private class SampleBarcodeCombination{
        
        private Sample sample;
        private String enzymeCutSite;
        private int location;
        private int mismatches;
        private int lengthFoundBarcode;
        private int lengthFoundEnzyme;
        /**
         * 
         * @param sample Sample | the sample of this combination
         * @param enzymeCutSite String | the used enzyme cutsite
         * @param location int | the start location of the barcodeEnzyme
         * @param mismatches int | the amount of mismatches
         * @param length int | the length of the combination of the sample and barcode
         */
        public SampleBarcodeCombination(Sample sample, String enzymeCutSite, int location, int mismatches, int lengthFoundBarcode, int lengthFoundEnzyme){
            this.sample = sample;
            this.enzymeCutSite = enzymeCutSite;
            this.location = location;
            this.mismatches = mismatches;
            this.lengthFoundBarcode = lengthFoundBarcode;
            this.lengthFoundEnzyme = lengthFoundEnzyme;
        }
        
        /**
         * 
         * @return Sample | the sample of this combination
         */
        public Sample getSample(){
            return this.sample;
        }
        
        /**
         * 
         * @return String | the enzymeCutsite
         */
        public String getEnzymeCutsite(){
            return this.enzymeCutSite;
        }
        
        /**
         * 
         * @return int | the location of the barcode + enzyme in a sequence
         */
        public int getLocation(){
            return this.location;
        }
        
        /**
         * 
         * @return int | the amount of mismatches between the barcode + enzyme and the sequence
         */
        public int getMismatches(){
            return this.mismatches;
        }
        
        /**
         * 
         * @return int | the length of the found enzyme
         */
        public int getLengthFoundEnzyme(){
            return this.lengthFoundEnzyme;
        }
        
        /**
         * 
         * @return int | the lenght of the found barcode
         */
        public int getLengthFoundBarcode(){
            return this.lengthFoundBarcode;
        }
        
    }
    
    /**
     * all information of the processed fragment:
     * <br> the sample of the read, 
     * <br> read1 (as HashMap of FastqParts and String)
     * <br> read2 (only pair-end) (as HashMap of FastqParts and String)
     * <br> mismatch occured in finding the barcode/enzyme
     */
    private class ProcessedFragment{
        
        private Sample sample;
        private FastqRead read1;
        private FastqRead read2;
        private int mismatch;
        private String sequenceComment = "";
        
        /**
         * create a new ProcessedFragment (pair-end)
         * @param sample Sample | the sample of the fragment
         * @param read1 FastqRead | the first read of the fragment
         * @param read2 FastqRead | the second read of the fragment (only pair-end)
         * @param mismatch int | number of mismatches in the barcode/enzyme
         * @param sequenceComment String | the comment on the cut of the sequence
         */
        public ProcessedFragment(Sample sample, FastqRead read1, FastqRead read2, int mismatch, String sequenceComment){
            this.sample = sample;
            this.read1 = read1;
            this.read2 = read2;
            this.mismatch = mismatch;
            this.sequenceComment = sequenceComment;
        }
        
        /**
         * create a new ProcessedFragment (pair-end)
         * @param sample Sample | the sample of the fragment
         * @param read1 FastqRead | the first read of the fragment
         * @param read2 FastqRead | the second read of the fragment (only pair-end)
         * @param mismatch int | number of mismatches in the barcode/enzyme
         */
        public ProcessedFragment(Sample sample, FastqRead read1, FastqRead read2, int mismatch){
            this.sample = sample;
            this.read1 = read1;
            this.read2 = read2;
            this.mismatch = mismatch;
        }
        
        /**
         * create a new ProcessedFragment (single read)
         * @param sample Sample | the sample of the fragment
         * @param read1 FastqRead | the only read of the fragment
         * @param mismatch int | number of mismatches in the barcode/enzyme
         */
        public ProcessedFragment(Sample sample, FastqRead read1, int mismatch){
            this.sample = sample;
            this.read1 = read1;
            this.read2 = null;
            this.mismatch = mismatch;
        }
        
        /**
         * 
         * @return Sample | the sample of this fragment
         */
        public Sample getSample(){
            return this.sample;
        }
        
        /**
         * 
         * @return HashMap of FastqParts and String | the first read
         */
        public FastqRead getRead1(){
            return this.read1;
        }
        
        /**
         * 
         * @return HashMap of FastqParts and String | the second read (only by pair-end)
         */
        public FastqRead getRead2(){
            return this.read2;
        }
        
        /**
         * 
         * @return int | number of mismatches in barcode/enzyme
         */
        public int getMismatch(){
            return this.mismatch;
        }
        
        /**
         * 
         * @return String | the comment on the cut of the reads
         */
        public String getComment(){
            return this.sequenceComment;
        }
        
    }
    
    
}
