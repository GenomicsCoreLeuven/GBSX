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
import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.ProcessedFragment;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.SampleBarcodeCombination;
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
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
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
     * @throws be.uzleuven.gc.logistics.GBSX.utils.exceptions.StopExcecutionException 
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
            System.err.println("Mix of single and double barcodes found.");
            throw new RuntimeException("Mix of single and double barcodes found.");
        }
        //double barcodes only possible with paired end sequencing
        if (doubleBarcodeCount == this.sampleList.size()){
            this.parameters.setDoubleBarcodes();
            if(this.parameters.getFastqFile2() == null){
                //single end demultiplex
                System.err.println("Use Double Barcodes, but no paired end fastq file found.");
            throw new RuntimeException("Use Double Barcodes, but no paired end fastq file found.");
            }
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
            
            ProgressTracker progressTracker = new ProgressTracker();
            
            ArrayList<DemultiplexThread> threadlist = new ArrayList<DemultiplexThread>();
            
            for (int i=0; i < this.parameters.getThreadNumber(); i++){
                DemultiplexThread demultiplexThread = new DemultiplexThread(this.parameters, statsFile, loggerFile, progressTracker, fastq1Reader, 
                    undeterminedFastqFile, sampleFiles, longestBarcodeLength, findingDistanceAlgorithm, correctionLog, sampleList);
                threadlist.add(demultiplexThread);
            }
                        
            ExecutorService service = Executors.newFixedThreadPool(this.parameters.getThreadNumber());
            List<Future<Runnable>> futures = new ArrayList<Future<Runnable>>();
            
            for (DemultiplexThread demultiplexThread : threadlist){
                Future f = service.submit(demultiplexThread);
                futures.add(f);
            }
            
            
            for (Future<Runnable> f : futures){
                try {
                    f.get();
                } catch (InterruptedException ex) {
                    Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
                } catch (ExecutionException ex) {
                    Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
            service.shutdownNow();
            
            progressTracker.showProgress();
                        
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
        System.out.println("USE DUAL BARCODING: " + this.parameters.useDoubleBarcodes());
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
            
            
            
            ProgressTracker progressTracker = new ProgressTracker();
            
            ArrayList<DemultiplexThread> threadlist = new ArrayList<DemultiplexThread>();
            
            for (int i=0; i < this.parameters.getThreadNumber(); i++){
                DemultiplexThread demultiplexThread = new DemultiplexThread(this.parameters, statsFile, loggerFile, progressTracker, fastqReader, 
                    undeterminedFastqFile, sampleFiles, longestBarcodeLength, findingDistanceAlgorithm, correctionLog, sampleList);
                threadlist.add(demultiplexThread);
            }
                        
            ExecutorService service = Executors.newFixedThreadPool(this.parameters.getThreadNumber());
            List<Future<Runnable>> futures = new ArrayList<Future<Runnable>>();
            
            for (DemultiplexThread demultiplexThread : threadlist){
                Future f = service.submit(demultiplexThread);
                futures.add(f);
            }
            
            
            for (Future<Runnable> f : futures){
                try {
                    f.get();
                } catch (InterruptedException ex) {
                    Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
                } catch (ExecutionException ex) {
                    Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
            service.shutdownNow();
            
            progressTracker.showProgress();
            
            
            
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
    
    
}
