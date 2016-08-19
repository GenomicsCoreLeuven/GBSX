/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package be.uzleuven.gc.logistics.GBSX.demultiplexer.application;

import be.uzleuven.gc.logistics.GBSX.demultiplexer.application.distanceAlgorithms.FindingDistanceAlgorithm;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.ErrorInLogException;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.InvalidReadEnum;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.InvalidReadException;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.infrastructure.fileInteractors.LoggerFile;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.DemultiplexParameters;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.ProcessedFragment;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqBufferedReader;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqBufferedWriter;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqPairBufferedReader;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqPairBufferedWriter;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqRead;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author koen
 */
public class DemultiplexThread implements Runnable{

    private Thread t;
    private final ProgressTracker progressTracker;
    private final FastqReadParser fastqReadParser;
    private final DemultiplexParameters parameters;
    private final LoggerFile loggerFile;
    private final DemultiplexStats statsFile;
    private FastqBufferedReader srFastqReader;
    private FastqBufferedWriter srUndeterminedFastqFile;
    private FastqPairBufferedReader peFastqReader;
    private FastqPairBufferedWriter peUndeterminedFastqFile;
    private HashMap<Sample, FastqBufferedWriter> srSampleFiles;
    private HashMap<Sample, FastqPairBufferedWriter> peSampleFiles;
    
    public DemultiplexThread(DemultiplexParameters parameters,
            DemultiplexStats demultiplexStats, LoggerFile loggerFile, ProgressTracker progressTracker,
            FastqBufferedReader fastq1Reader, FastqBufferedWriter undeterminedFastqFile, 
            HashMap<Sample, FastqBufferedWriter> sampleFiles, int longestBarcode,
            FindingDistanceAlgorithm findingDistanceAlgorithm, CorrectionLog correctionLog,
            ArrayList<Sample> sampleList){
        this.parameters = parameters;
        this.statsFile = demultiplexStats;
        this.loggerFile = loggerFile;
        this.progressTracker = progressTracker;
        this.srUndeterminedFastqFile = undeterminedFastqFile;
        this.srFastqReader = fastq1Reader;
        this.srSampleFiles = sampleFiles;
        this.fastqReadParser = new FastqReadParser(this.parameters, longestBarcode, findingDistanceAlgorithm, correctionLog, sampleList);
        this.peFastqReader = null;
        this.peSampleFiles = null;
        this.peUndeterminedFastqFile = null;
    }
    
    public DemultiplexThread(DemultiplexParameters parameters,
            DemultiplexStats demultiplexStats, LoggerFile loggerFile, ProgressTracker progressTracker,
            FastqPairBufferedReader fastq1Reader, FastqPairBufferedWriter undeterminedFastqFile, 
            HashMap<Sample, FastqPairBufferedWriter> sampleFiles, int longestBarcode,
            FindingDistanceAlgorithm findingDistanceAlgorithm, CorrectionLog correctionLog,
            ArrayList<Sample> sampleList){
        this.parameters = parameters;
        this.statsFile = demultiplexStats;
        this.loggerFile = loggerFile;
        this.progressTracker = progressTracker;
        this.peUndeterminedFastqFile = undeterminedFastqFile;
        this.peFastqReader = fastq1Reader;
        this.peSampleFiles = sampleFiles;
        this.fastqReadParser = new FastqReadParser(this.parameters, longestBarcode, findingDistanceAlgorithm, correctionLog, sampleList);
        this.srFastqReader = null;
        this.srSampleFiles = null;
        this.srUndeterminedFastqFile = null;
    }
    
    public void run() {
        try{
            if (this.peFastqReader == null){
                //Single end
                this.singleDemultiplex();
            }else{
                //Pair end
                this.pairedDemultiplex();
            }
        }catch(IOException e){
            
        }
    }
    
    private void pairedDemultiplex() throws IOException{
        FastqRead fastq1 = null;
        FastqRead fastq2 = null;
        HashMap<String, FastqRead> readMap;
        while (((readMap = peFastqReader.next()) != null) && ((fastq1 = readMap.get("r1")) != null) && ((fastq2 = readMap.get("r2")) != null)){
            progressTracker.addProgress();
            //read the next fastq line
            try {
                //try to parse the fastq
                ProcessedFragment newReads = this.fastqReadParser.parseFastqRead(fastq1, fastq2);
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
                    peUndeterminedFastqFile.write(fastq1, fastq2);
                    statsFile.addUndeterminedStat();
                }else{
                    //copies the read to the corresponding sample
                    peSampleFiles.get(newReads.getSample()).write(newReads.getRead1(), newReads.getRead2());
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
                    peUndeterminedFastqFile.write(fastq1, fastq2);
                    statsFile.addUndeterminedStat();
                } catch (IOException ioex) {
                    this.writeToLog("ERROR in writing the fastq files for " + fastq1.getDescription());
                    Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ioex);
                }
            }
        }
    }
    
    private void singleDemultiplex() throws IOException{
        //init vars needed in the loop to go furter
        FastqRead fastq1;
        while ((fastq1 = srFastqReader.next()) != null) {
            progressTracker.addProgress();
            //read the next fastq line
            try {
                //try to parse the fastq
                ProcessedFragment newReads = this.fastqReadParser.parseFastqRead(fastq1);
                //check if a sequence is empty: correct it
                if (newReads.getRead1().getSequence().equals("")){
                    newReads = new ProcessedFragment(newReads.getSample(), new FastqRead(newReads.getRead1().getDescription(), "N", "#"), newReads.getMismatch());
                }
                //check if a sequence must be rejected (to short)
                if ((newReads.getRead1().getSequence().length() < this.parameters.getMinimumSequenceLength())
                        || (! this.parameters.keepSequencesWithN() && newReads.getRead1().getSequence().contains("N"))){
                    statsFile.addRejectedRead(newReads.getSample());
                    srUndeterminedFastqFile.write(fastq1);
                    statsFile.addUndeterminedStat();
                }else{
                    //copies the read to the corresponding sample
                    srSampleFiles.get(newReads.getSample()).write(newReads.getRead1());
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
                    srUndeterminedFastqFile.write(fastq1);
                    statsFile.addUndeterminedStat();
                } catch (IOException ioex) {
                    this.writeToLog("ERROR in writing the fastq files for " + fastq1.getDescription());
                    Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ioex);
                }
            }  
        }
    }
    
    public void start(){
        if (this.t == null){
            this.t = new Thread(this);
            this.t.start();
            System.out.println("Thread started");
        }
    }
    
    /**
     * 
     * @param log String | the log to be written to the logfile
     * @see LoggerFile#addToLog(java.lang.String) 
     */
    private void writeToLog(String log){
        try {
            this.loggerFile.addToLog(log + "\n");
        } catch (ErrorInLogException ex) {
            Logger.getLogger(FastqDemultiplex.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}


/*
class RunnableDemo implements Runnable {
   private Thread t;
   private String threadName;
   
   RunnableDemo( String name){
       threadName = name;
       System.out.println("Creating " +  threadName );
   }
   public void run() {
      System.out.println("Running " +  threadName );
      try {
         for(int i = 4; i > 0; i--) {
            System.out.println("Thread: " + threadName + ", " + i);
            // Let the thread sleep for a while.
            Thread.sleep(50);
         }
     } catch (InterruptedException e) {
         System.out.println("Thread " +  threadName + " interrupted.");
     }
     System.out.println("Thread " +  threadName + " exiting.");
   }
   
   public void start ()
   {
      System.out.println("Starting " +  threadName );
      if (t == null)
      {
         t = new Thread (this, threadName);
         t.start ();
      }
   }

}*/