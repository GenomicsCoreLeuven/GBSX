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

import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.ErrorInLogException;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.infrastructure.fileInteractors.LoggerFile;
import be.uzleuven.gc.logistics.GBSX.utils.FileLocker;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqScores;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class DemultiplexStats {
    
    /**
     * stats per sample, then per mismatch the count
     */
    private HashMap<Sample, HashMap<Integer, Integer>> stats;
    /**
     * used log file
     */
    private LoggerFile logFile;
    /**
     * the count of the undetermined sequences
     */
    private int undetermined;
    /**
     * count of all bases in the sequence (quality)
     */
    private HashMap<Sample, Long> basecall_count;
    /**
     * total count of the quality of all sequences of that sample
     */
    private HashMap<Sample, Long> basecall_qual;
    /**
     * the count of all basepairs with a score above 30
     */
    private HashMap<Sample, Long> basecall_above_30;
    /**
     * the count of all rejected reads of that sample
     */
    private HashMap<Sample, Integer> rejected_reads;
    /**
     * the FastqScore used for the quality
     */
    private FastqScores fastqScore;
    
    private FileLocker fileLocker = new FileLocker();
    
    /**
     * creates a new DemultiplexStats
     * @param sampleList List of Sample | a list of all samples that the stats will use
     * @param logFile LoggerFile | the file to log errors
     */
    public DemultiplexStats(List<Sample> sampleList, LoggerFile logFile, FastqScores fastqScore){
        this.logFile = logFile;
        this.undetermined = 0;
        this.stats = new HashMap();
        this.basecall_count = new HashMap();
        this.basecall_qual = new HashMap();
        this.basecall_above_30 = new HashMap();
        this.rejected_reads = new HashMap();
        for (Sample sample : sampleList){
            HashMap<Integer, Integer> mismatchMap = new HashMap();
            this.stats.put(sample, mismatchMap);
            this.basecall_count.put(sample, 0L);
            this.basecall_qual.put(sample, 0L);
            this.basecall_above_30.put(sample, 0L);
            this.rejected_reads.put(sample, 0);
        }
        if (fastqScore == null){
            this.fastqScore = FastqScores.getStandard();
        }
        this.fastqScore = fastqScore;
    }
    
    /**
     * Adds a stat to the file (adds the mismatch, basequality, and other statistics)
     * @param sample Sample | the sample of this read
     * @param numberOfMismatches int | the number of mismatches occured in the barcode
     * @param quality String | the quality string of the sample (for paired end: concatination of both quality strings)
     */
    public void addStat(Sample sample, int numberOfMismatches, String quality){
        if (this.fileLocker.lock()){
            HashMap<Integer, Integer> misMap = this.stats.get(sample);
            if (misMap == null){
                //create everything for the new sample
                HashMap<Integer, Integer> mismatchMap = new HashMap();
                this.stats.put(sample, mismatchMap);
                this.basecall_count.put(sample, 0L);
                this.basecall_qual.put(sample, 0L);
                this.basecall_above_30.put(sample, 0L);
                misMap = this.stats.get(sample);
            }
            //get the map of mismatches and add the mismatch
            Integer mismatch = misMap.get(numberOfMismatches);
            if (mismatch == null){
                misMap.put(numberOfMismatches, 1);
            }else{
                misMap.put(numberOfMismatches, ++mismatch);
            }
            this.stats.put(sample, misMap);
            //calculate the phred score, the basecount, base quality, ...
            for (char ascii_score : quality.toCharArray()){
                int phred_score = (int) ascii_score - this.fastqScore.getStartScore() + this.fastqScore.getMinScore();
                this.basecall_count.put(sample, this.basecall_count.get(sample) + 1);
                this.basecall_qual.put(sample, this.basecall_qual.get(sample) + phred_score);
                if (phred_score >= 30){
                    this.basecall_above_30.put(sample, this.basecall_above_30.get(sample) + 1);
                }
            }
            this.fileLocker.unlock();
        }
    }
    
    /**
     * add 1 undetermined read
     */
    public void addUndeterminedStat(){
        if (this.fileLocker.lock()){
            this.undetermined++;
            this.fileLocker.unlock();
        }
    }
    
    /**
     * Adds a rejected Read
     * @param sample Sample | the sample of this read
     */
    public void addRejectedRead(Sample sample){
        if (this.fileLocker.lock()){
            this.rejected_reads.put(sample, this.rejected_reads.get(sample) + 1);
            this.fileLocker.unlock();
        }
    }
    
    /**
     * saves the generated stats to gbsDemultiplex.stats
     * <br> a new file found in the given outputFile (directory)
     * @param maxMismatch int | the number of the max mismatches
     * @param outputFile String | the directory of the file
     */
    public void saveStats(int maxMismatch, String outputFile){
        this.fileLocker.waitTillCompleteUnlock();
        boolean rejected = false;
        for (int numberOfReject : this.rejected_reads.values()){
            if (numberOfReject > 0){
                rejected = true;
            }
        }
        //file header
        String statsString = "sampleID" + "\t" + "barcode" + "\t" + "enzyme";
        statsString += "\t" + "total.count" + "\t" + "total.perc";
        for (int index=0; index <= maxMismatch; index++){
            statsString += "\t" + "mismatch." + index + ".count" + "\t" + "mismatch." + index + ".perc";
        }
        statsString += "\t" + "basecall.count" + "\t" + "basecall.above.30.perc" + "\t" + "basecall.qual.avg";
        if (rejected){
            statsString += "\t" + "rejected.count";
        }
        //stats par sample
        int reads_total = this.undetermined;
        for (Sample sample : this.stats.keySet()){
            HashMap<Integer, Integer> misMap = this.stats.get(sample);
            for (int index=0; index <= maxMismatch; index++){
                if (misMap.containsKey(index)){
                    reads_total += misMap.get(index);
                }
            }
        }
        ArrayList<Sample> sampleList = new ArrayList(this.stats.keySet());
        Collections.sort(sampleList);
        for (Sample sample : sampleList){
            //init
            HashMap<Integer, Integer> misMap = this.stats.get(sample);
            int[] misMatches = new int[maxMismatch + 1];
            int sampleTotal = 0;
            //go over each possible mismatch
            for (int index=0; index <= maxMismatch; index++){
                if (misMap.containsKey(index)){
                    misMatches[index] = misMap.get(index);
                    sampleTotal += misMap.get(index);
                }else{
                    misMatches[index] = 0;
                }
            }
            //write to statsString, and calculates possible stats
            statsString += "\n";
            statsString += sample.getSampleID() + "\t" + sample.getBarcode() + "\t" + sample.getEnzyme().getName();
            statsString += "\t" + sampleTotal + "\t" + (sampleTotal / (double) reads_total);
            for (int index=0; index <= maxMismatch; index++){
                double perc = 0.0;
                if (sampleTotal != 0){
                    perc = (misMatches[index] * 1.0) / (sampleTotal * 1.0);
                }
                statsString += "\t" + misMatches[index] + "\t" + perc;
            }
            statsString += "\t" + this.basecall_count.get(sample) + "\t" + (this.basecall_above_30.get(sample) / (double) this.basecall_count.get(sample)) + "\t" + (this.basecall_qual.get(sample) / (double) this.basecall_count.get(sample));
            if (rejected){
                statsString += "\t" + this.rejected_reads.get(sample);
            }
            
        }
        statsString += "\n" + "undetermined" + "\t" + "\t";
        statsString += "\t" + this.undetermined + "\t" + (this.undetermined / (double) reads_total);
        
        //write the string to a file
        File file = new File(outputFile + System.getProperty("file.separator") + "gbsDemultiplex.stats");
        try {
            file.createNewFile();
            BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file));
            bufferedWriter.write(statsString);
            bufferedWriter.close();
        } catch (IOException ex) {
            Logger.getLogger(DemultiplexStats.class.getName()).log(Level.SEVERE, null, ex);
            try {
                this.logFile.addToLog("Couldn't write the stats");
            } catch (ErrorInLogException ex1) {
                Logger.getLogger(DemultiplexStats.class.getName()).log(Level.SEVERE, null, ex1);
            }
        }
        
    }
    
}
