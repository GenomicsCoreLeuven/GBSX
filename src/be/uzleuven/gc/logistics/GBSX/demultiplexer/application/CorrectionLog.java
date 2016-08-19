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
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class CorrectionLog {
    
    /**
     * used log file
     */
    private LoggerFile logFile;
    /**
     * the used samples
     */
    private HashSet<Sample> sampleSet;
    private HashMap<Sample, Integer> correctMapTrimOk;
    private HashMap<Sample, Integer> correctMapR1Corrected;
    private HashMap<Sample, Integer> correctMapR2Corrected;
    private HashMap<Sample, Integer> correctMapR1Not;
    private HashMap<Sample, Integer> correctMapR2Not;
    
    private HashMap<Sample, Integer> trimMapNotR1NotR2;
    private HashMap<Sample, Integer> trimMapNotR1TrimR2butOK;
    private HashMap<Sample, Integer> trimMapNotR1TrimR2butFail;
    private HashMap<Sample, Integer> trimMapTrimR1NotR2Fail;
    private HashMap<Sample, Integer> trimMapTrimR1trimR2longR1;
    private HashMap<Sample, Integer> trimMapTrimR1trimR2longR2;
    private HashMap<Sample, Integer> trimMapTrimR1trimR2ok;
    
    private ReentrantLock lock = new ReentrantLock();
    
    /**
     * creates a new correction log
     * @param logFile 
     */
    public CorrectionLog(LoggerFile logFile){
        this.logFile = logFile;
        this.sampleSet = new HashSet();
        this.correctMapTrimOk = new HashMap();
        this.correctMapR1Corrected = new HashMap();
        this.correctMapR2Corrected = new HashMap();
        this.correctMapR1Not = new HashMap();
        this.correctMapR2Not = new HashMap();
        this.trimMapNotR1NotR2 = new HashMap();
        this.trimMapNotR1TrimR2butFail = new HashMap();
        this.trimMapNotR1TrimR2butOK = new HashMap();
        this.trimMapTrimR1NotR2Fail = new HashMap();
        this.trimMapTrimR1trimR2longR1 = new HashMap();
        this.trimMapTrimR1trimR2longR2 = new HashMap();
        this.trimMapTrimR1trimR2ok = new HashMap();
    }
    
    /**
     * creates a new correction log with the given samples
     * @param samples
     * @param logFile 
     */
    public CorrectionLog(Collection<Sample> samples, LoggerFile logFile){
        this.logFile = logFile;
        this.sampleSet = new HashSet(samples);
        this.correctMapTrimOk = new HashMap();
        this.correctMapR1Corrected = new HashMap();
        this.correctMapR2Corrected = new HashMap();
        this.correctMapR1Not = new HashMap();
        this.correctMapR2Not = new HashMap();
        this.trimMapNotR1NotR2 = new HashMap();
        this.trimMapNotR1TrimR2butFail = new HashMap();
        this.trimMapNotR1TrimR2butOK = new HashMap();
        this.trimMapTrimR1NotR2Fail = new HashMap();
        this.trimMapTrimR1trimR2longR1 = new HashMap();
        this.trimMapTrimR1trimR2longR2 = new HashMap();
        this.trimMapTrimR1trimR2ok = new HashMap();
    }
    
    /**
     * adds all samples to the correction log
     * @param samples 
     */
    public void addSamples(Collection<Sample> samples){
        this.sampleSet.addAll(samples);
    }
    
    /**
     * adds a correct trim for the sample
     * @param sample 
     */
    public void addCorrecterTrimOk(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.correctMapTrimOk.containsKey(sample)){
                this.correctMapTrimOk.put(sample, this.correctMapTrimOk.get(sample) + 1);
            }else{
                this.correctMapTrimOk.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a corrected R1 for the sample
     * @param sample 
     */
    public void addCorrecterR1Corrected(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.correctMapR1Corrected.containsKey(sample)){
                this.correctMapR1Corrected.put(sample, this.correctMapR1Corrected.get(sample) + 1);
            }else{
                this.correctMapR1Corrected.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a corrected R2 for the sample
     * @param sample 
     */
    public void addCorrecterR2Corrected(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.correctMapR2Corrected.containsKey(sample)){
                this.correctMapR2Corrected.put(sample, this.correctMapR2Corrected.get(sample) + 1);
            }else{
                this.correctMapR2Corrected.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a R1 that isn't corrected for the sample
     * @param sample 
     */
    public void addCorrecterR1NotCorrected(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.correctMapR1Not.containsKey(sample)){
                this.correctMapR1Not.put(sample, this.correctMapR1Not.get(sample) + 1);
            }else{
                this.correctMapR1Not.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a R2 that isn't corrected for the sample
     * @param sample 
     */
    public void addCorrecterR2NotCorrected(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.correctMapR2Not.containsKey(sample)){
                this.correctMapR2Not.put(sample, this.correctMapR2Not.get(sample) + 1);
            }else{
                this.correctMapR2Not.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a trim where not is trimmed for the sample
     * @param sample 
     */
    public void addTrimNotR1NotR2(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.trimMapNotR1NotR2.containsKey(sample)){
                this.trimMapNotR1NotR2.put(sample, this.trimMapNotR1NotR2.get(sample) + 1);
            }else{
                this.trimMapNotR1NotR2.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a trim where R1 is not trimmed, R2 is but a fail for the sample
     * @param sample 
     */
    public void addTrimNotR1TrimR2butFail(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.trimMapNotR1TrimR2butFail.containsKey(sample)){
                this.trimMapNotR1TrimR2butFail.put(sample, this.trimMapNotR1TrimR2butFail.get(sample) + 1);
            }else{
                this.trimMapNotR1TrimR2butFail.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a trim where R1 is not trimmed, R2 is but its ok for the sample
     * @param sample 
     */
    public void addTrimNotR1TrimR2butOk(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.trimMapNotR1TrimR2butOK.containsKey(sample)){
                this.trimMapNotR1TrimR2butOK.put(sample, this.trimMapNotR1TrimR2butOK.get(sample) + 1);
            }else{
                this.trimMapNotR1TrimR2butOK.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a trim where R1 is trimmed, not R2 so a fail for the sample
     * @param sample 
     */
    public void addTrimTrimR1NotR2Fail(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.trimMapTrimR1NotR2Fail.containsKey(sample)){
                this.trimMapTrimR1NotR2Fail.put(sample, this.trimMapTrimR1NotR2Fail.get(sample) + 1);
            }else{
                this.trimMapTrimR1NotR2Fail.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a trim where R1 is trimmed, so is R2, but R1 is longer for the sample
     * @param sample 
     */
    public void addTrimTrimR1TrimR2longR1(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.trimMapTrimR1trimR2longR1.containsKey(sample)){
                this.trimMapTrimR1trimR2longR1.put(sample, this.trimMapTrimR1trimR2longR1.get(sample) + 1);
            }else{
                this.trimMapTrimR1trimR2longR1.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a trim where R1 is trimmed, so is R2, but R2 is longer for this sample
     * @param sample 
     */
    public void addTrimTrimR1TrimR2longR2(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.trimMapTrimR1trimR2longR2.containsKey(sample)){
                this.trimMapTrimR1trimR2longR2.put(sample, this.trimMapTrimR1trimR2longR2.get(sample) + 1);
            }else{
                this.trimMapTrimR1trimR2longR2.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * adds a trim where R1 is trimmed, so is R2, and it is correct for the sample
     * @param sample 
     */
    public void addTrimTrimR1TrimR2ok(Sample sample){
        try{
            lock.lock();
            if (! this.sampleSet.contains(sample)){
                this.sampleSet.add(sample);
            }
            if (this.trimMapTrimR1trimR2ok.containsKey(sample)){
                this.trimMapTrimR1trimR2ok.put(sample, this.trimMapTrimR1trimR2ok.get(sample) + 1);
            }else{
                this.trimMapTrimR1trimR2ok.put(sample, 1);
            }
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * write the correction log to the given output
     * @param outputFile 
     */
    public void writeToFile(String outputFile){
        try{
            lock.lock();
            String statsString = "Corrections:" + "\n";
            statsString += "SampleID" + "\t" + "CorrectTrim" + "\t" + "CorrectedR1" + "\t" + "CorrectedR2" + "\t" + "NotCorrectR1" + "\t" + "NotCorrectR2" + "\n";
            for (Sample sample : this.sampleSet){
                statsString += sample.getSampleID() + "\t";
                if (this.correctMapTrimOk.containsKey(sample)){
                    statsString += this.correctMapTrimOk.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.correctMapR1Corrected.containsKey(sample)){
                    statsString += this.correctMapR1Corrected.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.correctMapR2Corrected.containsKey(sample)){
                    statsString += this.correctMapR2Corrected.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.correctMapR1Not.containsKey(sample)){
                    statsString += this.correctMapR1Not.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.correctMapR2Not.containsKey(sample)){
                    statsString += this.correctMapR2Not.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\n";
            }

            statsString += "\n\n\n";
            statsString += "Trimming:" + "\n";
            statsString += "SampleID" + "\t" + "not_R1_not_R2" + "\t" + "not_R1_but_OK" + "\t" + "not_R1_but_fail" + "\t" + "trim_R1_not_R2" + "\t" 
                    + "trim_R1_trim_R2_long_R1" + "\t" + "trim_R1_trim_R2_long_R2" + "\t" + "trim_R1_trim_R2_ok" + "\n";
            for (Sample sample : this.sampleSet){
                statsString += sample.getSampleID() + "\t";
                if (this.trimMapNotR1NotR2.containsKey(sample)){
                    statsString += this.trimMapNotR1NotR2.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.trimMapNotR1TrimR2butOK.containsKey(sample)){
                    statsString += this.trimMapNotR1TrimR2butOK.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.trimMapNotR1TrimR2butFail.containsKey(sample)){
                    statsString += this.trimMapNotR1TrimR2butFail.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.trimMapTrimR1NotR2Fail.containsKey(sample)){
                    statsString += this.trimMapTrimR1NotR2Fail.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.trimMapTrimR1trimR2longR1.containsKey(sample)){
                    statsString += this.trimMapTrimR1trimR2longR1.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.trimMapTrimR1trimR2longR2.containsKey(sample)){
                    statsString += this.trimMapTrimR1trimR2longR2.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\t";
                if (this.trimMapTrimR1trimR2ok.containsKey(sample)){
                    statsString += this.trimMapTrimR1trimR2ok.get(sample);
                }else{
                    statsString += "0";
                }
                statsString += "\n";
            }
            //write the string to a file
            File file = new File(outputFile + System.getProperty("file.separator") + "correction.stats");
            try {
                file.createNewFile();
                BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file));
                bufferedWriter.write(statsString);
                bufferedWriter.close();
            } catch (IOException ex) {
                Logger.getLogger(DemultiplexStats.class.getName()).log(Level.SEVERE, null, ex);
                try {
                    this.logFile.addToLog("Couldn't write the correction stats");
                } catch (ErrorInLogException ex1) {
                    Logger.getLogger(DemultiplexStats.class.getName()).log(Level.SEVERE, null, ex1);
                }
            }
        }finally{
            lock.unlock();
        }
    }
    
}
