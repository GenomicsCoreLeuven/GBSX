/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package be.uzleuven.gc.logistics.GBSX.demultiplexer.model;

import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqRead;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;

/**
 *
 * @author koen
 */
public class ProcessedFragment {

    /**
     * all information of the processed fragment:
     * <br> the sample of the read, 
     * <br> read1 (as HashMap of FastqParts and String)
     * <br> read2 (only pair-end) (as HashMap of FastqParts and String)
     * <br> mismatch occured in finding the barcode/enzyme
     */
        
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