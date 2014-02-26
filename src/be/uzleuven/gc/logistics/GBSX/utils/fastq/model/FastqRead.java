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
package be.uzleuven.gc.logistics.GBSX.utils.fastq.model;

import java.util.HashMap;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FastqRead {
    
    HashMap<FastqParts, String> parts;
    FastqScores fastqScore = FastqScores.getStandard();
    
    /**
     * creates a new empty read
     */
    private FastqRead(){
        this.parts = new HashMap();
    }
    
    /**
     * creates a new read with only a description
     * @param description String | the description of the fastq read
     */
    private FastqRead(String description){
        this.parts = new HashMap();
        if (! description.startsWith("@")){
            description = "@" + description;
        }
        this.parts.put(FastqParts.DESCRIPTION, description);
    }
    
    /**
     * creates a new complete fastq read
     * @param description String | the description of the fastq read
     * @param sequence String | the sequence of the read
     * @param quality String | the quality of the sequence
     */
    public FastqRead(String description, String sequence, String quality){
        this.parts = new HashMap();
        if (! description.startsWith("@")){
            description = "@" + description;
        }
        this.parts.put(FastqParts.DESCRIPTION, description);
        this.parts.put(FastqParts.SEQUENCE, sequence);
        this.parts.put(FastqParts.QUALITY, quality);
    }
    
    /**
     * creates a new fastq read based on the old fastq map with fastq parts system
     * @param fastqMap 
     */
    public FastqRead(HashMap<FastqParts, String> fastqMap){
        this.parts = new HashMap(fastqMap);
    }
    
    /**
     * @param score FastqScores
     * creates a new empty read
     */
    private FastqRead(FastqScores score){
        this.parts = new HashMap();
        this.fastqScore = score;
    }
    
    /**
     * creates a new read with only a description
     * @param description String | the description of the fastq read
     * @param score FastqScores
     */
    private FastqRead(String description, FastqScores score){
        this.parts = new HashMap();
        if (! description.startsWith("@")){
            description = "@" + description;
        }
        this.parts.put(FastqParts.DESCRIPTION, description);
        this.fastqScore = score;
    }
    
    /**
     * creates a new complete fastq read
     * @param description String | the description of the fastq read
     * @param sequence String | the sequence of the read
     * @param quality String | the quality of the sequence
     * @param score FastqScore
     */
    public FastqRead(String description, String sequence, String quality, FastqScores score){
        this.parts = new HashMap();
        if (! description.startsWith("@")){
            description = "@" + description;
        }
        this.parts.put(FastqParts.DESCRIPTION, description);
        this.parts.put(FastqParts.SEQUENCE, sequence);
        this.parts.put(FastqParts.QUALITY, quality);
        this.fastqScore = score;
    }
    
    /**
     * creates a new fastq read based on the old fastq map with fastq parts system
     * @param fastqMap 
     * @param score
     */
    public FastqRead(HashMap<FastqParts, String> fastqMap, FastqScores score){
        this.parts = new HashMap(fastqMap);
        this.fastqScore = score;
    }
    
    /**
     * 
     * @return String | the sequence
     */
    public String getSequence(){
        return this.parts.get(FastqParts.SEQUENCE);
    }
    
    /**
     * 
     * @return String | the description 
     */
    public String getDescription(){
        return this.parts.get(FastqParts.DESCRIPTION);
    }
    
    /**
     * 
     * @return String | the quality
     */
    public String getQuality(){
        return this.parts.get(FastqParts.QUALITY);
    }
    
    /**
     * set the given sequence as new
     * @param newSequence String | may not be null or empty
     */
    private void setSequence(String newSequence){
        if (newSequence != null || ! newSequence.isEmpty()){
            this.parts.put(FastqParts.SEQUENCE, newSequence);
        }
    }
    
    /**
     * set the given quality as new
     * @param newQuality String | may not be null or empty
     */
    private void setQuality(String newQuality){
        if (newQuality != null || ! newQuality.isEmpty()){
            this.parts.put(FastqParts.QUALITY, newQuality);
        }
    }
    
    /**
     * set the givend description as new
     * @param newDescription String | may not be null or empty
     */
    private void setDescription(String newDescription){
        if (newDescription != null || newDescription.isEmpty()){
            if (! newDescription.startsWith("@")){
                newDescription = "@" + newDescription;
            }
            this.parts.put(FastqParts.DESCRIPTION, newDescription);
        }
    }
    
    /**
     * this method trims at the left of both sequence and quality
     * the given distance is the index of the first base to keep (first is 0)
     * @param distance int | the first base to keep
     */
    public void ltrim(int distance){
        this.parts.put(FastqParts.SEQUENCE, this.parts.get(FastqParts.SEQUENCE).substring(distance));
        this.parts.put(FastqParts.QUALITY, this.parts.get(FastqParts.QUALITY).substring(distance));
    }
    
    /**
     * this method trims at the right of both sequence and quality
     * the given distance is the index of the last base exclusive (first is 0)
     * @param distance int | the first base to delete
     */
    public void rtrim(int distance){
        this.parts.put(FastqParts.SEQUENCE, this.parts.get(FastqParts.SEQUENCE).substring(0, distance));
        this.parts.put(FastqParts.QUALITY, this.parts.get(FastqParts.QUALITY).substring(0, distance));
    }
    
    
}
