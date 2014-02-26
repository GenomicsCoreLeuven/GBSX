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
package be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions;

import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqRead;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class InvalidReadException extends Exception{
    
    private final FastqRead read1;
    private final FastqRead read2;
    private final InvalidReadEnum invalidRead;
    
    /**
     * creates a new invalid read (pair-end)
     * @param read1 FastqRead | the first read
     * @param read2 FastqRead | the second read
     * @param invalidRead InvalidReadEnum | the invalid read
     */
    public InvalidReadException(FastqRead read1, FastqRead read2, InvalidReadEnum invalidRead){
        this.read1 = read1;
        this.read2 = read2;
        this.invalidRead = invalidRead;
    }
    
    /**
     * creates a new invalid read (single end)
     * @param read1 FastqRead | the first read
     * @param invalidRead InvalidReadEnum | the invalid read
     */
    public InvalidReadException(FastqRead read1, InvalidReadEnum invalidRead){
        this.read1 = read1;
        this.read2 = null;
        this.invalidRead = invalidRead;
    }
    
    /**
     * 
     * @return InvalidReadEnum | returns witch read was invalid 
     */
    public InvalidReadEnum getInvalidRead(){
        return this.invalidRead;
    }
    
    /**
     * 
     * @return FastqRead | the first read (read from R1)
     */
    public FastqRead getRead1(){
        return this.read1;
    }
    
    /**
     * 
     * @return FastqRead | the second read (read from R2)
     */
    public FastqRead getRead2(){
        return this.read2;
    }
    
    /**
     * 
     * @return String | Read1 as full String
     */
    public String getRun1(){
        return this.read1.getDescription() + "\n" + this.read1.getSequence() + "\n" + "+" + "\n" + this.read1.getQuality() + "\n";
    }
    
    /**
     * 
     * @return String | Read2 as full String
     */
    public String getRun2(){
        return this.read2.getDescription() + "\n" + this.read2.getSequence() + "\n" + "+" + "\n" + this.read2.getQuality() + "\n";
    }
    
}
