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
package be.uzleuven.gc.logistics.GBSX.utils.fasta.model;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FastaRead {
    
    private String sequence;
    private String description;
    
    /**
     * creates a new fasta read with "no description" and the given sequence
     * @param sequence | String the sequence
     */
    private FastaRead(String sequence){
        this.description = "> no description";
        this.sequence = sequence;
    }
    
    /**
     * creates a new fasta read with the given description and the given sequence
     * @param sequence | String the sequence
     */
    public FastaRead(String description, String sequence){
        if (! description.startsWith(">")){
            description = ">" + description;
        }
        this.description = description;
        this.sequence = sequence;
    }
    
    /**
     * 
     * @return String | the description of the fasta read
     */
    public String getDescription(){
        return this.description;
    }
    
    /**
     * 
     * @return String | the sequence of the fasta read
     */
    public String getSequence(){
        return this.sequence;
    }
    
}
