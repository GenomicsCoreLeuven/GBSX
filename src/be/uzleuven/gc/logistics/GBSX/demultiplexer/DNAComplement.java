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
package be.uzleuven.gc.logistics.GBSX.demultiplexer;

import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.BasePair;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class DNAComplement {
    
    public final static String VERSION = "DNAComplement v1.0";
    public final static String LICENCE = "GPLv3";
    public final static boolean DEBUG = false;
    
    public static void main (String[] args){
        if (DEBUG){
            args = new String[1];
            args[0] = "CAGCCTGCCCAGTCGAGGGGGAGGCGGGGCGAGGGCTGCGGGCAAGGCCTGGCCGGGGGGTCCCGGCTGGGGCGG";
        }
        if (args.length == 0 || args[0].equals("help") || args[0].equals("-help") || args[0].equals("-h")){
            System.out.println("This is the DNA complement creator.");
            System.out.println("The only parameter is a string of DNA.");
        }else if (args[0].equals("version") || args[0].equals("-version") || args[0].equals("-v")){
            System.out.println(DNAComplement.VERSION);
        }else{
            String DNA = args[0];
            String complement = BasePair.getComplementSequence(DNA);
            System.out.println(complement);
        }
    }
    
    public static String complement(String sequence){
        return BasePair.getComplementSequence(sequence);
    }
    
    
    
}
