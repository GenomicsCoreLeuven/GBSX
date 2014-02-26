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
package be.uzleuven.gc.logistics.GBSX.demultiplexer.model;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public enum FindingsAlgorithms {
    
    /**
     * use hammings distance (exact mismatches no indels)
     */
    HAMMINGS_DISTANCE ("hammings", "Checks for mismatches (no indels)"),
    
    /**
     * use hammings distance based algorithm (faster than hammings distance, but less accurate)
     */
    KNUTH_MORRIS_PRATT_DISTANCE ("knuth", "Faster than hammings, but can miss some locations"),
    
    /**
     * search for indels and mismatches and takes everytime the most accurate method (prefering the mismatch, than the insert, than the delete)
     */
    INDEL_MISMATCH_DISTANCE ("indelmis", "Checks for mismatches and indels, the barcode/enzyme/adaptor with the least errors (mismatches or indels) is taken"),
    
    /**
     * search first for mismatches, if no found, look for indels
     */
    MISMATCH_INDEL_DISTANCE ("misindel", "Checks for mismatches and indels, the mismatches are supperior to the indels (faster than indelmis, but errors can be higher)");
    
    /**
     * the name of the argument (-name)
     */
    private final String algorithmName;
    
    private final String description;
    
    /**
     * returns the asked algorithm, or HAMMINGS DISTANCE when not found (standard algorithm)
     * @param algorithmName String | the name of the algorithm
     * @return the algorithm or the standard when not found
     * @see FindingsAlgorithms#getStandard() 
     */
    public static FindingsAlgorithms getAlgorithm(String algorithmName){
        for (FindingsAlgorithms arg : FindingsAlgorithms.values()){
            if (algorithmName.equals(arg.getAlgorithmName())){
                return arg;
            }
        }
        //standard algorithm
        return FindingsAlgorithms.getStandard();
    }
    
    /**
     * 
     * @return the standard algorithm
     */
    public static FindingsAlgorithms getStandard(){
        return FindingsAlgorithms.HAMMINGS_DISTANCE;
    }
    
    /**
     * 
     * @param sortname String | name of the argument
     */
    private FindingsAlgorithms(String algorithmName, String description){
        this.algorithmName = algorithmName;
        this.description = description;
    }
    
    /**
     * 
     * @return String | the name of the algorithm
     */
    public String getAlgorithmName(){
        return this.algorithmName;
    }
    
    /**
     * 
     * @return String | a short description of the algorithm
     */
    public String getDescription(){
        return this.description;
    }
    
}
