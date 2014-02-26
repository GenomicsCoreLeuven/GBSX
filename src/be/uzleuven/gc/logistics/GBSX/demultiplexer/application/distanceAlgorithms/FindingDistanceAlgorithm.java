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
package be.uzleuven.gc.logistics.GBSX.demultiplexer.application.distanceAlgorithms;

import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.FindingsAlgorithms;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FindingDistanceAlgorithm {
    
    
    private final FindingsAlgorithms algorithm;
    
    public FindingDistanceAlgorithm(FindingsAlgorithms findingsAlgorithm){
        this.algorithm = findingsAlgorithm;
    }
    
    /**
     * 
     * @return FindingsAlgorithms | the used algorithm
     */
    public FindingsAlgorithms getFindingsAlgorithm(){
        return this.algorithm;
    }
    
    /**
     * returns false if the length of both sequences is different
     * <br> returns false if the distance is higher than the maxMismatch
     * <br> (not the whole distance is calculated, except when true is returned, or the last char is also a mismatch)
     * @param sequence1 | String the sequence to search in
     * @param sequence2 | String the original (like barcode, enzyme site)
     * @param maxMismatch | int the max amount of mismatches, where 0 is equals to
     * @return true if the distance of the search algorithm is lower or equals the maxMismatch, false otherwise
     */
    public boolean isEquivalent(String sequence1, String sequence2, int maxMismatch){
        if (this.algorithm == FindingsAlgorithms.HAMMINGS_DISTANCE){
            return HammingsDistance.isEquivalent(sequence1, sequence2, maxMismatch);
        }
        if (this.algorithm == FindingsAlgorithms.INDEL_MISMATCH_DISTANCE){
            return IndelMismatchDistance.isEquivalent(sequence1, sequence2, maxMismatch);
        }
        if (this.algorithm == FindingsAlgorithms.KNUTH_MORRIS_PRATT_DISTANCE){
            return KnuthMorrisPrattDistance.isEquivalent(sequence1, sequence2, maxMismatch);
        }
        if (this.algorithm == FindingsAlgorithms.MISMATCH_INDEL_DISTANCE){
            return MismatchIndelDistance.isEquivalent(sequence1, sequence2, maxMismatch);
        }
        return false;
    }
    
    /**
     * returns -1 if the length of both sequences is different
     * <br> returns -1 if the distance is higher than the maxMismatch
     * <br> the distance is calculated, but when the current distance is higher then the maxMismatch, -1 is returned
     * <br> (not the whole distance is calculated, except when a positive int is returned, or the last char is also a mismatch)
     * @param sequence1 | String the sequence to search in
     * @param sequence2 | String the original (like barcode, enzyme site)
     * @param maxMismatch | int the max amount of mismatches, where 0 is equals to
     * @return int[distance, length] | distance: the distance if the distance is lower or equals the maxMismatch, -1 otherwise, length: length of sequence2 (with indels)
     */
    public int[] calculateEquivalentDistance(String sequence1, String sequence2, int maxMismatch){
        if (this.algorithm == FindingsAlgorithms.HAMMINGS_DISTANCE){
            return new int[] {HammingsDistance.calculateEquivalentDistance(sequence1.substring(0, sequence2.length()), sequence2, maxMismatch), sequence2.length()};
        }
        if (this.algorithm == FindingsAlgorithms.INDEL_MISMATCH_DISTANCE){
            return IndelMismatchDistance.calculateEquivalentDistance(sequence1, sequence2, maxMismatch);
        }
        if (this.algorithm == FindingsAlgorithms.KNUTH_MORRIS_PRATT_DISTANCE){
            return new int[] {KnuthMorrisPrattDistance.calculateEquivalentDistance(sequence1.substring(0, sequence2.length()), sequence2, maxMismatch), sequence2.length()};
        }
        if (this.algorithm == FindingsAlgorithms.MISMATCH_INDEL_DISTANCE){
            return MismatchIndelDistance.calculateEquivalentDistance(sequence1, sequence2, maxMismatch);
        }
        return new int[]{-1, 0};
    }
    
    /**
     * checks if sequence2 is a substring of sequence1, but allows mismatches.
     * <br> this function is equivalent to the String indexOf, but allows mismatches.
     * <br> the location of sequence2 in sequence1 is returned, if sequence2 (with eventual mismatches) is found in sequence1
     * <br> the mismatches must be under the allowed mismatches (maxMismatch)
     * <br> -1 is returned if sequence2 isn't found
     * @param sequence1 String | the longest string, where sequence2 must be found in
     * @param sequence2 String | the shortest string, where is looked for in sequence1
     * @param maxMismatch int | the maximum allowed mismatches in the comparison
     * @return int[location, length] | location: the location of the sequence2 in sequence1 (first char is 0), length: length of sequence2 found (with indels)
     */
    public int[] indexOf(String sequence1, String sequence2, int maxMismatch){
        if (this.algorithm == FindingsAlgorithms.HAMMINGS_DISTANCE){
            return new int[] {HammingsDistance.indexOf(sequence1, sequence2, maxMismatch), sequence2.length()};
        }
        if (this.algorithm == FindingsAlgorithms.INDEL_MISMATCH_DISTANCE){
            return IndelMismatchDistance.indexOf(sequence1, sequence2, maxMismatch);
        }
        if (this.algorithm == FindingsAlgorithms.KNUTH_MORRIS_PRATT_DISTANCE){
            return new int[] {KnuthMorrisPrattDistance.indexOf(sequence1, sequence2, maxMismatch), sequence2.length()};
        }
        if (this.algorithm == FindingsAlgorithms.MISMATCH_INDEL_DISTANCE){
            return MismatchIndelDistance.indexOf(sequence1, sequence2, maxMismatch);
        }
        return new int[] {-1, 0};
    }
    
    /**
     * checks if sequence2 is a substring of sequence1, but allows mismathces
     * and if sequence3 is a substring of sequence1, right afther sequence2, allowing mismatches
     * @param sequence1 String | original sequence
     * @param sequence2 String | enzyme cut site
     * @param sequence3 String | adaptor or barcode
     * @param maxMismatch2 int | mismatch for sequence2 (enzyme cut site)
     * @param maxMismatch3 int | mismatch for sequence3 (adaptor)
     * @return int[location, length] | location: -1 if not found, the location if sequence2 + sequence3 are found or the location if sequence2 is found and sequence1 is to short to find sequence3 length: the length of both sequences when found (or with indels)
     */
    public int[] indexOf(String sequence1, String sequence2, String sequence3, int maxMismatch2, int maxMismatch3){
        if (this.algorithm == FindingsAlgorithms.HAMMINGS_DISTANCE){
            return new int[]{HammingsDistance.indexOf(sequence1, sequence2, sequence3, maxMismatch2, maxMismatch3), sequence2.length() + sequence3.length()};
        }
        if (this.algorithm == FindingsAlgorithms.INDEL_MISMATCH_DISTANCE){
            return IndelMismatchDistance.indexOf(sequence1, sequence2, sequence3, maxMismatch2, maxMismatch3);
        }
        if (this.algorithm == FindingsAlgorithms.KNUTH_MORRIS_PRATT_DISTANCE){
            return new int[] {KnuthMorrisPrattDistance.indexOf(sequence1, sequence2, sequence3, maxMismatch2, maxMismatch3), sequence2.length() + sequence3.length()};
        }
        if (this.algorithm == FindingsAlgorithms.MISMATCH_INDEL_DISTANCE){
            return MismatchIndelDistance.indexOf(sequence1, sequence2, sequence3, maxMismatch2, maxMismatch3);
        }
        return new int[]{-1, 0};
    }
    
    
}
