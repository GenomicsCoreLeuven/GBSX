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

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class KnuthMorrisPrattDistance {
    
    /**
     * THIS METHOD is based on the HammingsDistance
     * returns false if the length of both sequences is different
     * <br> returns false if the Distance is higher than the maxMismatch
     * <br> this method is faster then HammingsDistance.calculateHammingsDistance(sequence1, sequence2) smaller then maxMismatch
     * <br> the distance is calculated, but when the current distance is higher then the maxMismatch, false is returned
     * <br> (not the whole distance is calculated, except when true is returned, or the last char is also a mismatch)
     * @param sequence1 | String
     * @param sequence2 | String
     * @param maxMismatch | int the max amount of mismatches, where 0 is equals to
     * @return true if the distance is lower or equals the maxMismatch, false otherwise
     */
    public static boolean isEquivalent(String sequence1, String sequence2, int maxMismatch){
        if (sequence1.length() != sequence2.length()){
            return false;
        }
        int distance = 0;
        for (int index = 0; index < sequence1.length(); index++){
            if (sequence1.charAt(index) != sequence2.charAt(index)){
                distance++;
            }
            if (distance > maxMismatch){
                return false;
            }
        }
        return true;
    }
    
    
    
    
    
    /**
     * THIS method is based on the Hammings distance
     * returns -1 if the length of both sequences is different
     * <br> returns -1 if the HammingsDistance is higher than the maxMismatch
     * <br> this method is faster then HammingsDistance.calculateHammingsDistance(sequence1, sequence2) smaller then or equals maxMismatch
     * <br> the distance is calculated, but when the current distance is higher then the maxMismatch, -1 is returned
     * <br> (not the whole distance is calculated, except when a positive int is returned, or the last char is also a mismatch)
     * @param sequence1 | String
     * @param sequence2 | String
     * @param maxMismatch | int the max amount of mismatches, where 0 is equals to
     * @return the hammings distance if the hammings distance is lower or equals the maxMismatch, -1 otherwise
     */
    public static int calculateEquivalentDistance(String sequence1, String sequence2, int maxMismatch){
        if (sequence1.length() != sequence2.length()){
            return -1;
        }
        int distance = 0;
        for (int index = 0; index < sequence1.length(); index++){
            if (sequence1.charAt(index) != sequence2.charAt(index)){
                distance++;
            }
            if (distance > maxMismatch){
                return -1;
            }
        }
        return distance;
    }
    
    
    /**
     * This method is based on the Hammings Distance, Instead of checking each possible possition, it resumes where the previous failed (faster but less acurate Hammings Distance)
     * checks if sequence2 is a substring of sequence1, but allows mismatches.
     * <br> this function is equivalent to the String indexOf, but allows mismatches.
     * <br> the location of sequence2 in sequence1 is returned, if sequence2 (with eventual mismatches) is found in sequence1
     * <br> the mismatches must be under the allowed mismatches (maxMismatch)
     * <br> -1 is returned if sequence2 isn't found
     * @param sequence1 String | the longest string, where sequence2 must be found in
     * @param sequence2 String | the shortest string, where is looked for in sequence1
     * @param maxMismatch int | the maximum allowed mismatches in the comparison
     * @return int the location of the sequence2 in sequence1 (first char is 0),
     */
    public static int indexOf(String sequence1, String sequence2, int maxMismatch){
        //go over every possible location in the first string
        for (int location = 0; location <= sequence1.length() - sequence2.length(); location++){
            //reset the hammings distance
            int distance = 0;
            //check the possible substring
            for (int index = 0; index < sequence2.length(); index++){
                if (sequence1.charAt(location + index) != sequence2.charAt(index)){
                    //if the chars are not equal => add distance
                    distance++;
                }
                if (distance > maxMismatch){
                    //if the distance is higher than the allowed mismatches, stop the comparison soner
                    location += index - 1;
                    break;
                }
                if (index == sequence2.length() - 1 && distance <= maxMismatch){
                    //if the index is at the end of the string, and the distance is lower than or equals the allowed mismatches
                    //than this is a valid subsequence, so return the location
                    return location;
                }
            }
        }
        return -1;
    }
    
    
    
    
    /**
     * This method is based on the Hammings Distance, Instead of checking each possible possition, it resumes where the previous failed (faster but less acurate Hammings Distance)
     * checks if sequence2 is a substring of sequence1, but allows mismathces
     * and if sequence3 is a substring of sequence1, right afther sequence2, allowing mismatches
     * @param sequence1 String | original sequence
     * @param sequence2 String | enzyme cut site
     * @param sequence3 String | adaptor
     * @param maxMismatch2 int | mismatch for sequence2 (enzyme cut site)
     * @param maxMismatch3 int | mismatch for sequence3 (adaptor)
     * @return -1 if not found, the location if sequence2 + sequence3 are found or the location if sequence2 is found and sequence1 is to short to find sequence3
     */
    public static int indexOf(String sequence1, String sequence2, String sequence3, int maxMismatch2, int maxMismatch3){
        //go over every possible location in the first string
        for (int location = 0; location <= sequence1.length() - sequence2.length(); location++){
            //reset the hammings distance
            int distance = 0;
            //check the possible substring
            boolean firstFound = false;
            for (int index = 0; index < sequence2.length(); index++){
                if (sequence1.charAt(location + index) != sequence2.charAt(index)){
                    //if the chars are not equal => add distance
                    distance++;
                }
                if (distance > maxMismatch2){
                    //if the distance is higher than the allowed mismatches, stop the comparison soner
                    location += index - 1;
                    break;
                }
                if (index == sequence2.length() - 1 && distance <= maxMismatch2){
                    //if the index is at the end of the string, and the distance is lower than or equals the allowed mismatches
                    //than this is a valid subsequence, so return the location
                    //return location;
                    firstFound = true;
                }
            }
            if (firstFound){
                if (sequence1.length() - (location + sequence2.length()) <= sequence3.length()){
                    //rest is to smaller then sequence3
                    return location;
                }
                int distance3 = 0;
                for (int index = 0; index < sequence3.length(); index++){
                    if (sequence1.charAt(location + sequence2.length() + index) != sequence3.charAt(index)){
                        distance3++;
                    }
                    if (distance3 > maxMismatch3){
                        location += sequence2.length() + index -1;
                        break;
                    }
                    if (index == sequence3.length() - 1 && distance <= maxMismatch3){
                        return location;
                    }
                }
            }
        }
        return -1;
    }
    
    
}
