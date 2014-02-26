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
public class MismatchIndelDistance {
   
    
    
    
    
    /**
     * Calculates the distance between 2 sequences, a mismatch is a distancepoint, same as an insert or deletion.
     * @param sequence1 String | first sequence
     * @param sequence2 String | second sequence
     * @return int the distance between the 2 sequences
     */
    public static int calculateIndelMismatchDistance(String sequence1, String sequence2){
        int distance = 0;
        if (sequence1.length() == 1 || sequence2.length() == 1){
            return Math.max(sequence1.length(), sequence2.length());
        }
        int regularDistance = IndelMismatchDistance.calculateIndelMismatchDistance(sequence1.substring(1), sequence2.substring(1));
        if (sequence1.charAt(0) != sequence2.charAt(0)){
            regularDistance++;
        }
        int insertDistance = IndelMismatchDistance.calculateIndelMismatchDistance(sequence1, sequence2.substring(1)) + 1;
        int deleteDistance = IndelMismatchDistance.calculateIndelMismatchDistance(sequence1.substring(1), sequence2) + 1;
        distance = Math.min(regularDistance, insertDistance);
        distance = Math.min(distance, deleteDistance);
        return distance;
    }
      
    /**
     * returns true if the calculateEquivalentDistance function doesn't return -1
     * <br> so sequence2 is completely found in sequence1 (starting from the begining)
     * <br> with maximum mismatch times a mismatch
     * @param sequence1 String | the original sequence (for the indels, longer than the sequence2)
     * @param sequence2 String | the search sequence (like the barcode, enzyme or adaptor)
     * @param mismatch int | number of mismatches that may occure
     * @return true if sequence2 can be found in the beginning of sequence1
     */
    public static boolean isEquivalent(String sequence1, String sequence2, int mismatch){
        if (IndelMismatchDistance.calculateEquivalentDistance(sequence1, sequence2, mismatch)[0] != -1){
            return true;
        }else{
            return false;
        }
    }
    
    
    /**
     * Compare the begining of sequence1 with sequence2, in sequence1 can occure indels to match sequence2
     * <br> there can only be mismatch number of mismatches, then -1 is returned.
     * @param sequence1 String | the original sequence (possible much longer then sequence2)
     * @param sequence2 String | the search sequence (like the barcode, enzyme, adaptor)
     * @param mismatch int | number of mismatches that may occure
     * @return int [distance, length] | for distance: -1 if the number of mismatches > mismatch, else the number of mismatches found for length: the found length
     * @see FindingAlgorithm#calculateEquivalentDistance(java.lang.String, java.lang.String, int) 
     */
    public static int[] calculateEquivalentDistance(String sequence1, String sequence2, int mismatch){
        if (mismatch < 0){
            //mismatch is lower then zero
            return new int[]{-1, 0};
        }
        if (sequence2.length() == 0){
            //end of search sequence
            return new int[]{0, 0};
        }
        if (sequence1.length() == 0){
            if (sequence2.length() < mismatch){
                return new int[] {sequence2.length(), 0};
            }
            //end of original sequence
            return new int[] {-1, 0};
        }
        //end of strings
        if (sequence2.length() == 1){
            //end of strings
            if (sequence1.charAt(0) == sequence2.charAt(0)){
                //equal => distance is 0
                return new int[]{0, 1};
            }else{
                //not equal => mismatch
                mismatch--;
                if (mismatch < 0){
                    //more mismatches than may occure
                    return new int[]{-1, 0};
                }else{
                    //one mismatch
                    return new int[]{1, 1};
                }
            }
        }
        //not end of string
        int distance = -1;
        int length = 0;
        //calculate regular distance
        int regularMismatch = mismatch;
        int regularDistance = 0;
        if (sequence1.charAt(0) != sequence2.charAt(0)){
            regularDistance = 1;
            regularMismatch--;
        }
        int[] regularOutput = IndelMismatchDistance.calculateEquivalentDistance(sequence1.substring(1), sequence2.substring(1), regularMismatch);
        if (regularOutput[0] == -1){
            regularDistance = -1;
        }else{
            regularDistance += regularOutput[0];
        }
        int regularLength = regularOutput[1] + 1;
        
        
        //if a regular distance is found, do not look for a better indel distance
        if (regularDistance != -1 && regularDistance <= mismatch){
            distance = regularDistance;
            length = regularLength;
        }else{
            //calculate insert distance
            int insertDistance = -1;
            int insertMismatch = mismatch - 1;
            int[] insertOuput = IndelMismatchDistance.calculateEquivalentDistance(sequence1.substring(1), sequence2, insertMismatch);
            int insertLenght = insertOuput[1] + 1;
            if (insertOuput[0] != -1) insertDistance = 1 + insertOuput[0];

            //calculate delete distance
            int deleteDistance = -1;
            int deleteMismatch = mismatch - 1;
            int[] deleteOutput = IndelMismatchDistance.calculateEquivalentDistance(sequence1, sequence2.substring(1), deleteMismatch);
            int deleteLength = deleteOutput[1];
            if (deleteOutput[0] != -1) deleteDistance = deleteOutput[0] + 1;
            
            if (insertDistance != -1){
                distance = insertDistance;
                length = insertLenght;
            }
            if (deleteDistance != -1){
                if (distance == -1){
                    distance = deleteDistance;
                    length = deleteLength;
                }else{
                    if (distance > deleteDistance){
                        distance = deleteDistance;
                        length = deleteLength;
                    }
                }
            }
        }
        
        if (distance > mismatch){
            return new int[]{-1,0};
        }
        return new int[]{distance, length};
    }
    
    
    /**
     * 
     * @param sequence1 String | the original sequence
     * @param sequence2 String | the sequence where is searched for
     * @param mismatch int | the number of mismatches that may occure
     * @return int[distance, length] | distance: -1 if no possible match is found, else the index of sequence2 in sequence1 length: the length of the found piece
     */
    public static int[] indexOf(String sequence1, String sequence2, int mismatch){
        //go over every possible location in the first string
        for (int location = 0; location <= sequence1.length() - sequence2.length(); location++){
            //calculate distance
            int[] distanceLenght = IndelMismatchDistance.calculateEquivalentDistance(sequence1.substring(location), sequence2, mismatch);
            if (distanceLenght[0] != -1){
                return new int[]{location, distanceLenght[1]};
            }
        }
        return new int[]{-1, 0};
    }
    
    /**
     * 
     * @param sequence1 String | the original sequence
     * @param sequence2 String | the first sequence to find in the original
     * @param sequence3 String | the second sequence to find right after sequence2 in the original sequence
     * @param mismatch2 int | number of mismatches possible in sequence2
     * @param mismatch3 int | number of mismatches possible in sequence3
     * @return int[location, lenght] | location: -1 if not found, else start location, length: the length of the found sequences (sequence2.length + sequence3.length, longer when inserts, shorter when deletes)
     */
    public static int[] indexOf(String sequence1, String sequence2, String sequence3, int mismatch2, int mismatch3){
        //go over every possible location in the first string
        for (int location = 0; location <= sequence1.length() - sequence2.length(); location++){
            //calculate distance
            int[] distanceLength = IndelMismatchDistance.calculateEquivalentDistance(sequence1.substring(location), sequence2, mismatch2);
            if (distanceLength[0] != -1){
                //sequence2 found
                if ((sequence1.length() - (location + sequence2.length())) < sequence3.length()){
                    //sequence1 is not long enough to find sequence3
                    return new int[]{location, distanceLength[1]};
                }else{
                    int[] distanceLength3 = IndelMismatchDistance.calculateEquivalentDistance(sequence1.substring(location + distanceLength[1]), sequence3, mismatch3);
                    if (distanceLength3[0] != -1){
                        //sequence3 found
                        return new int[]{location, distanceLength[1] + distanceLength3[1]};
                    }
                }
            }
        }
        return new int[]{-1,0};
    }
    
    
    
    
}
