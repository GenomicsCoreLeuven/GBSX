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
package be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model;

import java.util.ArrayList;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BasePair {
    
    /**
     * returns true if the the found base can be the espected base
     * @param base1 char | the espected base (A, C, G or T)
     * @param base2 char | the found base (A,C,G,T,R,Y,K,M,S,W,B,D,H,V,N)
     * @return true if base2 possible is base1
     */
    public static boolean canBeEqual(char base1, char base2){
        base1 = Character.toUpperCase(base1);
        for (char possibleBase : BasePair.getPossibleBasePairs(base2)){
            if (possibleBase == base1){
                return true;
            }
        }
        return false;
    }
    
    /**
     * transforms the given base to a list of possible bases (usefull for R,Y,...)
     * @param base char | the base to transform
     * @return arraylist of character | the possible corresponding bases.
     */
    private static ArrayList<Character> getPossibleBasePairs(char base){
        base = Character.toUpperCase(base);
        ArrayList<Character> base1List = new ArrayList();
        switch (base){
            case 'A' : base1List.add('A');
                break;
            case 'C' : base1List.add('C');
                break;
            case 'G' : base1List.add('G');
                break;
            case 'T' : base1List.add('T');
                break;
            case 'R' : base1List.add('A');
                base1List.add('G');
                break;
            case 'Y' : base1List.add('T');
                base1List.add('G');
                break;
            case 'K' : base1List.add('G');
                base1List.add('T');
                break;
            case 'M' : base1List.add('A');
                base1List.add('C');
                break;
            case 'S' : base1List.add('G');
                base1List.add('C');
                break;
            case 'W' : base1List.add('A');
                base1List.add('T');
                break;
            case 'B' : base1List.add('G');
                base1List.add('T');
                base1List.add('C');
                break;
            case 'D' : base1List.add('G');
                base1List.add('A');
                base1List.add('T');
                break;
            case 'H' : base1List.add('A');
                base1List.add('C');
                base1List.add('T');
                break;
            case 'V' : base1List.add('G');
                base1List.add('C');
                base1List.add('A');
                break;
            case 'N' : base1List.add('A');
                base1List.add('G');
                base1List.add('C');
                base1List.add('T');
                break;
            default : base1List.add('A');
                base1List.add('G');
                base1List.add('C');
                base1List.add('T');
                break;
        }
        return base1List;
    }
    
    /**
     * creates the complement of the given DNA
     * @param sequence String | the sequence of DNA where the complement must be searched
     * @return the complement of the given DNA
     */
    public static String getComplementSequence(String sequence){
        String complementDNApiece = "";
        //sequence = sequence.toUpperCase();
        for (char base : sequence.toCharArray()){
            switch (base){
                case 'A': complementDNApiece = 'T' + complementDNApiece;
                    break;
                case 'T' : complementDNApiece = 'A' + complementDNApiece;
                    break;
                case 'G' : complementDNApiece = 'C' + complementDNApiece;
                    break;
                case 'C': complementDNApiece = 'G' + complementDNApiece;
                    break;
                case 'a': complementDNApiece = 't' + complementDNApiece;
                    break;
                case 't' : complementDNApiece = 'a' + complementDNApiece;
                    break;
                case 'g' : complementDNApiece = 'c' + complementDNApiece;
                    break;
                case 'c': complementDNApiece = 'g' + complementDNApiece;
                    break;
                case '[' : complementDNApiece = ']' + complementDNApiece;
                    break;
                case ']' : complementDNApiece = '[' + complementDNApiece;
                    break;
                default : complementDNApiece = base + complementDNApiece;
                    break;
            }
        }
        return complementDNApiece;
    }
}
