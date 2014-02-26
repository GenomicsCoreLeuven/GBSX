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
package be.uzleuven.gc.logistics.GBSX.barcodeCorrector.application;

import java.util.HashMap;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeCorrectingAlgorithm {
    
    
    /**
     * 
     * @param barcode String | the barcode to check
     * @return true if the hamming code is correct, false otherwise
     */
    public boolean isCorrectBarcode(String barcode){
        char[] charArray = barcode.toCharArray();
        char d1 = charArray[2];
        char d2 = charArray[4];
        char d3 = charArray[5];
        char d4 = charArray[6];
        char d5 = charArray[8];
        char d6 = charArray[9];
        char d7 = charArray[10];
        char d8 = charArray[11];
        char d9 = charArray[12];
        char d10 = charArray[13];
        char d11 = charArray[14];
        char p1 = charArray[0];
        char p2 = charArray[1];
        char p3 = charArray[3];
        char p4 = charArray[7];
        HashMap<Character, Integer> basesIntMap = new HashMap();
        basesIntMap.put('A', 0);
        basesIntMap.put('C', 1);
        basesIntMap.put('G', 2);
        basesIntMap.put('T', 3);
        int predictedP1 = (4 - (basesIntMap.get(d1) + basesIntMap.get(d2) + basesIntMap.get(d4) + basesIntMap.get(d5) 
                + basesIntMap.get(d7) + basesIntMap.get(d9) + basesIntMap.get(d11)) % 4) % 4;
        int predictedP2 = (4 - (basesIntMap.get(d1) + basesIntMap.get(d3) + basesIntMap.get(d4) + basesIntMap.get(d6) 
                + basesIntMap.get(d7) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        int predictedP3 = (4 - (basesIntMap.get(d2) + basesIntMap.get(d3) + basesIntMap.get(d4) + basesIntMap.get(d8) 
                + basesIntMap.get(d9) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        int predictedP4 = (4 - (basesIntMap.get(d5) + basesIntMap.get(d6) + basesIntMap.get(d7) + basesIntMap.get(d8) 
                + basesIntMap.get(d9) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        char base1 = 'N';
        char base2 = 'N';
        char base3 = 'N';
        char base4 = 'N';
        for (char base : basesIntMap.keySet()){
            if (basesIntMap.get(base) == predictedP1) base1 = base;
            if (basesIntMap.get(base) == predictedP2) base2 = base;
            if (basesIntMap.get(base) == predictedP3) base3 = base;
            if (basesIntMap.get(base) == predictedP4) base4 = base;
        }
        if (base1 == p1 && base2 == p2 && base3 == p3 && base4 == p4){
            return true;
        }
        return false;
    }
    
    /**
     * tries to correct the barcode
     * <br> if the barcode is correct, the barcode is returned
     * <br> if an error is found, this is tried to be fixed
     * <br> barcodes with 1 error are fixed correctly
     * <br> barcodes with more then 1 errors can't be fixed (a wrong fix is returned)
     * @param barcode String
     * @return String | the corrected barcode
     */
    public String correctBarcode(String barcode){
        char[] charArray = barcode.toCharArray();
        char d1 = charArray[2];
        char d2 = charArray[4];
        char d3 = charArray[5];
        char d4 = charArray[6];
        char d5 = charArray[8];
        char d6 = charArray[9];
        char d7 = charArray[10];
        char d8 = charArray[11];
        char d9 = charArray[12];
        char d10 = charArray[13];
        char d11 = charArray[14];
        char p1 = charArray[0];
        char p2 = charArray[1];
        char p3 = charArray[3];
        char p4 = charArray[7];
        HashMap<Character, Integer> basesIntMap = new HashMap();
        basesIntMap.put('A', 0);
        basesIntMap.put('C', 1);
        basesIntMap.put('G', 2);
        basesIntMap.put('T', 3);
        int predictedP1 = (4 - (basesIntMap.get(d1) + basesIntMap.get(d2) + basesIntMap.get(d4) + basesIntMap.get(d5) 
                + basesIntMap.get(d7) + basesIntMap.get(d9) + basesIntMap.get(d11)) % 4) % 4;
        int predictedP2 = (4 - (basesIntMap.get(d1) + basesIntMap.get(d3) + basesIntMap.get(d4) + basesIntMap.get(d6) 
                + basesIntMap.get(d7) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        int predictedP3 = (4 - (basesIntMap.get(d2) + basesIntMap.get(d3) + basesIntMap.get(d4) + basesIntMap.get(d8) 
                + basesIntMap.get(d9) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        int predictedP4 = (4 - (basesIntMap.get(d5) + basesIntMap.get(d6) + basesIntMap.get(d7) + basesIntMap.get(d8) 
                + basesIntMap.get(d9) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        char base1 = 'N';
        char base2 = 'N';
        char base3 = 'N';
        char base4 = 'N';
        for (char base : basesIntMap.keySet()){
            if (basesIntMap.get(base) == predictedP1) base1 = base;
            if (basesIntMap.get(base) == predictedP2) base2 = base;
            if (basesIntMap.get(base) == predictedP3) base3 = base;
            if (basesIntMap.get(base) == predictedP4) base4 = base;
        }
        if (base1 == p1 && base2 == p2 && base3 == p3 && base4 == p4){
            return barcode;
        }else{
            //calculates the position of the error
            int place = 0;
            if (predictedP4 != 0){
                place += 8;
            }
            if (predictedP3 != 0){
                place += 4;
            }
            if (predictedP2 != 0){
                place += 2;
            }
            if (predictedP1 != 0){
                place += 1;
            }
            //calculates the right value at that position
            int errorType = Math.max(Math.max(predictedP1, predictedP2), Math.max(predictedP3, predictedP4));
            int correcting = (basesIntMap.get(charArray[place - 1]) - errorType) % 4;
            if (correcting < 0){
                correcting = 4 + correcting;
            }
            char correctBase = 'N';
            for (char base : basesIntMap.keySet()){
                if (basesIntMap.get(base) == correcting) correctBase = base;
            }
            charArray[place - 1] = correctBase;
            String newBarcode = "";
            for (int i = 0; i < charArray.length; i++){
                newBarcode += charArray[i];
            }
            barcode = newBarcode;
        }
        return barcode;
    }
}
