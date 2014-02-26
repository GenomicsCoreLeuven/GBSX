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
package be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.infrastructure;

import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeCollection;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class InfoFileParser {
    
    private EnzymeCollection enzymeCollection;
    
    /**
     * creates a new InfoFileParser for the parsing of the info file. Uses the given enzyme collection to detect the possible enzymes
     * @param enzymeCollection 
     */
    public InfoFileParser(EnzymeCollection enzymeCollection){
        this.enzymeCollection = enzymeCollection;
    }
    
    
    /**
     * opens the given file and search all samples
     * <br> the file has no header, and the following columns:
     * <br> sampleID
     * <br> barcode
     * <br> enzyme
     * @param filename String | the name of the file
     * @return an arrayList of all samples found in the file. Non found enzymes are skiped
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public ArrayList<Sample> parseInfoFile(String filename) throws FileNotFoundException, IOException{
        ArrayList<Sample> sampleArrayList = new ArrayList();
        File file = new File(filename);
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null){
            String[] splitedLine = line.split("\t");
            Enzyme enzyme = this.enzymeCollection.getEnzyme(splitedLine[2]);
            if (enzyme != null){
                if (splitedLine.length > 3){
                    int mismatch = Integer.parseInt(splitedLine[3]);
                    Sample newSample = new Sample(splitedLine[0], enzyme, splitedLine[1], mismatch);
                    sampleArrayList.add(newSample);
                }else{
                    Sample newSample = new Sample(splitedLine[0], enzyme, splitedLine[1]);
                    sampleArrayList.add(newSample);
                }
            }
        }
        return sampleArrayList;
    }
    
}
