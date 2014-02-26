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
package be.uzleuven.gc.logistics.GBSX.barcodeGenerator.infrastructure.fileInteractors;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeWriter {
    
    private String outputDirectory;
    
    public BarcodeWriter(String outputDirectory){
        File outDir = new File(outputDirectory);
        if (! outDir.exists()){
            outDir.mkdirs();
        }
        this.outputDirectory = outputDirectory;
    }
    
    /**
     * 
     * @param barcodeCollection collection of String | all barcodes that must be writen to a file
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void write(Collection<String> barcodeCollection) throws FileNotFoundException, IOException{
        ArrayList<String> barcodes = new ArrayList(barcodeCollection);
        Collections.sort(barcodes);
        //make buffered writer
        File file = new File(this.outputDirectory + System.getProperty("file.separator") + "barcode_list.txt");
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(file))));
        //write every barcode
        for (String b : barcodes){
            writer.write(b + "\n");
        }
        writer.close();
    }
    
    
}
