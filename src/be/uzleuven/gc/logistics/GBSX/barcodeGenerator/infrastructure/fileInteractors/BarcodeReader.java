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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeReader {
    
    
    private final BufferedReader bufferedReader;
    private String[] sampleNames;
    
    /**
     * creates a new barcode reader of the given file.
     * @param file | the file that must be read
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public BarcodeReader(File file) throws FileNotFoundException, IOException{
        this.bufferedReader = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(file))));
    }
    
    /**
     * creates a new barcode reader of the given URL.
     * @param url | the URL of the file that must be read
     * @param ziped | boolean if the fastq file is ziped (gz)
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public BarcodeReader(URL url) throws FileNotFoundException, IOException{
        this.bufferedReader = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(url.getPath()))));
    } 
    
    
    /**
     * reads the next lines in the barcode file. and returns it as a string
     * @return null if there are no barcodes anymore, else a barcode
     * @throws IOException | if any error occures while reading the file
     */
    public String next() throws IOException{
        String line = this.bufferedReader.readLine();
        if (line == null){
            return null;
        }
        return line;
    }
    
    /**
     * reads all lines in the barcode file. and returns it as a set of barcodes
     * @return hashset of strings | all barcodes in a set
     * @throws IOException 
     */
    public HashSet<String> readAll() throws IOException{
        HashSet<String> barcodes = new HashSet();
        for (String barcode = this.next(); barcode != null; barcode = this.next()){
            barcodes.add(barcode);
        }
        return barcodes;
    }
    
    /**
     * closes this buffered reader
     * @throws IOException 
     */
    public void close() throws IOException{
        this.bufferedReader.close();
    }
    
    /**
     * closes first the file
     * then execute finilize()
     * @throws Throwable
     */
    @Override
    public void finalize() throws Throwable{
        try {
            this.bufferedReader.close();
        } catch (IOException ex) {
            Logger.getLogger(BarcodeReader.class.getName()).log(Level.SEVERE, null, ex);
        }finally{
            super.finalize();
        }
    }
    
    
}
