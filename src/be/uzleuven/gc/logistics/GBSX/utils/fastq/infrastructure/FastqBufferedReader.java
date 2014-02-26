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
package be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure;

import be.uzleuven.gc.logistics.GBSX.utils.FileLocker;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqRead;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FastqBufferedReader {
    
    private final BufferedReader fastqBufferedReader;
    private FileLocker fileLocker = new FileLocker();
    
    /**
     * creates a new FastqBufferedReader of the given file.
     * @param file | the file that must be read
     * @param ziped | boolean if the fastq file is ziped (gz)
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastqBufferedReader(File file, boolean ziped) throws FileNotFoundException, IOException{
        if (ziped){
            this.fastqBufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
        }else{
            this.fastqBufferedReader = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(file))));
        }
    }
    
    
    /**
     * reads the next lines in the fastq file. and returns it as an fastqread
     * @return null if there are no fastq files anymore, else a FastqRead
     * @throws IOException | if any error occures while reading the file
     */
    public FastqRead next() throws IOException{
        if (this.fileLocker.lock()){
            String descriptionLine = this.fastqBufferedReader.readLine();
            if (descriptionLine == null){
                this.fileLocker.unlock();
                return null;
            }

            String sequenceLine = this.fastqBufferedReader.readLine();
            if (sequenceLine == null){
                this.fileLocker.unlock();
                return null;
            }
            //the + line
            this.fastqBufferedReader.readLine();
            String qualityLine = this.fastqBufferedReader.readLine();
            if (qualityLine == null){
                this.fileLocker.unlock();
                return null;
            }
            FastqRead fastqRead = new FastqRead(descriptionLine, sequenceLine, qualityLine);
            this.fileLocker.unlock();
            return fastqRead;
        }
        return null;
    }  
    
    /**
     * closes this buffered reader
     * @throws IOException 
     */
    public void close() throws IOException{
        this.fileLocker.waitTillCompleteUnlock();
        this.fastqBufferedReader.close();
    }
    
    
    /**
     * closes first the file
     * then execute finilize()
     * @throws Throwable
     */
    @Override
    public void finalize() throws Throwable{
        try {
            this.fastqBufferedReader.close();
        } catch (IOException ex) {
            Logger.getLogger(FastqBufferedReader.class.getName()).log(Level.SEVERE, null, ex);
        }finally{
            super.finalize();
        }
    }
    
}
