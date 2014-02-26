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
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FastqBufferedWriter {
    
    private final BufferedWriter fastqBufferedWriter;
    private FileLocker fileLocker = new FileLocker();
    
    /**
     * creates a new FastqBufferedWriter of the given file.
     * @param file | the file that must be writen
     * @param ziped | true if the file must be zipped
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastqBufferedWriter(File file, boolean ziped) throws FileNotFoundException, IOException{
        if (ziped){
            this.fastqBufferedWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file))));
        }else{
            this.fastqBufferedWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(file))));
        }
    }
    
    /**
     * creates a new FastqBufferedWriter of the given file.
     * @param file | the file that must be writen
     * @param ziped | true if the file must be zipped
     * @param append | true if the file must append an existing file
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastqBufferedWriter(File file, boolean ziped, boolean append) throws FileNotFoundException, IOException{
        if (ziped){
            this.fastqBufferedWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file, append))));
        }else{
            this.fastqBufferedWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(file, append))));
        }
    }
    
    
    /**
     * writes the next lines in the fastq file. 
     * @param fastq | a fastqRead
     * @throws IOException | if any error occures while writing the file
     */
    public void write(FastqRead fastq) throws IOException{
        if (this.fileLocker.lock()){
            this.fastqBufferedWriter.write(fastq.getDescription());
            this.fastqBufferedWriter.write("\n");
            this.fastqBufferedWriter.write(fastq.getSequence());
            this.fastqBufferedWriter.write("\n");
            this.fastqBufferedWriter.write("+");
            this.fastqBufferedWriter.write("\n");
            this.fastqBufferedWriter.write(fastq.getQuality());
            this.fastqBufferedWriter.write("\n");
            this.fileLocker.unlock();
        }else{
            throw new IOException("FastqBuffer closed");
        }
    }
    
    
    /**
     * closes this bufferedWriter
     * @throws IOException 
     */
    public void close() throws IOException{
        this.fileLocker.waitTillCompleteUnlock();
        this.fastqBufferedWriter.close();
    }
    
    
    /**
     * closes first the file
     * then execute finilize()
     * @throws Throwable
     */
    @Override
    public void finalize() throws Throwable{
        try {
            this.fastqBufferedWriter.close();
        } catch (IOException ex) {
            Logger.getLogger(FastqBufferedWriter.class.getName()).log(Level.SEVERE, null, ex);
        }finally{
            super.finalize();
        }
    }
    
}
