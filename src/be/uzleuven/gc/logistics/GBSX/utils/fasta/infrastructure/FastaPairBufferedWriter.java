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
package be.uzleuven.gc.logistics.GBSX.utils.fasta.infrastructure;

import be.uzleuven.gc.logistics.GBSX.utils.fasta.model.FastaRead;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FastaPairBufferedWriter {
    private FastaBufferedWriter fastaBufferedWriter1;
    private FastaBufferedWriter fastaBufferedWriter2;
    private ReentrantLock lock = new ReentrantLock();
    
    /**
     * creates 2 new FastaBufferedWriter of the given files.
     * @param fileRead1 | the file for read 1
     * @param fileRead2 | the file for read 2
     * @param ziped | boolean if the fasta file is ziped (gz)
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastaPairBufferedWriter(File fileRead1, File fileRead2, boolean ziped) throws FileNotFoundException, IOException{
        this.fastaBufferedWriter1 = new FastaBufferedWriter(fileRead1, ziped);
        this.fastaBufferedWriter2 = new FastaBufferedWriter(fileRead2, ziped);
    }
    
    /**
     * reads the next lines in the fasta file. and returns it as an fastaread
     * @return null if there are no fasta files anymore, else a FastqRead
     * @throws IOException | if any error occures while reading the file
     */
    public void write(FastaRead fastaRead1, FastaRead fastaRead2) throws IOException{
        try{
            lock.lock();
            this.fastaBufferedWriter1.write(fastaRead1);
            this.fastaBufferedWriter2.write(fastaRead2);
        }finally{
            lock.unlock();
        }
    }  
    
    /**
     * closes this buffered reader
     * @throws IOException 
     */
    public void close() throws IOException{
        try{
            lock.lock();
            this.fastaBufferedWriter1.close();
            this.fastaBufferedWriter2.close();
        }finally{
            lock.unlock();
        }
    }
    
    
    /**
     * closes first the file
     * then execute finilize()
     * @throws Throwable
     */
    @Override
    public void finalize() throws Throwable{
        try {
            this.close();
        } catch (IOException ex) {
            Logger.getLogger(FastaBufferedReader.class.getName()).log(Level.SEVERE, null, ex);
        }finally{
            super.finalize();
        }
    }
}
