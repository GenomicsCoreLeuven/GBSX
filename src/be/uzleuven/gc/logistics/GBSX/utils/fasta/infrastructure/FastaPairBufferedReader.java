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
import java.util.HashMap;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FastaPairBufferedReader {
    
    private FastaBufferedReader fastaBufferedReader1;
    private FastaBufferedReader fastaBufferedReader2;
    private ReentrantLock lock = new ReentrantLock();
    
    /**
     * creates 2 new FastaBufferedReader of the given files.
     * @param fileRead1 | the file for read 1
     * @param fileRead2 | the file for read 2
     * @param ziped | boolean if the fastq file is ziped (gz)
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastaPairBufferedReader(File fileRead1, File fileRead2, boolean ziped) throws FileNotFoundException, IOException{
        this.fastaBufferedReader1 = new FastaBufferedReader(fileRead1, ziped);
        this.fastaBufferedReader2 = new FastaBufferedReader(fileRead2, ziped);
    }
    
    /**
     * reads the next lines in the fasta file. and returns it as an fastaread
     * @return null if there are no fasta files anymore, else a FastaRead
     * @throws IOException | if any error occures while reading the file
     */
    public HashMap<String, FastaRead> next() throws IOException{
        HashMap<String, FastaRead> tmpMap = new HashMap();
        try{
            lock.lock();
            tmpMap.put("r1", this.fastaBufferedReader1.next());
            tmpMap.put("r2", this.fastaBufferedReader2.next());
        }finally{
            lock.unlock();
        }
        return tmpMap;
    }  
    
    /**
     * closes this buffered reader
     * @throws IOException 
     */
    public void close() throws IOException{
        try{
            lock.lock();
            this.fastaBufferedReader1.close();
            this.fastaBufferedReader2.close();
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
