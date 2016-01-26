/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package be.uzleuven.gc.logistics.GBSX.utils.fasta.infrastructure;

import be.uzleuven.gc.logistics.GBSX.utils.FileLocker;
import be.uzleuven.gc.logistics.GBSX.utils.fasta.model.FastaRead;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FastaBufferedWriter {
    
    
    private final BufferedWriter fastaBufferedWriter;
    private FileLocker fileLocker = new FileLocker();
    
    /**
     * creates a new FastaBufferedWriter of the given file.
     * @param file | the file that must be writen
     * @param  ziped | true if the file must be zipped
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastaBufferedWriter(File file, boolean ziped) throws FileNotFoundException, IOException{
        if (ziped){
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file))));
        }else{
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(file))));
        }
    }
    
    /**
     * creates a new FastaBufferedWriter of the given file.
     * @param file | the file that must be writen
     * @param  ziped | true if the file must be zipped
     * @param append | true if append to existing file
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastaBufferedWriter(File file, boolean ziped, boolean append) throws FileNotFoundException, IOException{
        if (ziped){
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file, append))));
        }else{
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(file, append))));
        }
    }
    
    /**
     * creates a new FastaBufferedReader of the given URL.
     * @param url | the URL of the file that must be writen
     * @param  ziped | true if the file must be zipped
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastaBufferedWriter(URL url, boolean ziped) throws FileNotFoundException, IOException{
        if (ziped){
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(url.getPath()))));
        }else{
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(url.getPath()))));
        }
    } 
    /**
     * creates a new FastaBufferedReader of the given URL.
     * @param url | the URL of the file that must be writen
     * @param  ziped | true if the file must be zipped
     * @param append | true if append to existing file
     * @throws FileNotFoundException | if the given file isn't found
     * @throws IOException | if an error occures while opening the file
     */
    public FastaBufferedWriter(URL url, boolean ziped, boolean append) throws FileNotFoundException, IOException{
        if (ziped){
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(url.getPath(), append))));
        }else{
            this.fastaBufferedWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(url.getPath(), append))));
        }
    } 
    
    
    /**
     * writes the next lines in the fasta file. 
     * @param fasta | a fastaRead
     * @throws IOException | if any error occures while writing the file
     */
    public void write(FastaRead fasta) throws IOException{
        if (this.fileLocker.lock()){
            this.fastaBufferedWriter.write(fasta.getDescription());
            this.fastaBufferedWriter.write("\n");
            this.fastaBufferedWriter.write(fasta.getSequence());
            this.fastaBufferedWriter.write("\n");
            this.fileLocker.unlock();
        }else{
            throw new IOException("Fasta Buffer closed");
        }
    }
    
    
    /**
     * closes this bufferedWriter
     * @throws IOException 
     */
    public void close() throws IOException{
        this.fileLocker.waitTillCompleteUnlock();
        this.fastaBufferedWriter.close();
    }
    
    
    
    /**
     * closes first the file
     * then execute finilize()
     * @throws Throwable
     */
    @Override
    public void finalize() throws Throwable{
        try {
            this.fastaBufferedWriter.close();
        } catch (IOException ex) {
            Logger.getLogger(FastaBufferedWriter.class.getName()).log(Level.SEVERE, null, ex);
        }finally{
            super.finalize();
        }
    }
}
