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
package be.uzleuven.gc.logistics.GBSX.demultiplexer.infrastructure.fileInteractors;

import be.uzleuven.gc.logistics.GBSX.GBSX;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.GBSdemultiplex;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.ErrorInLogException;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.model.DemultiplexParameters;
import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Parameters;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Date;
import java.util.concurrent.locks.ReentrantLock;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class LoggerFile {
    
    private Writer logwriter;
    private ReentrantLock lock = new ReentrantLock();
    
    /**
     * creates and opens a logfile at the given location
     * @param outputFile String | the location of the logfile (directory)
     * @throws IOException | when a error occures
     */
    public LoggerFile(String outputFile) throws IOException{
        File file = new File(outputFile + System.getProperty("file.separator") + "gbsDemultiplex.log");
        file.createNewFile();
        this.logwriter = new BufferedWriter(new FileWriter(file));
    }
    
    public LoggerFile(DemultiplexParameters parameters) throws IOException, ErrorInLogException{
        File file = new File(parameters.getOutputDirectory() + System.getProperty("file.separator") + "gbsDemultiplex.log");
        file.createNewFile();
        this.logwriter = new BufferedWriter(new FileWriter(file));
        Date now = new Date();
        this.addToLog("Start GBSX demultiplex on " + now.toString());
        this.addToLog("\n");
        this.addToLog("Toolkit Version: " + GBSX.VERSION + "\t GBS demultiplexer Version: " + GBSdemultiplex.VERSION);
        this.addToLog("\n");
        this.addToLog(parameters.getParametersLogString());
    }
    
    public LoggerFile(Parameters parameters) throws IOException, ErrorInLogException{
        File file = new File(parameters.getOutputDirectory() + System.getProperty("file.separator") + "GBSX.log");
        file.createNewFile();
        this.logwriter = new BufferedWriter(new FileWriter(file));
        Date now = new Date();
        this.addToLog("Start tool on " + now.toString());
        this.addToLog("\n");
        this.addToLog("Toolkit Version: " + GBSX.VERSION);
        this.addToLog("\n");
        this.addToLog(parameters.getParametersLogString());
    }
    
    /**
     * add the given log to the logfile
     * @param log String | the log to be writen to the logfile
     * @throws ErrorInLogException 
     */
    public void addToLog(String log) throws ErrorInLogException{
        try {
            lock.lock();
            this.logwriter.write(log);
        } catch (IOException ex) {
            Logger.getLogger(LoggerFile.class.getName()).log(Level.SEVERE, null, ex);
            throw new ErrorInLogException(log, ex);
        }finally{
            lock.unlock();
        }
    }
    
    /**
     * close the logfile
     * @throws ErrorInLogException 
     */
    public void closeLog() throws ErrorInLogException{
        try {
            lock.lock();
            Date now2 = new Date();
            this.addToLog("Ended on " + now2.toString());
        } catch (ErrorInLogException e){
            
        }finally{
            lock.unlock();
            try {
                this.logwriter.close();
            } catch (IOException ex) {
                Logger.getLogger(LoggerFile.class.getName()).log(Level.SEVERE, null, ex);
                throw new ErrorInLogException("closed failed", ex);
            }
        }
    }
    
    @Override
    public void finalize() throws Throwable{
        this.logwriter.close();
        super.finalize();
    }
    
}
