/*
 * This is GBSX v1.0. A toolkit for experimental design and demultiplexing genotyping by sequencing experiments. 
 *  
 * Copyright (C) 2014 KU Leuven
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
package be.uzleuven.gc.logistics.GBSX.demultiplexer;

import be.uzleuven.gc.logistics.GBSX.demultiplexer.application.FastqDemultiplex;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions.ErrorInLogException;
import be.uzleuven.gc.logistics.GBSX.utils.exceptions.StopExcecutionException;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.infrastructure.fileInteractors.LoggerFile;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class GBSdemultiplex {
    
    public static final boolean DEBUG = true;
    public final static String VERSION = "GBSX v1.0";
    public final static String LICENCE = "GPLv3";
    
    /*
     * 
     * 
     */
    public static void main (String[] args){
//        if (GBSdemultiplex.DEBUG){
//            String input = "-help";
//            args = input.split(" ");
//        }
        GBSdemultiplex gbsdemultiplex = new GBSdemultiplex();
        gbsdemultiplex.runDemultiplex(args);
    }
            
    
    /**
     * The demultiplex
     * @param args 
     */
    public void runDemultiplex(String[] args){
        try{
            FastqDemultiplex demultiplex = new FastqDemultiplex(args);
            System.out.println("Start the demultiplexing.");
            demultiplex.processFastqFiles();
            System.out.println("Demultiplexing ended.");
        } catch (StopExcecutionException ex){
            //End excecution as help of version is given
            if (GBSdemultiplex.DEBUG){
                ex.printStackTrace();
            }
        } catch (Exception e){
            if (GBSdemultiplex.DEBUG){
                e.printStackTrace();
            }
            try {
                LoggerFile loggerFile = new LoggerFile(System.getProperty("user.dir"));
                loggerFile.addToLog("Arguments where: " + args.toString());
                loggerFile.addToLog("\n \n \n \n \n");
                loggerFile.addToLog(e.getMessage());
                loggerFile.closeLog();
            } catch (ErrorInLogException ex) {
                ex.printStackTrace();
                Logger.getLogger(GBSdemultiplex.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                ex.printStackTrace();
                Logger.getLogger(GBSdemultiplex.class.getName()).log(Level.SEVERE, null, ex);
            } finally{
                System.out.println("GBS demultiplexing ended with Major Errors");
            }
        }
    }
    
    
    
    
    
    
    
    
}
