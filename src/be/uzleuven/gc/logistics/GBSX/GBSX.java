/* 
 * This is GBSX v1.0. A toolkit for experimental design and demultiplexing genotyping by sequencing experiments.
 * 
 * Copyright 2014 KU Leuven
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
package be.uzleuven.gc.logistics.GBSX;

import be.uzleuven.gc.logistics.GBSX.demultiplexer.DNAComplement;
import be.uzleuven.gc.logistics.GBSX.barcodeDiscovery.BarcodeDiscovery;
import be.uzleuven.gc.logistics.GBSX.barcodeGenerator.BarcodeGenerator;
import be.uzleuven.gc.logistics.GBSX.demultiplexer.GBSdemultiplex;
import be.uzleuven.gc.logistics.GBSX.simulator.GBSsimulate;
import be.uzleuven.gc.logistics.GBSX.utils.exceptions.StopExcecutionException;

/**
 * Publication version
 * @author Koen Herten for the KU Leuven
 */
public class GBSX {
    
    public static final boolean DEBUG = false;
    public final static String VERSION = "GBSX v1.1.1";
    private final static String LICENCE = "GPLv3";
    
    /**
     * @param args 
     */
    public static void main (String[] args){
        if (GBSX.DEBUG){
            String input;
            input = "--DEBUG";
            args = input.split(" ");
        }
        try{
            if (GBSX.DEBUG){
                System.out.println("DEBUG MODUS");
            }
            if (args.length == 0) {
                System.out.println("\n\nNo arguments given.\n\n");
                throw new StopExcecutionException();
            }
            if (args[0].equals("version") || args[0].equals("-version") || args[0].equals("-v")){
                System.out.println(GBSX.VERSION);
                System.out.println("This is " + GBSX.VERSION + ". A toolkit for experimental design and demultiplexing genotyping by \n" +
"sequencing experiments.");
                throw new StopExcecutionException();
            }
            if (args[0].equals("-w")){
                System.out.println("GBSX is distributed in the hope that it will be useful,");
                System.out.println("but WITHOUT ANY WARRANTY; without even the implied warranty of");
                System.out.println("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the");
                System.out.println("GNU General Public License for more details.");
                throw new StopExcecutionException();
            }
            if (args[0].equals("-c")){
                System.out.println("GBSX is free software: you can redistribute it and/or modify");
                System.out.println("it under the terms of the GNU General Public License as published by");
                System.out.println("the Free Software Foundation, either version 3 of the License, or");
                System.out.println("(at your option) any later version.");
                throw new StopExcecutionException();
            }
            if (args[0].equals("help") || args[0].equals("-help") || args[0].equals("-h")){
                System.out.println("GBSX v1.0  Copyright (C) 2014 Genomics Core\n" +
                    "This program comes with ABSOLUTELY NO WARRANTY; for details type `perl GBSX_digest_v1.0.pl -w'.\n" +
                    "This is free software, and you are welcome to redistribute it under certain conditions; type `perl GBSX_digest_v1.0.pl -c' for details.\n" +
                    "\n" +
                    "####################################################################################");
                System.out.println("");
                System.out.println("GBSX is a collection of the following programs:");
                System.out.println("\t Demultiplexer  \t\t\t --Demultiplexer \n\t\t the inline barcode, gbs and rad demultiplexer \n \t\t " + GBSdemultiplex.VERSION);
                System.out.println("\t BarcodeDiscovery \t\t --BarcodeDiscovery \n\t\t a quick search for possible barcodes in a fastq file \n\t\t " + BarcodeDiscovery.VERSION);
                System.out.println("\t DNAComplement \t\t --DNAComplement \n\t\t a tool to make the complement of a given DNA sequence \n\t\t " + DNAComplement.VERSION);
                System.out.println("\t BarcodeGenerator \t\t --BarcodeGenerator \n\t\t a tool to generate random self-correcting barcodes (the barcodes has a hammings distance of at least 3)\n\t\t " + BarcodeGenerator.VERSION);
                System.out.println("\t GBSsimulator \t\t --GBSsimulator \n\t\t a tool to simulate GBS data (only purpose is testing)\n\t\t " + BarcodeGenerator.VERSION);
                System.out.println();
                System.out.println("For tool specific help use --tool -help");
                System.out.println();
                System.out.println(GBSX.VERSION);
                System.out.println("Developed by the KU Leuven 2014");
                System.out.println("Licenced under " + GBSX.LICENCE);
                System.out.println("For licence information use -licence");
                System.out.println();
                throw new StopExcecutionException();
            }
            if (args[0].equals("licence") || args[0].equals("-licence") || args[0].equals("-l")){
                System.out.println("This file is part of GBSX.\n" +
                    "\n" +
                    "GBSX is free software: you can redistribute it and/or modify\n" +
                    "it under the terms of the GNU General Public License as published by\n" +
                    "the Free Software Foundation, either version 3 of the License, or\n" +
                    "(at your option) any later version.\n" +
                    "\n" +
                    "GBSX is distributed in the hope that it will be useful,\n" +
                    "but WITHOUT ANY WARRANTY; without even the implied warranty of\n" +
                    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n" +
                    "GNU General Public License for more details.\n" +
                    "\n" +
                    "You should have received a copy of the GNU General Public License\n" +
                    "along with GBSX.  If not, see <http://www.gnu.org/licenses/>.");
                throw new StopExcecutionException();
            }
            String[] newArgs = new String[args.length -1];
            for (int i = 1; i < args.length; i++){
                newArgs[i-1] = args[i];
            }
            if (args[0].toLowerCase().equals("--Demultiplexer".toLowerCase())){
                GBSdemultiplex.main(newArgs);
            }else if (args[0].toLowerCase().equals("--BarcodeDiscovery".toLowerCase())){
                BarcodeDiscovery.main(newArgs);
            }else if (args[0].toLowerCase().equals("--DNAComplement".toLowerCase())){
                DNAComplement.main(newArgs);
            }else if (args[0].toLowerCase().equals("--BarcodeGenerator".toLowerCase())){
                BarcodeGenerator.main(newArgs);
            }else if (args[0].toLowerCase().equals("--GBSsimulator".toLowerCase())){
                GBSsimulate.main(newArgs);
            }else if (args[0].toLowerCase().equals("--DEBUG".toLowerCase())){
                System.out.println("DEBUG is on:");
                System.out.println("Demultiplexer: \t" + GBSdemultiplex.DEBUG);
                System.out.println("BarcodeDiscovery: \t" + BarcodeDiscovery.DEBUG);
                System.out.println("DNAComplement: \t" + DNAComplement.DEBUG);
                System.out.println("BarcodeGenerator: \t" + BarcodeGenerator.DEBUG);
                System.out.println("GBSsimulator: \t" + GBSsimulate.DEBUG);
            }else{
                System.out.println("No valid tool found.");
                System.out.println("Please use -help.");
            }
        } catch (StopExcecutionException ex){
            //End excecution as help of version is given
            if (GBSX.DEBUG){
                ex.printStackTrace();
            }
        } catch (IllegalArgumentException ex){
            //No argument was given
            System.out.println("" + ex.getMessage());
            if (GBSX.DEBUG){
                ex.printStackTrace();
            }
        }
        
        
    }
}
