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
package be.uzleuven.gc.logistics.GBSX.barcodeDiscovery;

import be.uzleuven.gc.logistics.GBSX.GBSX;
import be.uzleuven.gc.logistics.GBSX.barcodeDiscovery.model.BarcodeArguments;
import be.uzleuven.gc.logistics.GBSX.barcodeDiscovery.model.BarcodeParameters;
import be.uzleuven.gc.logistics.GBSX.utils.exceptions.StopExcecutionException;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.infrastructure.EnzymeFileParser;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqBufferedReader;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeComparator;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqRead;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeDiscovery {
    
    
    
    /**
     * 
     * @param args 
     */
    public static void main (String[] args){
        try{
            if (BarcodeDiscovery.DEBUG){
//                String input = "-help";
                String input ="-f1 /home/koen/Downloads/GBS_Example/dem/dem2/output.fastq";
                args = input.split(" ");
            }
            BarcodeDiscovery barcodeFinder = new BarcodeDiscovery(args);
            barcodeFinder.runBarcodeFinder();
        } catch (StopExcecutionException ex){
            //End excecution as help of version is given
            if (BarcodeDiscovery.DEBUG){
                ex.printStackTrace();
            }
        }
    }
    
    public static final boolean DEBUG = false;
    private BarcodeParameters parameters;
    public final static String VERSION = "Barcode Discovery v1.0";
    public final static String LICENCE = "GPLv3";
    
    public BarcodeDiscovery(String[] args) throws StopExcecutionException{
        this.parseParameters(args);
    }
    
    private void parseParameters(String[] args) throws StopExcecutionException{
        this.parameters = new BarcodeParameters();
        if (args.length == 0) {
            throw new IllegalArgumentException("\n\nNo arguments given.\n\n");
        }
        if (args[0].equals("version") || args[0].equals("-version") || args[0].equals("-v")){
            System.out.println(BarcodeDiscovery.VERSION);
            System.out.println("This is " + GBSX.VERSION + ". A toolkit for experimental design and demultiplexing genotyping by \n" +
"sequencing experiments.");
            throw new StopExcecutionException();
        }
        if (args[0].equals("help") || args[0].equals("-help") || args[0].equals("-h")){
            BarcodeDiscovery.getHelp();
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
        //parse parameters
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];

            // An option.
            if (arg.startsWith("--") || arg.startsWith("-")) {
                //check if the option is invalid
                if (i + 1 >= args.length || args[i + 1].startsWith("-")) {
                    //invalid option
                    throw new RuntimeException("Value required for option "
                            + arg);
                }else{
                    //get the argument
                    BarcodeArguments argument = BarcodeArguments.INVALID_ARGUMENT.getArgument(arg);
                    this.parameters.setParameter(argument, args[++i]);
                }
            }
        }
        //parse possible enzyme files
        if (this.parameters.mustAddEnzymes() || this.parameters.mustReplaceEnzymes()){
            String fileName;
            if (this.parameters.mustAddEnzymes()){
                fileName = this.parameters.getAddEnzymeFile();
            }else{
                fileName = this.parameters.getReplaceEnzymeFile();
            }
            try {
                EnzymeFileParser enzymeFileParser = new EnzymeFileParser();
                ArrayList<Enzyme> enzymeList = new ArrayList();
                enzymeList.addAll(enzymeFileParser.parseEnzymeFile(fileName));
                if (enzymeList.isEmpty()){
                    throw new RuntimeException("No valid enzymes found in the enzyme file");
                }
                if (this.parameters.mustAddEnzymes()){
                    this.parameters.getEnzymeCollection().addEnzymes(enzymeList);
                }else{
                    this.parameters.getEnzymeCollection().replaceEnzymes(enzymeList);
                }
            } catch (FileNotFoundException ex) {
                Logger.getLogger(BarcodeDiscovery.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Enzyme file.", ex);
            } catch (IOException ex) {
                Logger.getLogger(BarcodeDiscovery.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Enzyme file", ex);
            }
        }
        if (! this.parameters.areRequiredParametersSet()){
            throw new RuntimeException(this.parameters.getErrorRequiredParametersSet());
        }
        File outDir = new File(this.parameters.getOutputDirectory());
        if (! outDir.exists()){
            outDir.mkdirs();
        }
    }
    
    /**
     * prints the help of this barcode finder to the standard output (system.out.println)
     */
    public static void getHelp(){
        System.out.println("");
        System.out.println("");
        System.out.println("This is the help of the BarcodeFinder.");
        System.out.println();
        System.out.println("KU Leuven");
        System.out.println("Licenced under " + BarcodeDiscovery.LICENCE);
        System.out.println(BarcodeDiscovery.VERSION);
        System.out.println();
        System.out.println("This program search for possible barcodes, and barcode enzyme combinations.");
        System.out.println();
        System.out.println(new BarcodeParameters().getParametersHelp());
        System.out.println();
        System.out.println();
        System.out.println("Developed by the KU Leuven 2014");
        System.out.println("Licenced under " + BarcodeDiscovery.LICENCE);
        System.out.println("For licence information use -licence");
        System.out.println();
    }
    
    /**
     * the program for the running of the barcode finder
     */
    public void runBarcodeFinder(){
        int shortestBarcode = this.parameters.getMinimumLength();
        int longestBarcode = this.parameters.getMaximumLength();
        int longestEnzymeSite = 0;
        if (this.parameters.mustSearchEnzyme()){
            longestEnzymeSite = this.getLongestEnzymeLength();
        }
        System.out.println("Start the Barcode finder.");
        //open the files and the reader
        File fastqFile1 = new File(this.parameters.getFileName());
        try {
            //create buffered readers for the files
            FastqBufferedReader fastq1Reader = new FastqBufferedReader(fastqFile1, this.parameters.mustBeZipped());
            
            //create needed safe maps
            HashMap<Integer, HashMap> barcodeMap = new HashMap();
            for (int index = shortestBarcode; index <= (longestBarcode + longestEnzymeSite); index++){
                HashMap<String, Integer> sequenceCountMap = new HashMap();
                barcodeMap.put(index, sequenceCountMap);
            }
            
            //init vars needed in the loop to go furter
            FastqRead fastq1;
            while ((fastq1 = fastq1Reader.next()) != null) {
                //read the next fastq line
                //go over every possible barcode length
                for (int index = shortestBarcode; index <= (longestBarcode + longestEnzymeSite); index++){
                    HashMap<String, Integer> sequenceCountMap = barcodeMap.get(index);
                    //put into the map
                    if (sequenceCountMap.containsKey(fastq1.getSequence().substring(0, index))){
                        sequenceCountMap.put(fastq1.getSequence().substring(0, index), sequenceCountMap.get(fastq1.getSequence().substring(0, index)) + 1);
                    }else{
                        sequenceCountMap.put(fastq1.getSequence().substring(0, index), 1);
                    }
                }
                 
            }
            System.out.println("All sequences run");
            //go over every map
            //and write to file
            for (int index = shortestBarcode; index <= longestBarcode; index++){
                HashMap<String, Integer> sequenceCountMap = barcodeMap.get(index);
            
                System.out.println("File " + index + " begun");
                //write the string to a file
                File file = new File(this.parameters.getOutputDirectory() + System.getProperty("file.separator") +"counts." + index + ".cnt");
                try {
                    file.createNewFile();
                    BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file));
                    for (String sequence : sequenceCountMap.keySet()){
                        bufferedWriter.write(sequence + "\t" + sequenceCountMap.get(sequence) + "\n");
                    }
                    bufferedWriter.close();
                } catch (IOException ex) {
                    if (BarcodeDiscovery.DEBUG){
                        ex.printStackTrace();
                    }
                }
                System.out.println("File " + index + " saved");
            }
            System.out.println("File saved");
            
            int totalReads = 0;
            
            //clean the sort barcodes
            for (int index = shortestBarcode; index <= (longestBarcode + longestEnzymeSite); index++){
                HashMap<String, Integer> sequenceCountMap = barcodeMap.get(index);
                HashMap<String, Integer> cleanCountMap = new HashMap();
                for (String sequence : sequenceCountMap.keySet()){
                    totalReads += sequenceCountMap.get(sequence);
                    if (sequenceCountMap.get(sequence) > this.parameters.getMinimumOccuranceBarcodes()){
                        cleanCountMap.put(sequence, sequenceCountMap.get(sequence));
                    }
                }
                barcodeMap.put(index, cleanCountMap);
            }
            int lowestBarcodeOccurance = Integer.MAX_VALUE;
            HashMap<String, Integer> possibleBarcodes = new HashMap();
            HashMap<String, HashSet<Enzyme>> possibleBarcodeEnzymes = new HashMap();
            
            
            
            
            if (this.parameters.mustSearchEnzyme()){
                //enzyme finder
                EnzymeComparator enzymeComparator = new EnzymeComparator();
                for (int index = shortestBarcode; index <= longestBarcode; index++){
                    HashMap<String, Integer> sequenceCountMap = barcodeMap.get(index);
                    for (String sequence : sequenceCountMap.keySet()){
                        for (Enzyme enzyme : this.parameters.getAllEnzymes()){
                            if (enzymeComparator.compare(enzyme, this.parameters.getNeutralEnzyme()) == 0){
                                
                            }else{
                                for (String cutsites : enzyme.getInitialCutSiteRemnant()){
                                    HashMap<String, Integer> tempMap = barcodeMap.get(sequence.length() + cutsites.length());
                                    if (tempMap.containsKey(sequence + cutsites)){
                                        possibleBarcodes.put(sequence, sequenceCountMap.get(sequence));
                                        if (! possibleBarcodeEnzymes.containsKey(sequence)){
                                            possibleBarcodeEnzymes.put(sequence, new HashSet());
                                        }
                                        possibleBarcodeEnzymes.get(sequence).add(enzyme);
                                    }
                                }
                            }
                        }
                    }
                }
                
            }else{
                //no enzyme finder
                //create the first map of the shortest most occuring samples
                for (int index = shortestBarcode; index <= shortestBarcode; index++){
                    HashMap<String, Integer> sequenceCountMap = barcodeMap.get(index);
                    for (String sequence : sequenceCountMap.keySet()){
                        //test every sequence
                        if (possibleBarcodes.keySet().size() == this.parameters.getMaximumNumbersBarcodes()){
                            //map is full: check to replace
                            if (sequenceCountMap.get(sequence) > lowestBarcodeOccurance){
                                //replace lowest with new
                                int newLowest = Integer.MAX_VALUE;
                                String toRemove = "";
                                for (String barcodes : possibleBarcodes.keySet()){
                                    //find barcode with lowest occurance
                                    if (possibleBarcodes.get(barcodes) == lowestBarcodeOccurance){
                                        toRemove = barcodes;
                                    }else{
                                        //find new lowest occurance
                                        if (possibleBarcodes.get(barcodes) < newLowest){
                                            newLowest = possibleBarcodes.get(barcodes);
                                        }
                                    }
                                }
                                //check new lowest occurance, to last new
                                if (sequenceCountMap.get(sequence) < newLowest){
                                    newLowest = sequenceCountMap.get(sequence);
                                }
                                lowestBarcodeOccurance = newLowest;
                                //replace old with new
                                if (possibleBarcodes.containsKey(toRemove)){
                                    possibleBarcodes.remove(toRemove);
                                }
                                possibleBarcodes.put(sequence, sequenceCountMap.get(sequence));
                            }
                        }else{
                            //not yet full: fill map
                            possibleBarcodes.put(sequence, sequenceCountMap.get(sequence));
                            if (sequenceCountMap.get(sequence) < lowestBarcodeOccurance){
                                lowestBarcodeOccurance = sequenceCountMap.get(sequence);
                            }
                        }
                    }
                }
                //check longer barcodes: if longer barcode has a shorter in the found as substring, check if the occurance is more than x percent of the previous one
                for (int index = shortestBarcode + 1; index <= longestBarcode; index++){
                    HashMap<String, Integer> sequenceCountMap = barcodeMap.get(index);
                    for (String sequence : sequenceCountMap.keySet()){
                        String toRemove = "";
                        for (String barcodes : possibleBarcodes.keySet()){
                            if (sequence.substring(0, barcodes.length()).equals(barcodes)){
                                //possible new found, check if new one is less then x percent of the previous one (mismatches delete, but SNPs must remain)
                                if (possibleBarcodes.get(barcodes) - sequenceCountMap.get(sequence) < possibleBarcodes.get(barcodes) / 100.0 * this.parameters.getBarcodeMismatchPercentage()){
                                    toRemove = barcodes;
                                }
                            }
                        }
                        //if a new barcode is found, replace
                        if (possibleBarcodes.containsKey(toRemove)){
                            possibleBarcodes.remove(toRemove);
                            possibleBarcodes.put(sequence, sequenceCountMap.get(sequence));
                        }
                    }
                    
                }
                
                
            }
            
            File file = new File(this.parameters.getOutputDirectory() + System.getProperty("file.separator") +"possibleBarcodes.cnt");
            try {
                file.createNewFile();
                BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file));
                for (String sequence : possibleBarcodes.keySet()){
                    if (this.parameters.mustSearchEnzyme()){
                        for (Enzyme e : possibleBarcodeEnzymes.get(sequence)){
                            bufferedWriter.write(sequence + "\t" + possibleBarcodes.get(sequence));
                            bufferedWriter.write("\t" + e.getName());
                            bufferedWriter.write("\n");
                        }
                    }else{
                        bufferedWriter.write(sequence + "\t" + possibleBarcodes.get(sequence));
                        bufferedWriter.write("\n");
                    }
                }
                bufferedWriter.close();
            } catch (IOException ex) {
                if (BarcodeDiscovery.DEBUG){
                    ex.printStackTrace();
                }
            }
            System.out.println("File possible barcodes " + possibleBarcodes.size() + " saved");
            
            BufferedWriter logWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(new File(this.parameters.getOutputDirectory() + System.getProperty("file.separator") + "barcodeDiscovery.log")))));
            logWriter.write(this.parameters.getParametersLogString());
            logWriter.close();
        }catch (Exception e){
            if (BarcodeDiscovery.DEBUG){
                e.printStackTrace();
            }
        }
    }
    
    
    /**
     * Go over every enzyme and determine the length of the longest enzyme cutsite
     * @return 
     */
    private int getLongestEnzymeLength(){
        int longest = 0;
        for (Enzyme enzyme : this.parameters.getAllEnzymes()){
            for (String cutsite : enzyme.getInitialCutSiteRemnant()){
                if (cutsite.length() > longest){
                    longest = cutsite.length();
                }
            }
        }
        return longest;
    }
    
}
