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
package be.uzleuven.gc.logistics.GBSX.barcodeGenerator;

import be.uzleuven.gc.logistics.GBSX.GBSX;
import be.uzleuven.gc.logistics.GBSX.barcodeGenerator.infrastructure.fileInteractors.BarcodeReader;
import be.uzleuven.gc.logistics.GBSX.barcodeGenerator.infrastructure.fileInteractors.BarcodeSummaryWriter;
import be.uzleuven.gc.logistics.GBSX.barcodeGenerator.infrastructure.fileInteractors.BarcodeWriter;
import be.uzleuven.gc.logistics.GBSX.barcodeGenerator.model.BarcodeGeneratorArguments;
import be.uzleuven.gc.logistics.GBSX.barcodeGenerator.model.BarcodeGeneratorParameters;
import be.uzleuven.gc.logistics.GBSX.utils.exceptions.StopExcecutionException;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.infrastructure.EnzymeFileParser;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeGenerator {
    
    public static final boolean DEBUG = false;
    public final static String VERSION = "Barcode Generator v1.0";
    public final static String LICENCE = "GPLv3";
    
    private Random random = new Random();
    private BarcodeGeneratorParameters parameters;
    
    public static void main (String[] args){
        try{
            if (BarcodeGenerator.DEBUG){
                String input = "-help";
                args = input.split(" ");
            }
            BarcodeGenerator barcodeGenerator = new BarcodeGenerator(args);
            barcodeGenerator.runBarcodeGenerator();
        } catch (StopExcecutionException ex){
            //End excecution as help of version is given
            if (BarcodeGenerator.DEBUG){
                ex.printStackTrace();
            }
        }
    }
    
    public BarcodeGenerator(String[] args) throws StopExcecutionException{
        this.parseParameters(args);
    }
    
    private void parseParameters(String[] args) throws StopExcecutionException{
        this.parameters = new BarcodeGeneratorParameters();
        if (args.length == 0) {
            throw new IllegalArgumentException("\n\nNo arguments given.\n\n");
        }
        if (args[0].equals("version") || args[0].equals("-version") || args[0].equals("-v")){
            System.out.println(BarcodeGenerator.VERSION);
            System.out.println("This is " + GBSX.VERSION + ". A toolkit for experimental design and demultiplexing genotyping by \n" +
"sequencing experiments.");
            throw new StopExcecutionException();
        }
        if (args[0].equals("help") || args[0].equals("-help") || args[0].equals("-h")){
            BarcodeGenerator.getHelp();
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
                    BarcodeGeneratorArguments argument = BarcodeGeneratorArguments.INVALID_ARGUMENT.getArgument(arg);
                    this.parameters.setParameter(argument, args[++i]);
                }
            }
        }
        //parse possible enzyme files
        if (this.parameters.mustAddEnzymes()){
            String fileName = this.parameters.getEnzymeFile();
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
                Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Enzyme file.", ex);
            } catch (IOException ex) {
                Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Enzyme file", ex);
            }
        }
        //parse possible barcode files
        if (this.parameters.hasBasicSet()){
            try {
                BarcodeReader barcodeReader = new BarcodeReader(new File(this.parameters.getBasicSetFile()));
                this.parameters.getBasicSet().addAll(barcodeReader.readAll());
            } catch (FileNotFoundException ex) {
                Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Basic Barcode file.", ex);
            } catch (IOException ex) {
                Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Basic Barcode file.", ex);
            }
            
        }
        if (this.parameters.hasNotUseSet()){
            try {
                BarcodeReader barcodeReader = new BarcodeReader(new File(this.parameters.getNotUseSetFile()));
                this.parameters.getNotUseSet().addAll(barcodeReader.readAll());
            } catch (FileNotFoundException ex) {
                Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Not Use Barcode file.", ex);
            } catch (IOException ex) {
                Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
                throw new RuntimeException("Couldn't open the Not Use Barcode file.", ex);
            }
            this.parameters.getBasicSet().removeAll(this.parameters.getNotUseSet());
        }
        //check all parameters
        if (! this.parameters.areRequiredParametersSet()){
            throw new RuntimeException(this.parameters.getErrorRequiredParametersSet());
        }
        File outDir = new File(this.parameters.getOutputDirectory());
        if (! outDir.exists()){
            outDir.mkdirs();
        }
    }
    
    /**
     * print help to standard output
     */
    public static void getHelp(){
        System.out.println("");
        System.out.println("");
        System.out.println("This is the help of the BarcodeGenerator.");
        System.out.println();
        System.out.println("KU Leuven");
        System.out.println("Licenced under " + BarcodeGenerator.LICENCE);
        System.out.println(BarcodeGenerator.VERSION);
        System.out.println();
        System.out.println("This program generates a given number of random Barcodes.");
        System.out.println();
        System.out.println(new BarcodeGeneratorParameters().getParametersHelp());
        System.out.println();
        System.out.println();
        System.out.println("Developed by the KU Leuven 2014");
        System.out.println("Licenced under " + BarcodeGenerator.LICENCE);
        System.out.println("For licence information use -licence");
        System.out.println();
    }
    
    public void runBarcodeGenerator(){
        this.createBarcodes(this.parameters.getEnzyme(), this.parameters.getNumberOfBarcodes());
    }
    
    /**
     * creates the barcodes and writes them to a file
     * @param enzyme Enzyme | the enzyme used in the experiment
     * @param numberOfBarcodes int | number of barcodes that must be generated
     * @return collection of String | all generated barcodes
     */
    public Collection<String> createBarcodes(Enzyme enzyme, int numberOfBarcodes){
        //the possible bases
        HashMap<Character, Integer> basesIntMap = new HashMap();
        basesIntMap.put('A', 0);
        basesIntMap.put('C', 1);
        basesIntMap.put('G', 2);
        basesIntMap.put('T', 3);
        //set for all found barcodes
        HashSet<String> bestBarcodesTry = new HashSet();
        double[][] bestPositionsCount = new double[15 + this.parameters.getBigestEnzymeLength()][4];
        int currentBootstrap = 0;
        boolean mustTryAgain = true;
        while (currentBootstrap < this.parameters.getNumberOfBootstraps() && mustTryAgain){
            //multiple tries counter when working with complex enzymes, or lot of barcodes
            currentBootstrap++;
            HashSet<String> barcodes = new HashSet();
            //count of all positions
            double[][] positionsCount = new double[15 + this.parameters.getBigestEnzymeLength()][4];
            for (int pos1 = 0; pos1 < positionsCount.length; pos1++){
                for (int pos2 = 0; pos2 < 4; pos2++){
                    positionsCount[pos1][pos2] = 0.0;
                }
            }
            //go over every possible enzyme position and add there values to the count
            for (int number = 0; number < numberOfBarcodes; number++){
                int pos = 8 + (number % 8);
                for (String cutsite : enzyme.getInitialCutSiteRemnant()){
                    for (int pos1 = 0; pos1 < cutsite.length() && (pos + pos1) < positionsCount.length; pos1++){
                        positionsCount[pos + pos1][basesIntMap.get(cutsite.charAt(pos1))] += (1 / (enzyme.getInitialCutSiteRemnant().size() * 1.0));
                    }
                }
            }
            //go over every already used barcode and enzyme
            if (this.parameters.hasBasicSet()){
                for (String barcode : this.parameters.getBasicSet()){
                    int i;
                    for (i = 0; i < barcode.length(); i++){
                        positionsCount[i][basesIntMap.get(barcode.charAt(i))] += 1;
                    }
                    for (String cutsite : enzyme.getInitialCutSiteRemnant()){
                        for (int ii = 0; ii < cutsite.length() && ii < positionsCount.length; ii++){
                            positionsCount[i + ii][basesIntMap.get(cutsite.charAt(ii))] += (1 / (enzyme.getInitialCutSiteRemnant().size() * 1.0));
                        }
                    }
                }
            }
            //create a list of every possible base on that position
            ArrayList<ArrayList<Character>> positionBaseSet = new ArrayList();
            for (int pos1 = 0; pos1 < 15 + this.parameters.getBigestEnzymeLength(); pos1++){
                ArrayList<Character> baseset = new ArrayList();
                baseset.add('A');
                baseset.add('C');
                baseset.add('G');
                baseset.add('T');
                positionBaseSet.add(baseset);
            }
            //create all barcodes
            int currentNumberOfBarcodes = 0;
            //multiple trys
            int currentTrys = 0;
            int totalNumberOfBarcodes = numberOfBarcodes + this.parameters.getBasicSet().size();
            while (currentNumberOfBarcodes < numberOfBarcodes && currentTrys < this.parameters.getNumberOfBarcodeTries()){
                //first update the list
                for (int i = 0; i < 15; i++){
                    for (char base : basesIntMap.keySet()){
                        if (i < 8){
                            if (4 - positionBaseSet.get(i).size() < totalNumberOfBarcodes % 4){
                                if (positionsCount[i][basesIntMap.get(base)] >= Math.ceil(totalNumberOfBarcodes / 4.0)){
                                    if (positionBaseSet.get(i).contains(base)){
                                        positionBaseSet.get(i).remove(positionBaseSet.get(i).indexOf(base));
                                    }
                                }
                            }else{
                                if (positionsCount[i][basesIntMap.get(base)] >= Math.floor(totalNumberOfBarcodes / 4.0)){
                                    if (positionBaseSet.get(i).contains(base)){
                                        positionBaseSet.get(i).remove(positionBaseSet.get(i).indexOf(base));
                                    }
                                }
                            }
                        }else{
                            if (Math.floor(positionsCount[i][basesIntMap.get(base)]) >= Math.ceil(totalNumberOfBarcodes / 4.0)){
                                if (positionBaseSet.get(i).contains(base)){
                                    positionBaseSet.get(i).remove(positionBaseSet.get(i).indexOf(base));
                                }
                            }
                        }
                    }
                }
                //create the random barcode
                char d1 = positionBaseSet.get(2).get(this.random.nextInt(positionBaseSet.get(2).size()));
                char d2 = positionBaseSet.get(4).get(this.random.nextInt(positionBaseSet.get(4).size()));
                char d3 = positionBaseSet.get(5).get(this.random.nextInt(positionBaseSet.get(5).size()));
                char d4 = positionBaseSet.get(6).get(this.random.nextInt(positionBaseSet.get(6).size()));
                char d5 = positionBaseSet.get(8).get(this.random.nextInt(positionBaseSet.get(8).size()));
                char d6 = positionBaseSet.get(9).get(this.random.nextInt(positionBaseSet.get(9).size()));
                char d7 = positionBaseSet.get(10).get(this.random.nextInt(positionBaseSet.get(10).size()));
                char d8 = positionBaseSet.get(11).get(this.random.nextInt(positionBaseSet.get(11).size()));
                char d9 = positionBaseSet.get(12).get(this.random.nextInt(positionBaseSet.get(12).size()));
                char d10 = positionBaseSet.get(13).get(this.random.nextInt(positionBaseSet.get(13).size()));
                char d11 = positionBaseSet.get(14).get(this.random.nextInt(positionBaseSet.get(14).size()));
                char[] bases = new char[11];
                bases[0] = d1;
                bases[1] = d2;
                bases[2] = d3;
                bases[3] = d4;
                bases[4] = d5;
                bases[5] = d6;
                bases[6] = d7;
                bases[7] = d8;
                bases[8] = d9;
                bases[9] = d10;
                bases[10] = d11;
                int currentTotalNumberOfBarcodes = currentNumberOfBarcodes + this.parameters.getBasicSet().size();
                for (int i = 4 + (currentTotalNumberOfBarcodes % 8); i < bases.length; i++){
                    bases[i] = 'A';
                }
                String barcode = this.generateBarcode(bases).substring(0, 8 + (currentTotalNumberOfBarcodes % 8));
                //check if the barcode is ok
                if (this.isGoodBarcode(barcode, enzyme, barcodes)){
                    //barcode is ok: add to the list and update the positions count
                    char[] barcodeBases = barcode.toCharArray();
                    for (int i = 0; i < 8 + (currentTotalNumberOfBarcodes % 8); i++){
                        positionsCount[i][basesIntMap.get(barcodeBases[i])] += 1;
                    }
                    barcodes.add(barcode);
                    currentNumberOfBarcodes++;
                    currentTrys = 0;
                }else{
                    //try a new barcode
                    currentTrys++;
                }
            }
            //check if beter barcode: more barcodes?
            if (barcodes.size() > bestBarcodesTry.size()){
                bestBarcodesTry = new HashSet(barcodes);
                for (int i = 0; i < positionsCount.length; i++){
                    for (int ii = 0; ii < 4; ii++){
                        bestPositionsCount[i][ii] = positionsCount[i][ii];
                    }
                }
            }
            //check if beter barcodes: beter weights for the barcode bases distribution
            if (barcodes.size() == bestBarcodesTry.size()){
                double[] newScores = new double[8 + this.parameters.getSmallestEnzymeLength()];
                double[] oldScores = new double[8 + this.parameters.getSmallestEnzymeLength()];
                for (int i = 0; i < 8 + this.parameters.getSmallestEnzymeLength(); i++){
                    newScores[i] = this.getScore(positionsCount[i]);
                    oldScores[i] = this.getScore(bestPositionsCount[i]);
                }
                boolean replace = false;
                for (int i = 0; i < newScores.length && ! replace; i++){
                    double[] newnewScores = new double[newScores.length - i];
                    double[] oldoldScores = new double[oldScores.length - i];
                    for (int ii = 0; ii < newnewScores.length; ii++){
                        newnewScores[ii] = newScores[ii];
                        oldoldScores[ii] = oldScores[ii];
                    }
                    if (this.getScore(newnewScores) > this.getScore(oldoldScores)){
                        replace = true;
                    }
                }
                if (replace){
                    //save the new barcode
                    bestBarcodesTry = new HashSet(barcodes);
                    for (int i = 0; i < positionsCount.length; i++){
                        for (int ii = 0; ii < 4; ii++){
                            bestPositionsCount[i][ii] = positionsCount[i][ii];
                        }
                    }
                    boolean isUltimateScore = true;
                    for (double i : newScores){
                        double score = i * (Math.ceil(this.parameters.getNumberOfBarcodes() / 4.0) / Math.floor(this.parameters.getNumberOfBarcodes() / 4.0));
                        if (score < 0.5){
                            //0.5 is the score for each base a minimum of 12.5% occurance
                            isUltimateScore = false;
                        }
                    }
                    if (isUltimateScore){
                        mustTryAgain = false;
                    }
                }
            }
            //check if the condition can be stopped
            if (bestBarcodesTry.size() == numberOfBarcodes && !this.parameters.findUltimeBarcodeCombination()){
                mustTryAgain = false;
            }
            //End of bootstrap
        }
        
        
        //calculation of the positionQuality
        double[] quality = new double[15 + this.parameters.getBigestEnzymeLength()];
        for (int i = 0; i < 8 + this.parameters.getSmallestEnzymeLength(); i++){
            quality[i] = Math.min(Math.min(bestPositionsCount[i][0], bestPositionsCount[i][1]), Math.min(bestPositionsCount[i][2], bestPositionsCount[i][3])) / Math.max(Math.max(bestPositionsCount[i][0], bestPositionsCount[i][1]), Math.max(bestPositionsCount[i][2], bestPositionsCount[i][3]));
            quality[i] = quality[i] * (Math.ceil(this.parameters.getNumberOfBarcodes() / 4.0) / Math.floor(this.parameters.getNumberOfBarcodes() / 4.0));
        }
        for (int i = 8 + this.parameters.getSmallestEnzymeLength(); i < 15; i++){
            quality[i] = Double.NaN;
        }
        
        try {
            HashSet<String> allBarcodes = new HashSet(bestBarcodesTry);
            allBarcodes.addAll(this.parameters.getBasicSet());
            BarcodeWriter barcodeWriter = new BarcodeWriter(this.parameters.getOutputDirectory());
            barcodeWriter.write(allBarcodes);
            BarcodeSummaryWriter barcodeSummaryWriter = new BarcodeSummaryWriter(this.parameters);
            barcodeSummaryWriter.write(bestBarcodesTry, bestPositionsCount, basesIntMap, quality, currentBootstrap);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(BarcodeGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        System.out.println("");
        System.out.println("A\tC\tG\tT");
        for (int i = 0; i < bestPositionsCount.length; i++){
            for (int ii = 0; ii < 4; ii++){
                System.out.print(bestPositionsCount[i][ii] + "\t");
            }
            System.out.println("");
        }
        System.out.println("");
        System.out.println("");
        for (String barcode : bestBarcodesTry){
            System.out.println(barcode);
        }
        System.out.println("");
        System.out.println("");
        System.out.println("");
        System.out.println("");
        System.out.println("Number of barcodes asked: " + numberOfBarcodes + "\t Number of barcodes calculated: " + bestBarcodesTry.size());
        System.out.println("Used " + currentBootstrap + " of " + this.parameters.getNumberOfBootstraps() + " bootstraps.");
        
        return bestBarcodesTry;
    }
    
    /**
     * 
     * @param barcode String | the barcode
     * @param enzyme Enzyme | the enzyme used in the experiment
     * @param barcodes collection of String | all barcodes found until now
     * @return true if it's a good barcode, false otherwise
     */
    private boolean isGoodBarcode(String barcode, Enzyme enzyme, Collection<String> barcodes){
        if (barcode.contains("AAA") || barcode.contains("CCC") || barcode.contains("GGG") || barcode.contains("TTT")){
            return false;
        }
        for (String cutsites : enzyme.getInitialCutSiteRemnant()){
            if (barcode.contains(cutsites)){
                if (! cutsites.equals("")){
                    return false;
                }
            }
        }
        for (String otherBarcodes : barcodes){
            if (barcode.contains(otherBarcodes) || otherBarcodes.contains(barcode)){
                return false;
            }
        }
        if (this.parameters.getNotUseSet().contains(barcode)){
            return false;
        }
        return true;
    }
    
    /**
     * 
     * @param bases array of 11 char
     * @return String | the new barcode
     * @see BarcodeGenerator#generateBarcode(char, char, char, char, char, char, char, char, char, char, char) 
     */
    private String generateBarcode(char[] bases){
        return this.generateBarcode(bases[0], bases[1], bases[2], bases[3], bases[4], bases[5], bases[6], bases[7], bases[8], bases[9], bases[10]);
    }
    
    /**
     * Generates barcodes using the Hammings(14,11) principe
     * @param d1 char data base 1
     * @param d2 char data base 2
     * @param d3 char data base 3
     * @param d4 char data base 4
     * @param d5 char data base 5
     * @param d6 char data base 6
     * @param d7 char data base 7
     * @param d8 char data base 8
     * @param d9 char data base 9
     * @param d10 char data base 10
     * @param d11 char data base 11
     * @return String | the new barcode
     */
    private String generateBarcode(char d1, char d2, char d3, char d4, char d5, char d6, char d7, char d8, char d9, char d10, char d11){
        HashMap<Character, Integer> basesIntMap = new HashMap();
        basesIntMap.put('A', 0);
        basesIntMap.put('C', 1);
        basesIntMap.put('G', 2);
        basesIntMap.put('T', 3);
        int p1 = (4 - (basesIntMap.get(d1) + basesIntMap.get(d2) + basesIntMap.get(d4) + basesIntMap.get(d5) 
                + basesIntMap.get(d7) + basesIntMap.get(d9) + basesIntMap.get(d11)) % 4) % 4;
        int p2 = (4 - (basesIntMap.get(d1) + basesIntMap.get(d3) + basesIntMap.get(d4) + basesIntMap.get(d6) 
                + basesIntMap.get(d7) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        int p3 = (4 - (basesIntMap.get(d2) + basesIntMap.get(d3) + basesIntMap.get(d4) + basesIntMap.get(d8) 
                + basesIntMap.get(d9) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        int p4 = (4 - (basesIntMap.get(d5) + basesIntMap.get(d6) + basesIntMap.get(d7) + basesIntMap.get(d8) 
                + basesIntMap.get(d9) + basesIntMap.get(d10) + basesIntMap.get(d11)) % 4) % 4;
        char base1 = 'N';
        char base2 = 'N';
        char base3 = 'N';
        char base4 = 'N';
        for (char base : basesIntMap.keySet()){
            if (basesIntMap.get(base) == p1) base1 = base;
            if (basesIntMap.get(base) == p2) base2 = base;
            if (basesIntMap.get(base) == p3) base3 = base;
            if (basesIntMap.get(base) == p4) base4 = base;
        }
        String barcode = Character.toString(base1) + base2 + d1 + base3 + d2 + d3 + d4 + base4 + d5 + d6 + d7 + d8 + d9 + d10 + d11;
        return barcode;
    }
    
    /**
     * the score of this array is calculated
     * the score is the minimum diveded by the maximum
     * @param baseCounts array of doubles
     * @return int | the score of the array
     */
    private double getScore(double[] baseCounts){
        double minimum = Double.MAX_VALUE;
        double maximum = 0.0;
        for (double i : baseCounts){
            if (i < minimum){
                minimum = i;
            }
            if (i > maximum){
                maximum = i;
            }
        }
        return minimum / maximum;
    }
    
}
