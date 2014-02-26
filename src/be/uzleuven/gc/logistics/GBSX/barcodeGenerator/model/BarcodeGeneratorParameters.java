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
package be.uzleuven.gc.logistics.GBSX.barcodeGenerator.model;

import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Arguments;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeCollection;
import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Parameters;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeGeneratorParameters implements Parameters{

    private HashMap<BarcodeGeneratorArguments, String> arguments;
    private EnzymeCollection enzymeCollection;
    private HashSet<String> barcodesNotToUse;
    private HashSet<String> barcodesBaseSet;
    
    private final int STANDARD_NUMBER_OF_BOOTSTRAPS = 10000;
    private final int STANDARD_NUMBER_OF_BARCODE_TRIES = 20;
    
    public BarcodeGeneratorParameters(){
        this.arguments = new HashMap();
        //standard parameters
        this.arguments.put(BarcodeGeneratorArguments.OUTPUT_DIRECTORY, System.getProperty("user.dir"));
        this.arguments.put(BarcodeGeneratorArguments.BOOTSTRAPS, "" + this.STANDARD_NUMBER_OF_BOOTSTRAPS);
        this.arguments.put(BarcodeGeneratorArguments.BARCODE_TRIES, "" + this.STANDARD_NUMBER_OF_BARCODE_TRIES);
        this.arguments.put(BarcodeGeneratorArguments.ULTIME_SEARCH_BARCODES, "false");
        this.enzymeCollection = new EnzymeCollection();
        this.barcodesBaseSet = new HashSet();
        this.barcodesNotToUse = new HashSet();
    }
    
    @Override
    public boolean containsParameter(Arguments argument) {
        if (argument instanceof BarcodeGeneratorArguments){
            return this.containsParameter((BarcodeGeneratorArguments) argument);
        }else{
            return false;
        }
    }
    
    public boolean containsParameter(BarcodeGeneratorArguments argument) {
        return this.arguments.containsKey(argument);
    }

    @Override
    public String getParameter(Arguments argument) {
        if (argument instanceof BarcodeGeneratorArguments){
            return this.getParameter((BarcodeGeneratorArguments) argument);
        }else{
            return "ERROR";
        }
    }
    
    public String getParameter(BarcodeGeneratorArguments argument) {
        return this.getParameter((BarcodeGeneratorArguments) argument);
    }

    @Override
    public void setParameter(Arguments argument, String parameter) {
        if (argument instanceof BarcodeGeneratorArguments){
            this.setParameter((BarcodeGeneratorArguments) argument, parameter);
        }
    }
    
    public void setParameter(BarcodeGeneratorArguments argument, String parameter) {
        this.arguments.put(argument, parameter);
    }

    @Override
    public boolean areRequiredParametersSet() {
        boolean isok = true;
        if (! this.arguments.containsKey(BarcodeGeneratorArguments.ENZYME)){
            isok = false;
        }
        if (! this.arguments.containsKey(BarcodeGeneratorArguments.NUMBER_OF_BARCODES)){
            isok = false;
        }
        return isok;
    }

    @Override
    public String getErrorRequiredParametersSet() {
        String errorString = "";
        if (! this.arguments.containsKey(BarcodeGeneratorArguments.ENZYME)){
            errorString += "No enzyme is set.\n";
        }
        if (! this.arguments.containsKey(BarcodeGeneratorArguments.NUMBER_OF_BARCODES)){
            errorString += "No number of barcodes is given.\n";
        }
        return errorString;
    }

    @Override
    public String getParametersLogString() {
        String toLog = "";
        toLog += "Number of barcodes to generate: " + this.getNumberOfBarcodes() + "\n";
        toLog += "Enzyme used in the experiment: " + this.getEnzyme().getName() + "\n";
        if (this.mustAddEnzymes()){
            toLog += "Enzyme file used: " + this.getEnzymeFile() + "\n";
        }
        toLog += "Maximum number of bootstraps: " + this.getNumberOfBootstraps() + "\n";
        toLog += "Maximum number of barcode tries: " + this.getNumberOfBarcodeTries() + "\n";
        toLog += "Must use the ultime barcode search option: " + this.findUltimeBarcodeCombination() + "\n";
        toLog += "Output directory: " + this.getOutputDirectory() + "\n";
        toLog += "Used basic set: " + this.hasBasicSet() + "\n";
        if (this.hasBasicSet()){
            toLog += "Path to the basic set: " + this.getBasicSetFile() + "\n";
        }
        toLog += "Used not use barcode set: " + this.hasNotUseSet() + "\n";
        if (this.hasNotUseSet()){
            toLog += "Path to not use barcode set: " + this.getNotUseSetFile() + "\n";
        }
        return toLog;
    }

    @Override
    public String getParametersHelp() {
        String toHelp = "";
        toHelp += "Mandatory parameters:\n";
        toHelp += "\t -b \t the number of barcodes needed.\n";
        toHelp += "\t -e \t the enzyme used for the experiment.\n";
        toHelp += "Optional parameters:\n";
        toHelp += "\t -ef \t the enzyme file. This option adds new enzymes. The file must be tab delimited: First column the enzyme name, "
                + "second column the cutsites remains (comma separated). \n";
        toHelp += "\t -nb \t the number of bootstraps that maximum must be executed. (standard " + this.STANDARD_NUMBER_OF_BOOTSTRAPS + ")"
                + " By the start of a new bootstrap a complete new design is made. "
                + "The best scored design (most random barcodes and best scored bases distribution is kept as result) \n";
        toHelp += "\t -bt \t the number of barcode tries. (standard " + this.STANDARD_NUMBER_OF_BARCODE_TRIES + ")"
                + " If a random barcode does not fit into the current design try this number of times with a new random barcode before restarting the bootstrap.\n";
        toHelp += "\t -o \t the output directory (standard current working directory) \n";
        toHelp += "\t -us \t try tho find the ultime match: the best barcode combination with the best bases distribution (standard false) "
                + " true: continue even when the right number of barcodes is found.\n";
        toHelp += "\t -bf \t a file with all barcodes that are used as basic set (this file is one of the possible output files)\n";
        toHelp += "\t -nf \t a file with all barcodes that may not be used in the design. If this file contains barcodes that are also found"
                + "in the basic set file, these barcodes will be replaced in the design by new random barcodes\n";
        return toHelp;
    }
    
    
    public int getNumberOfBarcodes(){
        if (this.arguments.containsKey(BarcodeGeneratorArguments.NUMBER_OF_BARCODES)){
            return Integer.parseInt(this.arguments.get(BarcodeGeneratorArguments.NUMBER_OF_BARCODES));
        }else{
            return 0;
        }
    }
    
    /**
     * 
     * @return Enzyme | the enzyme used in the barcode generation
     */
    public Enzyme getEnzyme(){
        if (this.arguments.containsKey(BarcodeGeneratorArguments.ENZYME)){
            return this.enzymeCollection.getEnzyme(this.arguments.get(BarcodeGeneratorArguments.ENZYME));
        }else{
            return this.enzymeCollection.getNeutralEnzyme();
        }
    }
    
    /**
     * 
     * @return int | the length of the smallest enzyme cutsite
     */
    public int getSmallestEnzymeLength(){
        int enzymeLength = 0;
        for (String cutsite : this.getEnzyme().getInitialCutSiteRemnant()){
            if (cutsite.length() < enzymeLength || enzymeLength == 0){
                enzymeLength = cutsite.length();
            }
        }
        return enzymeLength;
    }
    
    /**
     * 
     * @return int  | the length of the largest enzyme cutsite
     */
    public int getBigestEnzymeLength(){
        int enzymeLength = 0;
        for (String cutsite : this.getEnzyme().getInitialCutSiteRemnant()){
            if (cutsite.length() > enzymeLength || enzymeLength == 0){
                enzymeLength = cutsite.length();
            }
        }
        return enzymeLength;
    }
    
    /**
     * 
     * @return EnzymeCollection
     */
    public EnzymeCollection getEnzymeCollection(){
        return this.enzymeCollection;
    }
    
    /**
     * 
     * @return int | the number of bootstraps to run
     */
    public int getNumberOfBootstraps(){
        return Integer.parseInt(this.arguments.get(BarcodeGeneratorArguments.BOOTSTRAPS));
    }
    
    /**
     * 
     * @return int | the number of barcodes that may be created if a barcode failes to be inserted in the collection
     */
    public int getNumberOfBarcodeTries(){
        return Integer.parseInt(this.arguments.get(BarcodeGeneratorArguments.BARCODE_TRIES));
    }
    
    /**
     * 
     * @return String | the output directory
     */
    public String getOutputDirectory(){
        return this.arguments.get(BarcodeGeneratorArguments.OUTPUT_DIRECTORY);
    }
    
    /**
     * 
     * @return true if the user gave an enzyme file
     */
    public boolean mustAddEnzymes(){
        return this.containsParameter(BarcodeGeneratorArguments.ENZYME_FILE);
    }
    
    /**
     * 
     * @return String | the path to the enzyme file
     */
    public String getEnzymeFile(){
        return this.arguments.get(BarcodeGeneratorArguments.ENZYME_FILE);
    }
    
    /**
     * 
     * @return true if the user gave a basic barcode set
     */
    public boolean hasBasicSet(){
        return this.containsParameter(BarcodeGeneratorArguments.BASIC_BARCODES);
    }
    
    /**
     * 
     * @return String | the path to the basic barcode set 
     */
    public String getBasicSetFile(){
        return this.arguments.get(BarcodeGeneratorArguments.BASIC_BARCODES);
    }
    
    /**
     * 
     * @return set of String | the basic barcode set
     */
    public HashSet<String> getBasicSet(){
        return this.barcodesBaseSet;
    }
    
    
    /**
     * 
     * @return true if the user gave a not use or remove barcode set
     */
    public boolean hasNotUseSet(){
        return this.containsParameter(BarcodeGeneratorArguments.NOT_USE_BARCODES);
    }
    
    /**
     * 
     * @return String | the path to the not use barcode set
     */
    public String getNotUseSetFile(){
        return this.arguments.get(BarcodeGeneratorArguments.NOT_USE_BARCODES);
    }
    
    /**
     * 
     * @return set of string | the not use barcode set
     */
    public HashSet<String> getNotUseSet(){
        return this.barcodesNotToUse;
    }
    
    /**
     * 
     * @return true if the ultime barcode combination must be found (calculate min values)
     */
    public boolean findUltimeBarcodeCombination(){
        if (this.arguments.get(BarcodeGeneratorArguments.ULTIME_SEARCH_BARCODES).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
}
