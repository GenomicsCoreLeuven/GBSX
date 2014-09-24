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
package be.uzleuven.gc.logistics.GBSX.barcodeDiscovery.model;

import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Arguments;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeCollection;
import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Parameters;
import java.util.Collection;
import java.util.HashMap;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeParameters implements Parameters{
    
    
    private HashMap<BarcodeArguments, String> arguments;
    private static final int MINIMUM_LENGTH = 6;
    private static final int MAXIMUM_LENGTH = 16;
    private final static int MIN_BARCODE_OCCURANCE = 200;
    private final static int MAX_NUMBER_OF_BARCODES = 100;
    private final static int PERCENTAGE_OF_MISMATCH = 10;
    private EnzymeCollection enzymeCollection;
    
    
    public BarcodeParameters(){
        this.arguments = new HashMap();
        this.arguments.put(BarcodeArguments.MINIMUM_LENGTH, "" + BarcodeParameters.MINIMUM_LENGTH);
        this.arguments.put(BarcodeArguments.MAXIMUM_LENGTH, "" + BarcodeParameters.MAXIMUM_LENGTH);
        this.arguments.put(BarcodeArguments.USE_GZIP_FILES, "false");
        this.arguments.put(BarcodeArguments.OUTPUT_DIRECTORY, System.getProperty("user.dir"));
        this.arguments.put(BarcodeArguments.FIND_ENZYME, "true");
        this.arguments.put(BarcodeArguments.MIN_BARCODE_OCCURANCE, "" + BarcodeParameters.MIN_BARCODE_OCCURANCE);
        this.arguments.put(BarcodeArguments.MAX_NUMBER_OF_BARCODES, "" + BarcodeParameters.MAX_NUMBER_OF_BARCODES);
        this.arguments.put(BarcodeArguments.PERCENTAGE_OF_MISMATCH, "" + BarcodeParameters.PERCENTAGE_OF_MISMATCH);
        this.arguments.put(BarcodeArguments.ENZYME_ADD, "false");
        this.arguments.put(BarcodeArguments.ENZYME_REPLACE, "false");
        this.enzymeCollection = new EnzymeCollection();
    }
    
    /**
     * 
     * @param argument BarcodeArgument
     * @return true if the given argument is a parameter
     */
    public boolean containsParameter(BarcodeArguments argument){
        return this.arguments.containsKey(argument);
    }
    
    /**
     * 
     * @param argument BarcodeArguments
     * @return the asked argument as a String
     */
    public String getParameter(BarcodeArguments argument){
        return this.arguments.get(argument);
    }
    
    /**
     * adds the given parameter for the given argument
     * @param argument BarcodeArgument | the argument
     * @param parameter String | the parameter
     * @throws RuntimeException when the given argument is an invalid argument
     */
    public void setParameter(BarcodeArguments argument, String parameter){
        if (argument.equals(BarcodeArguments.INVALID_ARGUMENT)){
            throw new RuntimeException("Contains invalid arguments");
        }else {
            this.arguments.put(argument, parameter);
        }
    }
    
    /**
     * checks if all required parameters are set
     * @return the String "All values are set" if all required parameters are set, else an error message
     */
    @Override
    public String getErrorRequiredParametersSet(){
        String error = "";
        if (! this.arguments.containsKey(BarcodeArguments.FILENAME1)){
            error += "Value required for option " + BarcodeArguments.FILENAME1.getSortName() + "\n";
        }
        if (error.equals("")){
            error = "All values are set";
        }
        return error;
    }
    
    /**
     * checks if all required parameters are set
     * @return the boolean true if all required parameters are set
     */
    @Override
    public boolean areRequiredParametersSet(){
        if (! this.arguments.containsKey(BarcodeArguments.FILENAME1)){
            return false;
        }
        return true;
    }

    @Override
    public boolean containsParameter(Arguments argument) {
        if (argument instanceof BarcodeArguments){
            return this.containsParameter((BarcodeArguments) argument);
        }else{
            return false;
        }
    }

    @Override
    public String getParameter(Arguments argument) {
        if (argument instanceof BarcodeArguments){
            return this.getParameter((BarcodeArguments) argument);
        }else{
            return "";
        }
    }

    @Override
    public void setParameter(Arguments argument, String parameter) {
        if (argument instanceof BarcodeArguments){
            this.setParameter((BarcodeArguments) argument, parameter);
        }else{
            throw new RuntimeException("Contains invalid arguments");
        }
    }
    
    /**
     * 
     * @return String | the name of the file
     */
    public String getFileName(){
        return this.getParameter(BarcodeArguments.FILENAME1);
    }
    
    /**
     * 
     * @return int | the minimum length of the barcode
     */
    public int getMinimumLength(){
        return Integer.parseInt(this.getParameter(BarcodeArguments.MINIMUM_LENGTH));
    }
    
    /**
     * 
     * @return int | the maximum length of the barcode
     */
    public int getMaximumLength(){
        return Integer.parseInt(this.getParameter(BarcodeArguments.MAXIMUM_LENGTH));
    }
    
    /**
     * 
     * @return true if the file is a gzip file 
     */
    public boolean mustBeZipped(){
        if (this.getParameter(BarcodeArguments.USE_GZIP_FILES).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * 
     * @return String | the location of the output directory
     */
    public String getOutputDirectory(){
        return this.arguments.get(BarcodeArguments.OUTPUT_DIRECTORY);
    }
    
    /**
     * 
     * @return true if the enzyme must also be found
     */
    public boolean mustSearchEnzyme(){
        if (this.getParameter(BarcodeArguments.FIND_ENZYME).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * 
     * @return a Collection of all known enzymes
     * @see EnzymeCollection#getAllEnzymes() 
     */
    public Collection<Enzyme> getAllEnzymes(){
        return this.enzymeCollection.getAllEnzymes();
    }
    
    /**
     * 
     * @return the neutral enzyme
     * @see EnzymeCollection#getNeutralEnzyme() 
     */
    public Enzyme getNeutralEnzyme(){
        return this.enzymeCollection.getNeutralEnzyme();
    }
    
    /**
     * returns the enzyme with the given name
     * @param enzymeName String | the name of the enzyme
     * @return Enzyme
     */
    public Enzyme getEnzyme(String enzymeName){
        return this.enzymeCollection.getEnzyme(enzymeName);
    }
    
    /**
     * returns the enzyme collection
     * @return EnzymeCollection
     */
    public EnzymeCollection getEnzymeCollection(){
        return this.enzymeCollection;
    }
    
    /**
     * 
     * @return false if there aren't enzymes to be added, otherwise true
     */
    public boolean mustAddEnzymes(){
        if (this.arguments.get(BarcodeArguments.ENZYME_ADD).toLowerCase().equals("false")){
            return false;
        }else{
            return true;
        }
    }
    
    /**
     * 
     * @return String | the file name where the enzymes must be added from to the enzymeCollection
     */
    public String getAddEnzymeFile(){
        return this.arguments.get(BarcodeArguments.ENZYME_ADD);
    }
    
    /**
     * 
     * @return false if there aren't enzymes to be replaced, otherwise true
     */
    public boolean mustReplaceEnzymes(){
        if (this.arguments.get(BarcodeArguments.ENZYME_REPLACE).toLowerCase().equals("false")){
            return false;
        }else{
            return true;
        }
    }
    
    /**
     * 
     * @return String | the file name where the enzymes must be added from to the enzymeCollection
     */
    public String getReplaceEnzymeFile(){
        return this.arguments.get(BarcodeArguments.ENZYME_REPLACE);
    }
    
    /**
     * 
     * @return Integer | the minimum number that a barcode must be found before it is reported
     */
    public Integer getMinimumOccuranceBarcodes(){
        return Integer.parseInt(this.arguments.get(BarcodeArguments.MIN_BARCODE_OCCURANCE));
    }
    
    /**
     * 
     * @return Integer | the maximum number of barcodes to show in the output 
     */
    public Integer getMaximumNumbersBarcodes(){
        return Integer.parseInt(this.arguments.get(BarcodeArguments.MAX_NUMBER_OF_BARCODES));
    }
    
    /**
     * 
     * @return Integer | the percentage of mismatches in the barcodes (between 1 and 100)
     */
    public Integer getBarcodeMismatchPercentage(){
        return Integer.parseInt(this.arguments.get(BarcodeArguments.PERCENTAGE_OF_MISMATCH));
    }

    @Override
    public String getParametersLogString() {
        String toLog = "";
        toLog += "Parameters: \n";
        toLog += "Fastq File name: \t" + this.getFileName() + "\n";
        toLog += "Output directory= \t" + this.getOutputDirectory() + "\n";
        toLog += "Minimum barcode length: \t" + this.getMinimumLength() + "\n";
        toLog += "Maximum barcode length: \t" + this.getMaximumLength() + "\n";
        toLog += "Use gzip files: \t" + this.mustBeZipped() + "\n";
        if (this.mustAddEnzymes()){
            toLog += "\t Added enzymes of file: \t" + this.getAddEnzymeFile() + "\n";
        }
        if (this.mustReplaceEnzymes()){
            toLog += "\t Replaced enzymes with file: \t" + this.getReplaceEnzymeFile() + "\n";
        }
        toLog += "\t Minimum Occurance of the barcodes: \t" + this.getMinimumOccuranceBarcodes() + "\n";
        toLog += "\t Maximum Number of barcodes reported: \t" + this.getMaximumNumbersBarcodes() + "\n";
        toLog += "\t Mismatch percentage of the barcodes (between 1 and 100):\t" + this.getBarcodeMismatchPercentage() + "\n";
        toLog += "\t Find possible enzymes: \t" + this.mustSearchEnzyme() + "\n";
        return toLog;
    }

    @Override
    public String getParametersHelp() {
        String toHelp = "";
        toHelp += "Mandatory paramters:\n";
        toHelp += "\t -f1 \t the name of the input file \n";
        toHelp += "Optional parameters:\n";
        toHelp += "\t -min \t the minimum length of the barcode (standard 6) \n";
        toHelp += "\t -max \t the maximum length of the barcode (standard 16) \n";
        toHelp += "\t -gzip \t use gzip files as input and output (standard false) \n";
        toHelp += "\t -o \t the output directory (standard the directory of execution) \n";
        toHelp += "\t -ea \t Add enzymes from the given file (keeps the standard enzymes, and add the new) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites are comma separeted)) (only use once, not use -er)" + "\n";
        toHelp += "\t -er \t Replace enzymes from the given file (don't keep the standard enzymes) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites are comma separeted)) (only use once, not use -ea)" + "\n";
        toHelp += "\t -barmin \t The minimum occurance of a barcode before it is shown in the results (standard: 200) \n";
        toHelp += "\t -barmax \t The maximum of barcodes shown in the output (increasing this number will increase ram usage, but gives a slightly better result) (standard: 100)\n";
        toHelp += "\t -barmis \t The percentage of mismatches that may occure between barcodes (integer between 1 and 10) (standard: 10)\n";
        toHelp += "\t -fe \t Find possible enzyme combinations (standard: true)\n";
        return toHelp;
    }
    
    
    
}
