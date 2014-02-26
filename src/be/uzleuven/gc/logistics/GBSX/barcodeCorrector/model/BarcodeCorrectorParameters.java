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
package be.uzleuven.gc.logistics.GBSX.barcodeCorrector.model;

import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Arguments;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeCollection;
import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Parameters;
import java.util.HashMap;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeCorrectorParameters implements Parameters{
    
    private HashMap<BarcodeCorrectorArguments, String> arguments;
    private EnzymeCollection enzymeCollection;
    
    public BarcodeCorrectorParameters(){
        this.arguments = new HashMap();
        this.enzymeCollection = new EnzymeCollection();
        this.arguments.put(BarcodeCorrectorArguments.IS_ZIPPED, "true");
    }

    @Override
    public boolean containsParameter(Arguments argument) {
        if (argument instanceof BarcodeCorrectorArguments){
            return this.containsParameter((BarcodeCorrectorArguments) argument);
        }else{
            return false;
        }
    }
    
    public boolean containsParameter(BarcodeCorrectorArguments argument){
        return this.arguments.containsKey(argument);
    }

    @Override
    public String getParameter(Arguments argument) {
        if (argument instanceof BarcodeCorrectorArguments){
            return this.getParameter((BarcodeCorrectorArguments) argument);
        }else{
            return "ERROR";
        }
    }
    
    public String getParameter(BarcodeCorrectorArguments argument){
        return this.arguments.get(argument);
    }

    @Override
    public void setParameter(Arguments argument, String parameter) {
        if (argument instanceof BarcodeCorrectorArguments){
            this.setParameter((BarcodeCorrectorArguments) argument, parameter);
        }
    }
    
    public void setParameter(BarcodeCorrectorArguments argument, String parameter){
        this.arguments.put(argument, parameter);
    }

    @Override
    public boolean areRequiredParametersSet() {
        boolean isOk = true;
        if (! this.containsParameter(BarcodeCorrectorArguments.ENZYME)){
            isOk = false;
        }
        if (! this.containsParameter(BarcodeCorrectorArguments.FILE)){
            isOk = false;
        }
        if (! this.containsParameter(BarcodeCorrectorArguments.INFO_FILE)){
            isOk = false;
        }
        return isOk;
    }
    
    

    @Override
    public String getErrorRequiredParametersSet() {
        String error = "";
        if (! this.containsParameter(BarcodeCorrectorArguments.ENZYME)){
            error = "No enzyme is given. \n";
        }
        if (! this.containsParameter(BarcodeCorrectorArguments.FILE)){
            error = "No fastq file is given.\n";
        }
        if (! this.containsParameter(BarcodeCorrectorArguments.INFO_FILE)){
            error = "No info file is given.\n";
        }
        return error;
    }
    
    /**
     * 
     * @return String | the path to the fastq file
     */
    public String getFastqFile(){
        return this.getParameter(BarcodeCorrectorArguments.FILE);
    }
    
    /**
     * 
     * @return String | the path to the info file
     */
    public String getInfoFile(){
        return this.getParameter(BarcodeCorrectorArguments.INFO_FILE);
    }
    
    /**
     * 
     * @return String | the path to the output directory
     */
    public String getOutputDirectory(){
        return this.getParameter(BarcodeCorrectorArguments.OUTPUT_DIRECTORY);
    }
    
    /**
     * 
     * @return String | the path to the enzyme file
     */
    public String getEnzymeFile(){
        return this.getParameter(BarcodeCorrectorArguments.ENZYME_FILE);
    }
    
    /**
     * 
     * @return EnzymeCollection | the enzyme collection
     */
    public EnzymeCollection getEnzymeCollection(){
        return this.enzymeCollection;
    }
    
    /**
     * 
     * @return String | the name of the enzyme given by the user
     */
    public String getEnzymeName(){
        return this.getParameter(BarcodeCorrectorArguments.ENZYME);
    }
    
    /**
     * 
     * @return Enzyme | the used enzyme
     */
    public Enzyme getEnzyme(){
        return this.enzymeCollection.getEnzyme(this.getEnzymeName());
    }
    
    /**
     * 
     * @return true if the input and output must be zipped
     */
    public boolean mustBeZipped(){
        if (this.getParameter(BarcodeCorrectorArguments.IS_ZIPPED).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }

    /**
     * 
     * @return true if the parameters contains an enzyme file
     */
    public boolean mustAddEnzymes(){
        if (this.containsParameter(BarcodeCorrectorArguments.ENZYME_FILE)){
            return true;
        }else{
            return false;
        }
    }
    
    @Override
    public String getParametersLogString() {
        String toLog = "";
        toLog += "Used input file: " + this.getFastqFile() + "\n";
        toLog += "Used info file: " + this.getInfoFile() + "\n";
        toLog += "The input and output are zipped: " + this.mustBeZipped() + "\n";
        toLog += "The used output directory: " + this.getOutputDirectory() + "\n";
        toLog += "Used enzyme: " + this.getEnzyme().getName() + "\n";
        toLog += "Added enzymes: " + this.mustAddEnzymes() + "\n";
        if (this.mustAddEnzymes()){
            toLog += "Used enzyme file: " + this.getEnzymeFile() + "\n";
        }
        return toLog;
    }

    @Override
    public String getParametersHelp() {
        String toHelp = "";
        toHelp += "\t -f \t The input file (fastq or fastq.gz) (mandatory) \n";
        toHelp += "\t -i \t The info file, the used barcodes (mandatory) \n";
        toHelp += "\t -gz \t The input and output file (fastq) are zipped (standard true) \n";
        toHelp += "\t -o \t The output directory (standard the current working directory) \n";
        toHelp += "\t -e \t The used enzyme (no standard multiple enzyme support, solution: execute twice or use custom enzyme) \n";
        toHelp += "\t -ef \t An extra enzyme file, with custom enzymes (tab delimited file with as the first column the enzyme names, the second"
                + " column comma seperated enzyme cutsites remains. \n";
        return toHelp;
    }
    
}
