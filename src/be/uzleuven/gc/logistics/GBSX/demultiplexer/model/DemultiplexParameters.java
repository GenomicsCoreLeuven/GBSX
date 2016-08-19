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
package be.uzleuven.gc.logistics.GBSX.demultiplexer.model;

import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqScores;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.EnzymeCollection;
import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Parameters;
import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Arguments;
import java.util.Collection;
import java.util.HashMap;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class DemultiplexParameters implements Parameters{
    
    private HashMap<DemultiplexArguments, String> arguments;
    //standard barcode adaptor (can be changed by using parameter -ca)
    private final static String COMMON_AND_BARCODED_ADAPTOR = "AGATCGGAAGAGCG";
    private EnzymeCollection enzymeCollection;
    private boolean doubleBarcodes;
    
    public DemultiplexParameters(){
        this.arguments = new HashMap();
        //standard parameters
        this.arguments.put(DemultiplexArguments.OUTPUT_DIRECTORY, System.getProperty("user.dir"));
        this.arguments.put(DemultiplexArguments.ALLOWED_MISMATCHES_BARCODE, "1");
        this.arguments.put(DemultiplexArguments.ALLOWED_MISMATCHES_ENZYME, "1");
        this.arguments.put(DemultiplexArguments.START_BARCODE_DISTANCE, "0");
        this.arguments.put(DemultiplexArguments.COMPLETE_CHECK, "true");
        this.arguments.put(DemultiplexArguments.QUALITY_SCORE, FastqScores.ILLUMINA_1_8.getType());
        this.arguments.put(DemultiplexArguments.ADAPTOR_LIGASE, "1");
        this.arguments.put(DemultiplexArguments.MINIMUM_SEQUENCE_LENGTH, "0");
        this.arguments.put(DemultiplexArguments.KEEP_SEQUENCES_WITH_N, "true");
        this.arguments.put(DemultiplexArguments.COMMON_ADAPTOR, DemultiplexParameters.COMMON_AND_BARCODED_ADAPTOR);
        this.arguments.put(DemultiplexArguments.USE_GZIP_FILES, "false");
        this.arguments.put(DemultiplexArguments.KEEP_ENZYME_CUTSITES, "true");
        this.arguments.put(DemultiplexArguments.IS_RAD, "false");
        this.arguments.put(DemultiplexArguments.ENZYME_ADD, "false");
        this.arguments.put(DemultiplexArguments.ENZYME_REPLACE, "false");
        this.arguments.put(DemultiplexArguments.LONG_FILE_NAMES, "false");
        this.arguments.put(DemultiplexArguments.USE_SELF_CORRECTING_BARCODES, "false");
        this.arguments.put(DemultiplexArguments.COMMON_ADAPTOR_COMPARE_SIZE, "14");
        this.arguments.put(DemultiplexArguments.THREADS, "1");
        //this will automaticaly return the standard algorithm
        this.arguments.put(DemultiplexArguments.MISMATCH_ALGORITHM, FindingsAlgorithms.getStandard().getAlgorithmName());
        this.enzymeCollection = new EnzymeCollection();
        this.doubleBarcodes = false;
    }
    
    /**
     * Changes the double barcodes setting to true
     */
    public void setDoubleBarcodes(){
        this.doubleBarcodes = true;
    }
    
    /**
     * 
     * @return the double barcodes status, true if double barcodes are used
     */
    public boolean useDoubleBarcodes(){
        return this.doubleBarcodes;
    }
    
    /**
     * 
     * @param argument DemultiplexArguments
     * @return true if the given argument is a parameter
     */
    public boolean containsParameter(DemultiplexArguments argument){
        return this.arguments.containsKey(argument);
    }
    
    /**
     * 
     * @param argument DemultiplexArguments
     * @return the asked argument as a String
     */
    public String getParameter(DemultiplexArguments argument){
        return this.arguments.get(argument);
    }
    
    /**
     * adds the given parameter for the given argument
     * @param argument GBSargument | the argument
     * @param parameter String | the parameter
     * @throws RuntimeException when the given argument is an invalid argument
     */
    public void setParameter(DemultiplexArguments argument, String parameter){
        if (argument.equals(DemultiplexArguments.INVALID_ARGUMENT)){
            throw new RuntimeException("Contains invalid arguments");
        }else if(argument == DemultiplexArguments.QUALITY_SCORE){
            FastqScores fastqScores = FastqScores.getFastqScores(parameter);
            this.arguments.put(DemultiplexArguments.QUALITY_SCORE, fastqScores.getType());
        }else if (argument == DemultiplexArguments.COMMON_ADAPTOR){
            if (parameter.length() < 0){
                throw new RuntimeException("Invalid Common Adaptor");
            }else{
                this.arguments.put(argument, parameter);
                this.arguments.put(DemultiplexArguments.COMMON_ADAPTOR_COMPARE_SIZE, "" + parameter.length());
            }
        }else{
            this.arguments.put(argument, parameter);
        }
    }
    
    /**
     * checks if all required parameters are set
     * @return true if all parameters are set, else false
     */
    @Override
    public boolean areRequiredParametersSet(){
        if (! this.arguments.containsKey(DemultiplexArguments.FILENAME1)){
            return false;
        }
        if (! this.arguments.containsKey(DemultiplexArguments.INFO_FILE)){
            return false;
        }
        return true;
    }
    
    /**
     * checks if all required parameters are set
     * @return the String 'All values are set' if all required parameters are set, else an error message
     */
    @Override
    public String getErrorRequiredParametersSet(){
        String error = "";
        if (! this.arguments.containsKey(DemultiplexArguments.FILENAME1)){
            error += "Value required for option " + DemultiplexArguments.FILENAME1.getSortName() + "\n";
        }
        if (! this.arguments.containsKey(DemultiplexArguments.INFO_FILE)){
            error += "Value required for option " + DemultiplexArguments.INFO_FILE.getSortName() + "\n";
        }
        if (error.equals("")){
            error = "All values are set";
        }
        return error;
    }
    
    
    /**
     * 
     * @return int | the allowed mismatches for the barcode
     */
    public int getAllowedMismatchesBarcode(){
        return Integer.parseInt(this.arguments.get(DemultiplexArguments.ALLOWED_MISMATCHES_BARCODE));
    }
    
    /**
     * returns the sample specific mismatches
     * <br> if the sample barcode is 0 or positif, the sample barcode is returned
     * <br> else the standard or universal defined mismatches are used
     * @param sample Sample | the sample for the mismatches
     * @return int | the allowed mismatches for the barcode
     * @see DemultiplexParameters#getAllowedMismatchesBarcode() 
     */
    public int getAllowedMismatchesBarcode(Sample sample){
        if (sample.getBarcodeMismatches() >= 0){
            return sample.getBarcodeMismatches();
        }else{
            return this.getAllowedMismatchesBarcode();
        }
    }
    
    /**
     * returns true if the input is a gziped file, and the output must be a gziped file
     * @return true if the input file is gziped (.gz) and the output must be ziped
     */
    public boolean mustBeZipped(){
        if (this.arguments.get(DemultiplexArguments.USE_GZIP_FILES).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * returns true if the enzyme cutsites remains must be kept, false otherwise
     * @return true if the enzyme cutsites must be kept
     */
    public boolean keepCutSites(){
        if (this.arguments.get(DemultiplexArguments.KEEP_ENZYME_CUTSITES).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * return the extension of the input and output files 
     * <br> this is .fastq.gz for gziped files
     * <br> this is .fastq for regular files
     * @return String | the extension of the file
     */
    public String getFileExtension(){
        if (this.mustBeZipped()){
            return ".fastq.gz";
        }else{
            return ".fastq";
        }
    }
    
    /**
     * 
     * @return true if long file names must be used
     */
    public boolean useLongFileNames(){
        if (this.arguments.get(DemultiplexArguments.LONG_FILE_NAMES).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * 
     * @return int | the allowed mismatches for the enzyme
     */
    public int getAllowedMismatchesEnzyme(){
        return Integer.parseInt(this.arguments.get(DemultiplexArguments.ALLOWED_MISMATCHES_ENZYME));
    }
    
    /**
     * 
     * @return int | the max distance of the startsite (distance between the start of the read sequence and the barcode)
     */
    public int getStartDistance(){
        return Integer.parseInt(this.arguments.get(DemultiplexArguments.START_BARCODE_DISTANCE));
    }
    
    /**
     * 
     * @return String | the location of the output directory
     */
    public String getOutputDirectory(){
        return this.arguments.get(DemultiplexArguments.OUTPUT_DIRECTORY);
    }
    
    /**
     * 
     * @return String | the name of the info file
     */
    public String getInfoFile(){
        return this.arguments.get(DemultiplexArguments.INFO_FILE);
    }
    
    /**
     * 
     * @return String | the name of the first fastq file
     */
    public String getFastqFile1(){
        return this.arguments.get(DemultiplexArguments.FILENAME1);
    }
    
    /**
     * 
     * @return String | the name of the second fastq file, or null if only 1 file was given
     */
    public String getFastqFile2(){
        return this.arguments.get(DemultiplexArguments.FILENAME2);
    }
    
    /**
     * 
     * @return true if the complete check option is used
     */
    public boolean completeCheck(){
        if (this.arguments.get(DemultiplexArguments.COMPLETE_CHECK).toLowerCase().equals("true")) {
            return true;
        }
        return false;
    }
    
    /**
     * 
     * @return true if the sequences where N occurs as a nucleotide must be kept. 
     */
    public boolean keepSequencesWithN(){
        if (this.arguments.get(DemultiplexArguments.KEEP_SEQUENCES_WITH_N).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * 
     * @return FastqScores | the used fastq score
     */
    public FastqScores getFastqQualityScore(){
        return FastqScores.getFastqScores(this.arguments.get(DemultiplexArguments.QUALITY_SCORE));
    }
    
    /**
     * 
     * @return int | the number of mismatches in the adaptor
     */
    public int getAdaptorLigaseMismatches(){
        if (this.arguments.get(DemultiplexArguments.ADAPTOR_LIGASE).equals("no")){
            return -1;
        }else{
            return Integer.parseInt(this.arguments.get(DemultiplexArguments.ADAPTOR_LIGASE));
        }
    }
    
    /**
     * 
     * @return String | the common adaptor that is used in the GBS
     */
    public String getCommonAdaptor(){
        return this.arguments.get(DemultiplexArguments.COMMON_ADAPTOR);
    }
    
    /**
     * 
     * @return int | the minimum amount of basepairs in a sequence
     */
    public int getMinimumSequenceLength(){
        return Integer.parseInt(this.arguments.get(DemultiplexArguments.MINIMUM_SEQUENCE_LENGTH));
    }
    
    /**
     * 
     * @return int | the size of the common adapter to compare to
     */
    public int getAdaptorCompareSize(){
        return Integer.parseInt(this.arguments.get(DemultiplexArguments.COMMON_ADAPTOR_COMPARE_SIZE));
    }
    
    /**
     * returns true if the data is RAD data, returns false if the data is gbs data
     * @return true if the data is RAD data
     */
    public boolean isRadData(){
        if (this.arguments.get(DemultiplexArguments.IS_RAD).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * returns the findings algorithm that must be used for the mismatches and indels
     * @return 
     */
    public FindingsAlgorithms getFindingsAlgorithm(){
        return FindingsAlgorithms.getAlgorithm(this.arguments.get(DemultiplexArguments.MISMATCH_ALGORITHM));
    }
    
    /**
     * returns all known enzymes
     * @return a collection of all enzymes
     * @see EnzymeCollection#getAllEnzymes() 
     */
    public Collection<Enzyme> getAllEnzymes(){
        return this.enzymeCollection.getAllEnzymes();
    }
    
    /**
     * returns the neutral enzyme
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
        if (this.arguments.get(DemultiplexArguments.ENZYME_ADD).toLowerCase().equals("false")){
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
        return this.arguments.get(DemultiplexArguments.ENZYME_ADD);
    }
    
    /**
     * 
     * @return false if there aren't enzymes to be replaced, otherwise true
     */
    public boolean mustReplaceEnzymes(){
        if (this.arguments.get(DemultiplexArguments.ENZYME_REPLACE).toLowerCase().equals("false")){
            return false;
        }else{
            return true;
        }
    }
    
    /**
     * 
     * @return true if the barcodes are self correcting
     */
    public boolean useSelfCorrectingBarcodes(){
        if (this.arguments.get(DemultiplexArguments.USE_SELF_CORRECTING_BARCODES).toLowerCase().equals("true")){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * 
     * @return String | the file name where the enzymes must be added from to the enzymeCollection
     */
    public String getReplaceEnzymeFile(){
        return this.arguments.get(DemultiplexArguments.ENZYME_REPLACE);
    }
    
    /**
     * 
     * @return integer | the number of threads
     */
    public int getThreadNumber(){
        return Integer.parseInt(this.arguments.get(DemultiplexArguments.THREADS));
    }
    
    /**
     * configures a log file for all known parameters
     * @return String to put in the log file
     */
    @Override
    public String getParametersLogString(){
        String toLog = "";
        toLog += "Started to open the first files, parsing of the parameters succeded." + "\n";
        toLog += "Parameters: " + "\n";
        //-f1
        toLog += "\t First fastq file: \t" + this.getFastqFile1() + "\n";
        if (this.getFastqFile2() != null){
            //-f2
            toLog += "\t Second fastq file: \t" + this.getFastqFile2() + "\n";
        }
        //-i
        toLog += "\t Info-file: \t" + this.getInfoFile() + "\n";
        //-o
        toLog += "\t Output directory: \t" + this.getOutputDirectory() + "\n";
        //-gzip
        toLog += "\t In- and output file are gziped: \t" + this.mustBeZipped() + "\n";
        //-lf
        toLog += "\t Use long file names: \t" + this.useLongFileNames()+ "\n";
        //-rad
        toLog += "\t Data type: ";
        if (this.isRadData()){
            toLog += "\t RAD" + "\n";
        }else{
            toLog += "\t GBS" + "\n";
        }
        //-malg
        toLog += "\t Used algorithm to detect mismatches and/or indels: \t" + this.getFindingsAlgorithm().getAlgorithmName() + "\n";
        //-scb
        toLog += "\t Used self correcting barcodes: \t" + this.useSelfCorrectingBarcodes() + "\n";
        //-m or -mb
        toLog += "\t Allowed mismatches in the barcode: \t" + this.getAllowedMismatchesBarcode() + "\n";
        //-m or -me
        toLog += "\t Allowed mismatches in the enzyme: \t" + this.getAllowedMismatchesEnzyme() + "\n";
        //-p
        toLog += "\t Must check reads completely: \t" + this.completeCheck() + "\n";
        //-kc
        toLog += "\t Must keep the cutsites: \t" + this.keepCutSites() + "\n";
        //-al
        if (this.getAdaptorLigaseMismatches() == -1){
            toLog += "\t Check for Adaptor ligase: \t" + "no" + "\n";
        }else{
            toLog += "\t Allowed mismatches to check the adaptor ligase: \t" + this.getAdaptorLigaseMismatches() + "\n";
            //-ca
            toLog += "\t Used common adaptor: \t" + this.getCommonAdaptor() + "\n";
        }
        //-s
        toLog += "\t Maximum distance between start of sequence and start of barcode: \t" + this.getStartDistance() + "\n";
        //-q
        toLog += "\t Used scale for the quality: \t" + this.getFastqQualityScore().getType() + "\n";
        //-minsl
        toLog += "\t Minimum length of the sequence: \t" + this.getMinimumSequenceLength() + "\n";
        //-n
        toLog += "\t Keep sequences with N as nucleotide: \t" + this.keepSequencesWithN() + "\n";
        if (this.mustAddEnzymes()){
            toLog += "\t Added enzymes of file: \t" + this.getAddEnzymeFile() + "\n";
        }
        if (this.mustReplaceEnzymes()){
            toLog += "\t Replaced enzymes with file: \t" + this.getReplaceEnzymeFile() + "\n";
        }
        toLog += "\n";
        if (this.getFastqFile2() == null){
            toLog += "Single read demultiplexing" + "\n";
        }else{
            toLog += "Paired end demultiplexing" + "\n";
        }
        toLog += "Use double barcodes " + this.useDoubleBarcodes() + "\n";
        toLog += "\n";
        return toLog;
    }
    
    
    /**
     * makes a help for all known parameters
     * @return String | the help page
     */
    @Override
    public String getParametersHelp(){
        String toHelp = "";
        toHelp += "These parameters are mandatory: " + "\n";
        toHelp += "\t -f1 \t the name and path of the fastq or fastq.gz file to demultiplex" + "\n";
        toHelp += "\t -i \t the name and path of the info file. This is a tab delimeted file without headings, "
                + "with three (or more) columns: sample, sequence of the barcode, name of the enzyme, "
                + "name of the second enzyme (optional, can be an empty string), the second barcode (optional, can be an empty string),"
                + "mismatches for the barcode (optional)" + "\n";
        toHelp += "\n";
        toHelp += "These parameters are optional: " + "\n";
        toHelp += "\t -f2 \t the name of the second fastq or fastq.gz file (only with paired-end sequencing)" + "\n";
        toHelp += "\t -o \t the name of the output directory (standard the directory of the call)" + "\n";
        toHelp += "\t -lf \t use long file names (standard false) filename is standard the sample name, long file names is sample name _ barcode _ enzyme" + "\n";
        toHelp += "\t -rad \t if the data is rad data or not (-rad true for RAD data, -rad false for GBS data) standard false (GBS)" + "\n";
        toHelp += "\t -gzip \t the input and output are/must be gziped (.gz) (standard false: input and output are .fastq, if true this is .fastq.gz)" + "\n";
        toHelp += "\t -t \t the number of threads to use (standard 1)" + "\n";
        toHelp += "\t -mb \t the allowed mismatches in the barcodes (overrides the option -m)" + "\n";
        toHelp += "\t -me \t the allowed mismatches in the enzymes (overrides the option -m)" + "\n";
        toHelp += "\t -minsl \t the minimum allowed length for the sequences (standard 0, rejected sequences are found in the stats for each sample in the rejected.count column. The sequences self are found untrimmed in the undetermined file.)" + "\n";
        toHelp += "\t -n \t keep sequences where N occurs as nucleotide (standard true)" + "\n";
        toHelp += "\t -ca \t the common adaptor used in the sequencing (standard (only first piece) AGATCGGAAGAGCG) currently only used for adaptor ligase see -al and when -rad is true) (minimum length is 10)" + "\n";
        toHelp += "\t -s \t the posible distance of the start. This is the distance count from the start of the read to the first basepair of the barcode or enzyme (standard 0, maximum 20)" + "\n";
        //toHelp += "\t -cc \t Checks the complete read for the enzyme (if false, stops at the first possible enzyme cutsite) (use values true or false, standard is true) if used, the sequence after the enzyme site is compared to the adaptors, if the first basepairs of the sequence are compaired to the first basepairs of the adaptor" + "\n";
        toHelp += "\t -kc \t Keep the enzyme cut-site remains (standard true)" + "\n";
        toHelp += "\t -ea \t Add enzymes from the given file (keeps the standard enzymes, and add the new) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites are comma separeted)) (only use once, not use -er)" + "\n";
        toHelp += "\t -er \t Replace enzymes from the given file (don't keep the standard enzymes) (enzyme file: no header, enzyme name tab cutsites (multiple cutsites are comma separeted)) (only use once, not use -ea)" + "\n";
        //toHelp += "\t -al \t check for adaptor ligase: no (for no check) or a positive integer (starts at 0), for the number of mismatches (only checks 10 basepairs of the adaptor), standard 1" + "\n";
        toHelp += "\t -scb \t Use self correcting barcodes (barcodes created by the barcodeGenerator) (standard false)" + "\n";
        toHelp += "\t -malg \t the used algorithm to find mismatches and indels, possible algorithms (see README): " + "\n";
        for (FindingsAlgorithms algorithm : FindingsAlgorithms.values()){
            toHelp += "\t \t \t" + algorithm.getAlgorithmName();
            if (algorithm == FindingsAlgorithms.getStandard()){
                toHelp += " (Standard)";
            }
            toHelp += "\t" + algorithm.getDescription();
            toHelp += "\n";
        }
        toHelp += "\n";
        toHelp += "\t -q \t the kind of quality scores used in the fastq file (including how phred scores are encoded): " + "\n";
        for (FastqScores fastqScores : FastqScores.values()){
            toHelp += "\t \t \t" + fastqScores.getType();
            if (fastqScores == FastqScores.getStandard()){
                toHelp += " (Standard)";
            }
            toHelp += "\n";
        }
        toHelp += "\n";
        return toHelp;
    }

    @Override
    public boolean containsParameter(Arguments argument) {
        if (argument instanceof DemultiplexArguments){
            return this.containsParameter((DemultiplexArguments) argument);
        }else{
            return false;
        }
    }

    @Override
    public String getParameter(Arguments argument) {
        if (argument instanceof DemultiplexArguments){
            return this.getParameter((DemultiplexArguments) argument);
        }else{
            return "ERROR";
        }
    }

    @Override
    public void setParameter(Arguments argument, String parameter) {
        if (argument instanceof DemultiplexArguments){
            this.setParameter((DemultiplexArguments) argument, parameter);
        }
    }

   
    
    
    
}
