/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package be.uzleuven.gc.logistics.GBSX.simulator;

import be.uzleuven.gc.logistics.GBSX.utils.fasta.infrastructure.FastaBufferedReader;
import be.uzleuven.gc.logistics.GBSX.utils.fasta.model.FastaRead;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqBufferedWriter;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqPairBufferedWriter;
import be.uzleuven.gc.logistics.GBSX.utils.fastq.model.FastqRead;
import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.BasePair;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the Genomics Core Leuven
 */
public class GBSsimulate {
    
    public final static String VERSION = "GBS Data Simulator v2.0";
    public final static String LICENCE = "GLPv3";
    public final static boolean DEBUG = false;
    
    public static void main (String[] args){
        if (GBSsimulate.DEBUG){
            String input;
            input = "-help";
            args = input.split(" ");
        }
            
        GBSsimulate gbsSimulator = new GBSsimulate(args);
        gbsSimulator.simulate();
    }
    
    
    public GBSsimulate(String[] args){
        this.parseParameters(args);
    }
    
    ArrayList<Barcodes> barcodeList = new ArrayList();
    
    /**
     * parses the parameters and generates the help
     * @param args
     */
    private void parseParameters(String[] args){
        if (args.length == 0) {
            throw new IllegalArgumentException("\n\nNo arguments given.\n\n");
        }
        if (args[0].equals("version") || args[0].equals("-version") || args[0].equals("-v")){
            System.out.println(GBSsimulate.VERSION);
            System.exit(0);
        }
        if (args[0].equals("help") || args[0].equals("-help") || args[0].equals("-h")){
            GBSsimulate.getHelp();
            System.exit(0);
        }
        //parse parameters
        this.isPairedEnd = true; 
        this.commonAdapter = "AGATCGGAAGAGCG";
        this.readLength = 100;
        this.readDepth = 6;
        this.generateErrors = true;
        this.outputDir = null;
        this.fastaFilePath = null;
        this.barcodeFilePath = null;
        this.outputDir = System.getProperty("user.dir");
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
                    if (args[i].equals("-o")){
                        this.outputDir = args[i + 1];
                    }
                    if (args[i].equals("-f")){
                        this.fastaFilePath = args[i + 1];
                    }
                    if (args[i].equals("-p")){
                        if (args[i + 1].toLowerCase().equals("true")){
                            this.isPairedEnd = true;
                        }else{
                            this.isPairedEnd = false;
                        }
                    }
                    if (args[i].equals("-a")){
                        this.commonAdapter = args[i + 1];
                    }
                    if (args[i].equals("-l")){
                        this.readLength = Integer.parseInt(args[i + 1]);
                    }
                    if (args[i].equals("-rpl")){
                        this.readDepth = Integer.parseInt(args[i + 1]);
                    }
                    if (args[i].equals("-e")){
                        if (args[i + 1].toLowerCase().equals("true")){
                            this.generateErrors = true;
                        }else{
                            this.generateErrors = false;
                        }
                    }
                    if (args[i].equals("-b")){
                        this.barcodeFilePath = args[i + 1];
                    }
                }
            }
        }
        if (this.outputDir == null || this.fastaFilePath == null || this.barcodeFilePath == null){
            throw new RuntimeException();
        }
        File outDir = new File(this.outputDir);
        if (! outDir.exists()){
            outDir.mkdirs();
        }
        try {
            BufferedReader barcodeFileReader = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(this.barcodeFilePath))));
            String line;
            while((line = barcodeFileReader.readLine()) != null){
                String[] linesplit = line.split("\t");
                String barcode1 = linesplit[1];
                String barcode2 = "";
                if (linesplit.length >= 5){
                    barcode2 = linesplit[4];
                }
                this.barcodeList.add(new Barcodes(barcode1, barcode2));
            }
        
        } catch (FileNotFoundException ex) {
            System.exit(0);
        } catch (IOException ex) {
        }
    }
    
    private String outputDir;
    private String fastaFilePath;
    private boolean isPairedEnd;
    private boolean generateErrors;
    private String commonAdapter;
    private int readLength;
    private int readDepth;
    private String barcodeFilePath;
    
    /**
     * prints the help to the standard output (system.out.println)
     */
    public static void getHelp(){
        System.out.println("");
        System.out.println("");
        System.out.println("This is the help of the GBSsimulator.");
        System.out.println();
        System.out.println("Genomics Core Leuven");
        System.out.println("Licenced under " + GBSsimulate.LICENCE);
        System.out.println(GBSsimulate.VERSION);
        System.out.println();
        System.out.println("This program simulates GBS data, given a fasta file and a barcode file.");
        System.out.println("The program will output one or two fastq files (single/paired-end), "
                + "a file with the errors per barcode");
        System.out.println();
        System.out.println("\t-o\tOutput directory");
        System.out.println("\t-f\tfastafile (sequences in the fasta file must be orientated as enzyme1 to enzyme2)");
        System.out.println("\t-b\tbarcode file (output of the Barcode Generator)");
        System.out.println("\t-p\tis paired end (optional, standard true)(dual barcodes are only possible in the paired end mode)");
        System.out.println("\t-a\tcommon adapter (optional, AGATCGGAAGAGCG)");
        System.out.println("\t-l\tread length (optional, standard 100)");
        System.out.println("\t-rpl\tread per locus (optional, standard 6)");
        System.out.println("\t-e\terrors (optional, standard true)");
        System.out.println();
        System.out.println();
        System.out.println("Developed by the Genomics Core Leuven 2014");
        System.out.println("Licenced under " + GBSsimulate.LICENCE);
        System.out.println();
    }
    
    private FastaBufferedReader fastaReader;
    private ErrorSamples errorSamples;
    
    /**
     * prepares the simulation
     */
    public void simulate(){
        try {
            File fastaFile = new File(this.fastaFilePath);
            this.fastaReader = new FastaBufferedReader(fastaFile, false);
            this.errorSamples = new ErrorSamples();
            if (this.isPairedEnd){
                this.simulatePaired();
            }else{
                this.simulateSingle();
            }
        } catch (Exception ex) {
            Logger.getLogger(GBSsimulate.class.getName()).log(Level.SEVERE, null, ex);
            System.out.println("ENDED with major errors");
        }
    }
    
    /**
     * Reads the fasta file one by one
     * creates fastqreads for both possible ends (forward, and reverse)
     * @throws FileNotFoundException
     * @throws IOException 
     */
    private void simulateSingle() throws FileNotFoundException, IOException{
        String filePathFastq = this.outputDir + System.getProperty("file.separator") + "output.fastq";
        FastqBufferedWriter outputWriter = new FastqBufferedWriter(new File(filePathFastq), false);
        FastaRead fastaRead;
        while ((fastaRead = this.fastaReader.next()) != null){
            for (Barcodes barcode : this.barcodeList){
                String originalSequence = fastaRead.getSequence();
                for (int i=0; i < this.readDepth; i++){
                    //normal sequence
                    String sequenceToUse = originalSequence;
                    String sequence = barcode.getBarcode1() + sequenceToUse + this.commonAdapter;
                    FastqRead fastq = this.createFastq(fastaRead.getDescription() + "-" + i + "|" + barcode.getBarcode1(), sequence, barcode.getBarcode1(), false);
                    outputWriter.write(fastq);
                }
            }
        }
        outputWriter.close();
        List<String> tmpBarcodeList = new ArrayList();
        for (Barcodes barcode : this.barcodeList){
            tmpBarcodeList.add(barcode.getBarcode1());
        }
        Collections.sort(tmpBarcodeList);
        String[] barcodeOrder = tmpBarcodeList.toArray(new String[tmpBarcodeList.size()]);
        BufferedWriter errorWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(new File(this.outputDir + System.getProperty("file.separator") + "output.error.txt")))));
        errorWriter.write(this.errorSamples.getCompleteErrorReport(barcodeOrder));
        errorWriter.close();
        
    }
    
    /**
     * Reads the fasta file one by one
     * creates fastqreads for both possible ends and pair end (forward, and reverse)
     * @throws FileNotFoundException
     * @throws IOException 
     */
    

    private void simulatePaired() throws FileNotFoundException, IOException{
        String filePathFastq1 = this.outputDir + System.getProperty("file.separator") + "output.R1.fastq";
        String filePathFastq2 = this.outputDir + System.getProperty("file.separator") + "output.R2.fastq";
        FastqPairBufferedWriter outputWriter = new FastqPairBufferedWriter(new File(filePathFastq1), new File(filePathFastq2), false);
        FastaRead fastaRead;
        while ((fastaRead = this.fastaReader.next()) != null){
            for (Barcodes barcode : this.barcodeList){
                String originalSequence = fastaRead.getSequence();
                for (int i=0; i < this.readDepth; i++){
                    String sequenceToUse = originalSequence;
                    String barcodeSeq = barcode.getBarcode1() + sequenceToUse + BasePair.getComplementSequence(barcode.getBarcode2());
                    String sequence = barcodeSeq + this.commonAdapter;
                    String sequence2 = BasePair.getComplementSequence(barcodeSeq) + this.commonAdapter;
                    
                    FastqRead fastq1 = this.createFastq(fastaRead.getDescription() + "-" + i + "|" + barcode.getBarcode1() + "_" + barcode.getBarcode2(), 
                            sequence, barcode.getBarcode1() + "_" + barcode.getBarcode2(), false);
                    FastqRead fastq2 = this.createFastq(fastaRead.getDescription() + "-" + i + "|" + barcode.getBarcode1() + "_" + barcode.getBarcode2(),
                            sequence2, barcode.getBarcode1() + "_" + barcode.getBarcode2(), true);
                    outputWriter.write(fastq1, fastq2);
                }
            }
        }
        outputWriter.close();
        List<String> tmpBarcodeList = new ArrayList();
        for (Barcodes barcode : this.barcodeList){
            tmpBarcodeList.add(barcode.getBarcode1() + "_" + barcode.getBarcode2());
        }
        Collections.sort(tmpBarcodeList);
        String[] barcodeOrder = tmpBarcodeList.toArray(new String[tmpBarcodeList.size()]);
        BufferedWriter errorWriter = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(new File(this.outputDir + System.getProperty("file.separator") + "output.error.txt")))));
        errorWriter.write(this.errorSamples.getCompleteErrorReport(barcodeOrder));
        errorWriter.close();
    }
    
    
    /**
     * cuts the given sequence on the disired length,
     * add possible errors in the sequence, according to the error rate of the sample
     * creates an artificial quality for each base
     * 
     * @param description String | the description of the read
     * @param sequence String | the sequence to generate the new FastqRead
     * @param barcode String | the barcode
     * @param isSecondRead boolean | true if this is the second read of pair end data
     * @return FastqRead | the new made fastq
     */
    private FastqRead createFastq(String description, String sequence, String barcode, boolean isSecondRead){
        if (sequence.length() < this.readLength){
            for (int i=sequence.length(); i < this.readLength; i++){
                sequence = sequence + "A";
            }
        }
        String newSeq = sequence.substring(0, this.readLength);
        String newQual = "";
        ArrayList<Integer> errorPositions = new ArrayList();
        for (int qi=0; qi < this.readLength; qi++){
            int qual = 30;
            Random r = new Random();
            //change on mutation
            if (this.generateErrors && r.nextInt(100) <= 1){
                //error
                char[] bases = {'A', 'C', 'G', 'T', 'N'};
                String correctedSeq = newSeq;
                char newbase = bases[r.nextInt(bases.length)];
                if (newbase == correctedSeq.charAt(qi)){
                }else{
                    errorPositions.add(qi);
                }
                if (qi == 0){
                    correctedSeq = newbase + newSeq.substring(1);
                }else if (qi == this.readLength - 1){
                    correctedSeq = newSeq.substring(0, newSeq.length() - 1) + newbase;
                }else{
                    correctedSeq = newSeq.substring(0, qi) + newbase + newSeq.substring(qi + 1);
                }
                newSeq = correctedSeq;
                int score = r.nextInt(qual);
                newQual += "" + Character.toChars((score - 0 + 33))[0];
            }else{
                int score = r.nextInt(41 - qual) + qual;
                newQual += "" + Character.toChars((score - 0 + 33))[0];
            }
        }
        this.errorSamples.addErrors(barcode, errorPositions, isSecondRead);
        return new FastqRead(description, newSeq.toUpperCase(), newQual);
    }
    
    
    
    /**
     * used to keep track of the errors of the samples
     */
    private class ErrorSamples{
        
        private HashMap<String, SampleError> errorMap;
        
        public ErrorSamples(){
            this.errorMap = new HashMap();
        }
        
        /**
         * add a new sample
         * @param barcode | String the barcode of the sample
         */
        public void addSample(String barcode){
            this.errorMap.put(barcode, new SampleError(barcode));
        }
        
        /**
         * add new errors to the given sample
         * @param barcode | String the barcode of the sample
         * @param errorPositions | ArrayList of Integers with the error positions 
         * @param isSecondRead | boolean true if second read of pair-end
         */
        public void addErrors(String barcode, ArrayList<Integer> errorPositions, boolean isSecondRead){
            if (! this.errorMap.containsKey(barcode)){
                this.addSample(barcode);
            }
            this.errorMap.get(barcode).addErrorPositions(errorPositions, isSecondRead);
        }
        
        /**
         * 
         * @param barcodeOrder | String[] of all barcodes (samples)
         * @return String with a summary of all samples and errors
         * @see SampleError#toSummaryLine() 
         */
        public String getCompleteErrorReport(String[] barcodeOrder){
            String line = new SampleError("").getHeading() + "\n";
            for (String b : barcodeOrder){
                if (! this.errorMap.containsKey(b)){
                    line += b + "\t" + "NOT FOUND\n";
                }else{
                    line += this.errorMap.get(b).toSummaryLine() + "\n";
                }
            }
            return line;
        }
        
        private class SampleError{
            
            private String barcode;
            private int correctCount;
            private int correctCountInBarcode;
            private int errorCountInBarcode;
            private int errorCountInSequences;
            private int multipleErrors;
            
            /**
             * create a new Sample error tracking
             * @param barcode | String the barcode of the sample
             */
            public SampleError(String barcode){
                this.barcode = barcode;
                this.correctCount = 0;
                this.correctCountInBarcode = 0;
                this.errorCountInBarcode = 0;
                this.errorCountInSequences = 0;
                this.multipleErrors = 0;
            }
            
            /**
             * adds new error positions
             * @param errorPos | arraylist of integers with all error positions
             * @param isSecondRead | boolean true if this is the second read of paired-end data
             */
            public void addErrorPositions(ArrayList<Integer> errorPos, boolean isSecondRead){
                if (errorPos.isEmpty() && !isSecondRead){
                    this.correctCount++;
                    this.correctCountInBarcode++;
                }else{
                    boolean barcodeErrorFound = false;
                    boolean multiple = false;
                    for (int pos : errorPos){
                        if (isSecondRead){
                            this.errorCountInSequences++;
                        }else{
                            if (pos < this.barcode.length()){
                                if (!barcodeErrorFound){
                                    this.errorCountInBarcode++;
                                    barcodeErrorFound = true;
                                }else if (! multiple){
                                    this.errorCountInBarcode--;
                                    this.multipleErrors++;
                                    multiple = true;
                                }
                            }
                            this.errorCountInSequences++;
                        }
                    }
                    if (!isSecondRead && !barcodeErrorFound){
                        this.correctCountInBarcode++;
                    }
                }
            }
            
            /**
             * 
             * @return String | the summary of the sample (tabdelimited): 
             * barcode, barcode, number of correct barcodes, number of barcodes with an error, number of reads with multiple errors, number of complete correct reads, total number of errors
             */
            public String toSummaryLine(){
                return this.barcode + "\t" + this.barcode + "\t" + this.correctCountInBarcode + "\t" + this.errorCountInBarcode + "\t" + this.multipleErrors + "\t" + this.correctCount + "\t" + this.errorCountInSequences;
            }
            
            public String getHeading(){
                return "SampleName\tBarcode\tCorrectBarcodes\tErrorBarcodes\tSamplesWithMultipleErrors\tCompleteCorrectSequences\tTotalErrors";
            }
        }
        
        
    }
    
    private class Barcodes{
        
        private final String barcode1;
        private final String barcode2;
        
        public Barcodes(String barcode){
            this.barcode1 = barcode;
            this.barcode2 = "";
        }
        
        public Barcodes(String barcode1, String barcode2){
            this.barcode1 = barcode1;
            this.barcode2 = barcode2;
        }
        
        public String getBarcode1(){
            return this.barcode1;
        }
        
        public String getBarcode2(){
            return this.barcode2;
        }
        
    }
    
}
