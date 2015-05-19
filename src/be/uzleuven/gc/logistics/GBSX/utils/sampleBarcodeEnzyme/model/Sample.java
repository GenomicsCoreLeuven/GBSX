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
package be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model;

import be.uzleuven.gc.logistics.GBSX.utils.enzyme.model.Enzyme;
import java.util.Collection;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class Sample implements Comparable<Sample> {
    
    private final String sampleID;
    private final Enzyme enzyme;
    private final Enzyme enzyme2;
    private final String barcode;
    private final String barcode2;
    private final String complementBarcode;
    private final HashSet<String> possibleEnzymeSites;
    private final HashSet<String> possibleEnzyme2Sites;
    private int barcodeMismatches;
    
    /**
     * 
     * @param sampleID String | the name of the sample
     * @param enzyme Enzyme | the enzyme used
     * @param barcode String | the barcode used
     */
    public Sample(String sampleID, Enzyme enzyme, String barcode){
        //sample ID
        this.sampleID = sampleID;
        //the used enzyme for this sample
        this.enzyme = enzyme;
        //the used second enzyme for this sample (for double digest)
        this.enzyme2 = enzyme;
        //the used barcode for this sample
        this.barcode = barcode;
        //ths used second barcode (here not used)
        this.barcode2 = null;
        //make the complement of the barcode
        this.complementBarcode = this.complementDNA(barcode);
        //saving the enzyme cut sites and the complement
        this.possibleEnzymeSites = new HashSet();
        for (String enzymeSite : enzyme.getInitialCutSiteRemnant()){
            this.possibleEnzymeSites.add(enzymeSite);
            this.possibleEnzymeSites.add(this.complementDNA(enzymeSite));
        }
        this.possibleEnzyme2Sites = new HashSet();
        for (String enzymeSite : enzyme.getInitialCutSiteRemnant()){
            this.possibleEnzyme2Sites.add(enzymeSite);
            this.possibleEnzyme2Sites.add(this.complementDNA(enzymeSite));
        }
        //save the mismatches of the barcode
        this.barcodeMismatches = -1;
    }
    
    /**
     * 
     * @param sampleID String | the name of the sample
     * @param enzyme Enzyme | the enzyme used
     * @param barcode String | the barcode used
     * @param barcode2 String | the second barcode used
     * @param barcodeMismatches int | the possible mismatches for the barcode for this sample
     */
    public Sample(String sampleID, Enzyme enzyme, Enzyme enzyme2, String barcode, String barcode2, int barcodeMismatches){
        //sample ID
        this.sampleID = sampleID;
        //the used enzyme for this sample
        this.enzyme = enzyme;
        //the used second enzyme for this sample (for double digest)
        this.enzyme2 = enzyme2;
        //the used barcode for this sample
        this.barcode = barcode;
        //ths used second barcode
        this.barcode2 = barcode2;
        //make the complement of the barcode
        this.complementBarcode = this.complementDNA(barcode);
        //saving the enzyme cut sites and the complement
        this.possibleEnzymeSites = new HashSet();
        for (String enzymeSite : enzyme.getInitialCutSiteRemnant()){
            this.possibleEnzymeSites.add(enzymeSite);
            this.possibleEnzymeSites.add(this.complementDNA(enzymeSite));
        }
        this.possibleEnzyme2Sites = new HashSet();
        for (String enzymeSite : enzyme2.getInitialCutSiteRemnant()){
            this.possibleEnzyme2Sites.add(enzymeSite);
            this.possibleEnzyme2Sites.add(this.complementDNA(enzymeSite));
        }
        //save the mismatches of the barcode
        this.barcodeMismatches = barcodeMismatches;
    }
    
    /**
     * 
     * @param sampleID String | the name of the sample
     * @param enzyme Enzyme | the enzyme used
     * @param barcode String | the barcode used
     * @param barcode2 String | the second barcode used
     */
    public Sample(String sampleID, Enzyme enzyme, Enzyme enzyme2, String barcode, String barcode2){
        //sample ID
        this.sampleID = sampleID;
        //the used enzyme for this sample
        this.enzyme = enzyme;
        //the used second enzyme for this sample (for double digest)
        this.enzyme2 = enzyme2;
        //the used barcode for this sample
        this.barcode = barcode;
        //ths used second barcode
        this.barcode2 = barcode2;
        //make the complement of the barcode
        this.complementBarcode = this.complementDNA(barcode);
        //saving the enzyme cut sites and the complement
        this.possibleEnzymeSites = new HashSet();
        for (String enzymeSite : enzyme.getInitialCutSiteRemnant()){
            this.possibleEnzymeSites.add(enzymeSite);
            this.possibleEnzymeSites.add(this.complementDNA(enzymeSite));
        }
        this.possibleEnzyme2Sites = new HashSet();
        for (String enzymeSite : enzyme2.getInitialCutSiteRemnant()){
            this.possibleEnzyme2Sites.add(enzymeSite);
            this.possibleEnzyme2Sites.add(this.complementDNA(enzymeSite));
        }
        //save the mismatches of the barcode
        this.barcodeMismatches = -1;
    }
    
    /**
     * 
     * @return String the id of the sample
     */
    public String getSampleID(){
        return this.sampleID;
    }
    
    /**
     * 
     * @return Enzyme | the enzyme of this sample 
     */
    public Enzyme getEnzyme(){
        return this.enzyme;
    }
    
    
    /**
     * 
     * @return Enzyme | the enzyme second of this sample 
     */
    public Enzyme getEnzyme2(){
        return this.enzyme2;
    }
    
    /**
     * 
     * @return String | the name of the enzyme
     */
    public String getEnzymeName(){
        return this.enzyme.getName();
    }
    
    
    /**
     * 
     * @return String | the name of the second enzyme
     */
    public String getEnzyme2Name(){
        return this.enzyme2.getName();
    }
    
    /**
     * 
     * @return Collection of String | the remains of the cutsite after the cut (place found in the DNA afther cut)
     * @see Enzyme#getInitialCutSiteRemnant() 
     */
    private Collection<String> getEnzymeCutSites(){
        return this.enzyme.getInitialCutSiteRemnant();
    }
    
    /**
     * 
     * @return Collection of String | the remains of the cutsite after the cut (place found in the DNA afther cut) of the second enzyme
     * @see Enzyme#getInitialCutSiteRemnant() 
     */
    private Collection<String> getEnzyme2CutSites(){
        return this.enzyme2.getInitialCutSiteRemnant();
    }
    
    /**
     * 
     * @return HashSet of String | all remains of the cutsite after the cut AND the complements of this site
     */
    private HashSet<String> getAllPossibleEnzymeCutSites(){
        return this.possibleEnzymeSites;
    }
    
    
    /**
     * 
     * @return HashSet of String | all remains of the cutsite after the cut AND the complements of this site
     */
    private HashSet<String> getAllPossibleEnzyme2CutSites(){
        return this.possibleEnzyme2Sites;
    }
    
    /**
     * 
     * @return the length of the first found enzyme site
     */
    public int getPossibleEnzymeCutSiteLength(){
        for (String site : this.getAllPossibleEnzymeCutSites()){
            return site.length();
        }
        return 0;
    }
    
    /**
     * 
     * @return the length of the first found enzyme site
     */
    public int getPossibleEnzyme2CutSiteLength(){
        for (String site : this.getAllPossibleEnzyme2CutSites()){
            return site.length();
        }
        return 0;
    }
    
    /**
     * 
     * @return String | the barcode used to identify the sample 
     */
    public String getBarcode(){
        return this.barcode;
    }
    
    /**
     * 
     * @return String | the second barcode used to identify the sample 
     */
    public String getBarcodeSecond(){
        return this.barcode2;
    }
    
    /**
     * 
     * @return boolean | true if there are 2 barcodes (on both ends), else false
     */
    public boolean has2barcodes(){
        return this.barcode2 != null;
    }
    
    /**
     * 
     * @return String | the complement of the barcode used to identify the sample
     */
    public String getComplementBarcode(){
        return this.complementBarcode;
    }
    
    /**
     * returns the mismatches the barcode can have, or -1 when the sample doesn't know how many mismatches (can be standard, or user defined in parameter)
     * @return int | the mismatches of the barcode
     */
    public int getBarcodeMismatches(){
        return this.barcodeMismatches;
    }
    
    
    /**
     * creates the complement of the given DNA
     * @param DNApiece String | the piece of DNA where the complement must be searched
     * @return the complement of the given DNA
     * @see BasePair#getComplementSequence(java.lang.String) 
     */
    private String complementDNA(String DNApiece){
        return BasePair.getComplementSequence(DNApiece);
    }
    
    /**
     * return true if the given barcode is valid.
     * <br> valid barcodes only exists of A, T, G and/or C
     * @param candidateBarcode String | the barcode to validate
     * @return true if the given candidate barcode is a valid barcode, false otherwise
     */
    public static boolean isValidBarcode(String candidateBarcode){
        Pattern pattern = Pattern.compile("[ATGCatgc]+");
        Matcher matcher = pattern.matcher(candidateBarcode);
        return matcher.matches();
    }

    /**
     * returns 0 if getSampleID == sample.getSampleID
     * returns positive int if this sampleID is alphabeticaly before the sample's sampleID
     * returns negative int if this sampleID is alphabeticly after the sample's sampleID
     * @param sample Sample | the sample to compare to
     * @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than the specified object.
     * @see String#compareTo(java.lang.String) 
     */
    public int compareTo(Sample sample) {
        return this.getSampleID().compareTo(sample.getSampleID());
    }
    
}
