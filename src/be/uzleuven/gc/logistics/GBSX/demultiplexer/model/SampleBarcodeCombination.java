/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package be.uzleuven.gc.logistics.GBSX.demultiplexer.model;

import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.Sample;
/**
 *
 * @author koen
 */
public class SampleBarcodeCombination {

    /**
     * a combination of a sample, a barcode and enzyme, the mismatches (between sequence and barcode/enzyme) and the start location of the barcode in the sequence
     */
        
    private final Sample sample;
    private final String enzymeCutSite;
    private final int location;
    private final int mismatches;
    private final int lengthFoundBarcode;
    private final int lengthFoundEnzyme;
    /**
     * 
     * @param sample Sample | the sample of this combination
     * @param enzymeCutSite String | the used enzyme cutsite
     * @param location int | the start location of the barcodeEnzyme
     * @param mismatches int | the amount of mismatches
     * @param lengthFoundBarcode int | length of the found barcode
     * @param lengthFoundEnzyme int | length of the found enzyme
     */
    public SampleBarcodeCombination(Sample sample, String enzymeCutSite, int location, int mismatches, int lengthFoundBarcode, int lengthFoundEnzyme){
        this.sample = sample;
        this.enzymeCutSite = enzymeCutSite;
        this.location = location;
        this.mismatches = mismatches;
        this.lengthFoundBarcode = lengthFoundBarcode;
        this.lengthFoundEnzyme = lengthFoundEnzyme;
    }

    /**
     * 
     * @return Sample | the sample of this combination
     */
    public Sample getSample(){
        return this.sample;
    }

    /**
     * 
     * @return String | the enzymeCutsite
     */
    public String getEnzymeCutsite(){
        return this.enzymeCutSite;
    }

    /**
     * 
     * @return int | the location of the barcode + enzyme in a sequence
     */
    public int getLocation(){
        return this.location;
    }

    /**
     * 
     * @return int | the amount of mismatches between the barcode + enzyme and the sequence
     */
    public int getMismatches(){
        return this.mismatches;
    }

    /**
     * 
     * @return int | the length of the found enzyme
     */
    public int getLengthFoundEnzyme(){
        return this.lengthFoundEnzyme;
    }

    /**
     * 
     * @return int | the lenght of the found barcode
     */
    public int getLengthFoundBarcode(){
        return this.lengthFoundBarcode;
    }

}
