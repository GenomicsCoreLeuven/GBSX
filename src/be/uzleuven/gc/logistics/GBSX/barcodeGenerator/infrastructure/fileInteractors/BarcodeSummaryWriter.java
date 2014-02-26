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
package be.uzleuven.gc.logistics.GBSX.barcodeGenerator.infrastructure.fileInteractors;

import be.uzleuven.gc.logistics.GBSX.barcodeGenerator.model.BarcodeGeneratorParameters;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class BarcodeSummaryWriter {
    
    
    private String outputDirectory;
    private BarcodeGeneratorParameters barcodeGeneratorParameters = null;
    
    public BarcodeSummaryWriter(String outputDirectory){
        File outDir = new File(outputDirectory);
        if (! outDir.exists()){
            outDir.mkdirs();
        }
        this.outputDirectory = outputDirectory;
    }
    
    public BarcodeSummaryWriter(BarcodeGeneratorParameters barcodeGeneratorParameters){
        this.outputDirectory = barcodeGeneratorParameters.getOutputDirectory();
        File outDir = new File(outputDirectory);
        if (! outDir.exists()){
            outDir.mkdirs();
        }
        this.barcodeGeneratorParameters = barcodeGeneratorParameters;
    }
    
    
    /**
     * 
     * @param barcodeCollection collection of String | all barcodes that must be writen to a file
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void write(Collection<String> barcodeCollection, double[][] positionsCount, HashMap<Character, Integer> basesMap, double[] quality, int bootstraps) throws FileNotFoundException, IOException{
        ArrayList<String> barcodes = new ArrayList(barcodeCollection);
        Collections.sort(barcodes);
        //make buffered writer
        File file = new File(this.outputDirectory + System.getProperty("file.separator") + "barcode_summary.txt");
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new DataOutputStream(new FileOutputStream(file))));
        //write every barcode
        writer.write("Barcode Generator\n");
        writer.write("-----------------\n");
        writer.write("\n");
        if (this.barcodeGeneratorParameters != null){
            writer.write(this.barcodeGeneratorParameters.getParametersLogString());
            writer.write("\n");
        }
        writer.write("Barcodes: \n");
        if (this.barcodeGeneratorParameters.hasBasicSet()){
            writer.write("Old Barcodes: \n");
            for (String barcode : this.barcodeGeneratorParameters.getBasicSet()){
                writer.write("\t" + barcode + "\n");
            }
            writer.write("New Barcodes: \n");
        }
        for (String b : barcodes){
            writer.write("\t" + b + "\n");
        }
        writer.write("\n\n");
        writer.write("Position matrix:\n");
        ArrayList<Character> basesList = new ArrayList(basesMap.keySet());
        Collections.sort(basesList);
        for (int i = 0; i < positionsCount.length; i++){
            writer.write("\t\t" + (i + 1));
        }
        writer.write("\n");
        for (char base : basesList){
            writer.write("" + base);
            for (int i = 0; i < positionsCount.length; i++){
                writer.write("\t\t" + positionsCount[i][basesMap.get(base)]);
            }
            writer.write("\n");
        }
        writer.write("Qual:");
        for (int i = 0; i < quality.length; i++){
            DecimalFormat df = new DecimalFormat("##0.00");
            writer.write("\t\t" + df.format(quality[i]));
        }
        writer.write("\n");
        writer.write("Created " + barcodeCollection.size() + " of " + this.barcodeGeneratorParameters.getNumberOfBarcodes() + " barcodes.\n");
        writer.write("Used " + bootstraps + " of " + this.barcodeGeneratorParameters.getNumberOfBootstraps() + " bootstraps to create this file.\n");
        writer.write("\n");
        writer.close();
    }
}
