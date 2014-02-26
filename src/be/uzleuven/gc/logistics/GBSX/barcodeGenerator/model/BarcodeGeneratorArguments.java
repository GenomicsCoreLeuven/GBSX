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

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public enum BarcodeGeneratorArguments implements Arguments{

    /**
     * the number of barcodes needed to generate
     */
    NUMBER_OF_BARCODES ("-b"),
    /**
     * the enzyme used
     */
    ENZYME ("-e"),
    /**
     * add extra enzymes
     */
    ENZYME_FILE ("-ef"),
    /**
     * number of bootstraps to use
     */
    BOOTSTRAPS ("-nb"),
    /**
     * number of tries for each barcode, before stopping this bootstrap
     */
    BARCODE_TRIES ("-bt"),
    /**
     * the output directory
     */
    OUTPUT_DIRECTORY("-o"),
    /**
     * find the ultime match: best barcode combination with the counts (for the partial bases)
     */
    ULTIME_SEARCH_BARCODES ("-us"),
    /**
     * the file with barcodes that must be used to start the design from
     */
    BASIC_BARCODES ("-bf"),
    /**
     * the file with barcodes that not may be used
     */
    NOT_USE_BARCODES ("-nf"),
    
    /**
     * invalid argument
     */
    INVALID_ARGUMENT("");
    
    private final String name;
    
    private BarcodeGeneratorArguments(String name){
        this.name = name;
    }

    @Override
    public BarcodeGeneratorArguments getArgument(String sortName) {
        for (BarcodeGeneratorArguments arg : BarcodeGeneratorArguments.values()){
            if (arg.getSortName().equals(sortName)){
                return arg;
            }
        }
        return BarcodeGeneratorArguments.INVALID_ARGUMENT;
    }

    @Override
    public String getSortName() {
        return this.name;
    }
    
}
