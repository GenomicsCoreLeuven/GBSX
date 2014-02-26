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

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public enum BarcodeCorrectorArguments implements Arguments{

    /**
     * the fastq file
     */
    FILE ("-f"),
    /**
     * the info file (with the barcodes)
     */
    INFO_FILE ("-i"),
    /**
     * the input file and output file must be zipped
     */
    IS_ZIPPED ("-gz"),
    /**
     * the output directory
     */
    OUTPUT_DIRECTORY ("-o"),
    /**
     * the used enzyme
     */
    ENZYME ("-e"),
    /**
     * the file with extra enzymes
     */
    ENZYME_FILE ("-ef"),
    /**
     * invalid argument
     */
    INVALID_ARGUMENT("");

    private final String name;
    
    private BarcodeCorrectorArguments(String name){
        this.name = name;
    }
    
    @Override
    public BarcodeCorrectorArguments getArgument(String sortName) {
        for (BarcodeCorrectorArguments arg : BarcodeCorrectorArguments.values()){
            if (arg.getSortName().equals(sortName)){
                return arg;
            }
        }
        return BarcodeCorrectorArguments.INVALID_ARGUMENT;
    }

    @Override
    public String getSortName() {
        return this.name;
    }
    
}
