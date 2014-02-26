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

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public enum BarcodeArguments implements Arguments{
    
    /**
     * Add enzymes from the given file
     */
    ENZYME_ADD ("-ea"),
    /**
     * Replace systems enzyme with the given file
     */
    ENZYME_REPLACE ("-er"),
    /**
     * The filename of the first read
     * Mandatory
     */
    FILENAME1 ("-f1"),
    /**
     * find also the possible enzymes
     */
    FIND_ENZYME ("-fe"),
    /**
     * minimum length of the barcode
     */
    MINIMUM_LENGTH ("-min"),
    /**
     * maximum length of the barcode
     */
    MAXIMUM_LENGTH ("-max"),
    /**
     * The place of the output
     * Optional
     */
    OUTPUT_DIRECTORY ("-o"),
    /**
     * Use gzip files
     * Optional
     */
    USE_GZIP_FILES ("-gzip"),
    /**
     * the minimum number of barcode occurances
     */
    MIN_BARCODE_OCCURANCE ("-barmin"),
    /**
     * the maximum number of barcodes to report
     */
    MAX_NUMBER_OF_BARCODES ("-barmax"),
    /**
     * the percentage of mismatch between barcodes (between 0-100)
     */
    PERCENTAGE_OF_MISMATCH ("-barmis"),
    /**
     * If the given argument was invalid
     */
    INVALID_ARGUMENT("ERROR");
    
    /**
     * returns the argument for the given name (-name), or invalid if no argument exists for the given name
     * @param sortName String | the name of the argument
     * @return GBSargument if any exists for the given name, else GBSaguments.INVALID_ARGUMENT
     */
    public BarcodeArguments getArgument(String sortName){
        for (BarcodeArguments arg : BarcodeArguments.values()){
            if (sortName.equals(arg.getSortName())){
                return arg;
            }
        }
        return BarcodeArguments.INVALID_ARGUMENT;
    }
    
    /**
     * the name of the argument (-name)
     */
    private final String sortname; 
    
    /**
     * 
     * @param sortname String | name of the argument
     */
    private BarcodeArguments(String sortname){
        this.sortname = sortname;
    }
    
    /**
     * 
     * @return String | the option name of the argument
     */
    public String getSortName(){
        return this.sortname;
    }
}
