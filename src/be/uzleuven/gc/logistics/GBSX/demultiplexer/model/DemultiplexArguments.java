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

import be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters.Arguments;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public enum DemultiplexArguments implements Arguments{
    
    
    /**
     * Looks if there occured an adapter ligase
     */
    ADAPTOR_LIGASE ("-al"),
    /**
     * The common adaptor wich must be used
     */
    COMMON_ADAPTOR ("-ca"),
    /**
     * The common adaptor compare size
     */
    COMMON_ADAPTOR_COMPARE_SIZE ("-cas"),
    /**
     * Checks the complete read for the enzyme (if false, stops at the first possible enzyme cutsite)
     */
    COMPLETE_CHECK ("-cc"),
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
     * The filename of the second read
     * Optional
     */
    FILENAME2 ("-f2"),
    /**
     * The info file (tab delimited file:
     * <br> first col: sample
     * <br> second col: barcode
     * <br> third col: enzyme
     * Mandatory
     */
    INFO_FILE ("-i"),
    /**
     * Keep the enzyme cut sites
     */
    KEEP_ENZYME_CUTSITES ("-kc"),
    /**
     * The used algorithm to search for mismatches
     */
    MISMATCH_ALGORITHM ("-malg"),
    /**
     * The allowed mismatches in a barcode enzyme combination
     */
    ALLOWED_MISMATCHES_BARCODE ("-mb"),
    /**
     * The allowed mismatches in the enzyme
     */
    ALLOWED_MISMATCHES_ENZYME ("-me"),
    /**
     * The minimum length of the sequence
     */
    MINIMUM_SEQUENCE_LENGTH ("-minsl"),
    /**
     * Keep sequences with N as nucleotide
     */
    KEEP_SEQUENCES_WITH_N ("-n"),
    /**
     * The place of the output
     * Optional
     */
    OUTPUT_DIRECTORY ("-o"),
    /**
     * The type of the quality score
     * @see FastqScores
     */
    QUALITY_SCORE ("-q"),
    /**
     * If the data is RAD data
     */
    IS_RAD ("-rad"),
    /**
     * The maximum number of basepairs between the start of the sequence and the barcode
     * Optional
     */
    START_BARCODE_DISTANCE ("-s"),
    /**
     * Use gzip files
     * Optional
     */
    USE_GZIP_FILES ("-gzip"),
    /**
     * use long file names
     */
    LONG_FILE_NAMES ("-lf"),
    /**
     * use self correcting barcodes (hamming code)
     */
    USE_SELF_CORRECTING_BARCODES ("-scb"),
    /**
     * number of threads
     */
    THREADS ("-t"),
    /**
     * If the given argument was invalid
     */
    INVALID_ARGUMENT ("ERROR");
    
    /**
     * returns the argument for the given name (-name), or invalid if no argument exists for the given name
     * @param sortName String | the name of the argument
     * @return GBSargument if any exists for the given name, else GBSaguments.INVALID_ARGUMENT
     */
    public DemultiplexArguments getArgument(String sortName){
        for (DemultiplexArguments arg : DemultiplexArguments.values()){
            if (sortName.equals(arg.getSortName())){
                return arg;
            }
        }
        return DemultiplexArguments.INVALID_ARGUMENT;
    }
    
    /**
     * the name of the argument (-name)
     */
    private final String sortname; 
    
    /**
     * 
     * @param sortname String | name of the argument
     */
    private DemultiplexArguments(String sortname){
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
