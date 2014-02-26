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
package be.uzleuven.gc.logistics.GBSX.utils.argumentsAndParameters;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public interface Parameters {
    
    
    /**
     * 
     * @param argument 
     * @return true if the given argument is a parameter
     */
    public boolean containsParameter(Arguments argument);
    
    /**
     * 
     * @param argument 
     * @return the asked argument as a String
     */
    public String getParameter(Arguments argument);
    
    /**
     * adds the given parameter for the given argument
     * @param argument  | the argument
     * @param parameter String | the parameter
     * @throws RuntimeException when the given argument is an invalid argument
     */
    public void setParameter(Arguments argument, String parameter);
    
    /**
     * checks if all required parameters are set
     * @return the String "All values are set" if all required parameters are set, else an error message
     */
    public boolean areRequiredParametersSet();
    
    /**
     * checks if all required parameters are set
     * @return the boolean true if all required parameters are set
     */
    public String getErrorRequiredParametersSet();
    
    /**
     * configures a log file for all known parameters
     * @return String to put in the log file
     */
    public String getParametersLogString();
    
    /**
     * makes a help for all known parameters
     * @return String | the help page
     */
    public String getParametersHelp();
    
    /**
     * 
     * @return String | the output directory of the tool
     */
    public String getOutputDirectory();
}
