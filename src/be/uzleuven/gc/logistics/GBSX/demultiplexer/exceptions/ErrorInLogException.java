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
package be.uzleuven.gc.logistics.GBSX.demultiplexer.exceptions;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class ErrorInLogException extends Exception{
    
    private final String log;
    private final Exception exeption;
    
    public ErrorInLogException(String log, Exception exeption){
        this.log = log;
        this.exeption = exeption;
    }
    
    /**
     * 
     * @return the log tried to insert
     */
    public String getLog(){
        return this.log;
    }
    
    /**
     * 
     * @return Exception | the exception throwed in the LoggerFile 
     */
    public Exception getException(){
        return this.exeption;
    }
    
}
