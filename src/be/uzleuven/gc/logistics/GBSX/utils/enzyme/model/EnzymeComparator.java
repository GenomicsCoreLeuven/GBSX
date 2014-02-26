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
package be.uzleuven.gc.logistics.GBSX.utils.enzyme.model;

import java.util.Comparator;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class EnzymeComparator implements Comparator<Enzyme>{

    
    /**
     * 
     * @param enzyme1 Enzyme
     * @param enzyme2 Enzyme
     * @return a negative int if the first enzyme is less than the second,
     *          zero when both enzymes are equal
     *          a positive int when the first enzyme is greater than the second
     */
    public int compare(Enzyme enzyme1, Enzyme enzyme2) {
        return enzyme1.getName().compareTo(enzyme2.getName());
    }
    
    
}
