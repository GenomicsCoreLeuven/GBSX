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

import java.util.Collection;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public interface Enzyme{    
    /**
     * 
     * @return String the name of the enzyme
     */
    public String getName();
    
    /**
     * 
     * @return collection of String the sites after the cut
     */
    public Collection<String> getInitialCutSiteRemnant();
    
    /**
     * 
     * @return collection of String the complement sites after the cut
     */
    public Collection<String> getComplementCutSiteRemnant();
    
}
