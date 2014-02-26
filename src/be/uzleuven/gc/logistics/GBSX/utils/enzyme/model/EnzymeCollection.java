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
import java.util.HashSet;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class EnzymeCollection {
    
    
    HashSet<Enzyme> enzymeSet = new HashSet();
    
    /**
     * creates a enzymeCollection from the known system enzymes
     */
    public EnzymeCollection(){
        for (Enzyme enzyme : EnzymeEnum.values()){
            this.enzymeSet.add(enzyme);
        }
    }
    
    /**
     * creates a new enzymeCollection from the given enzymes + the neutral enzyme
     * @param enzymeCollection collection of enzymes
     */
    public EnzymeCollection(Collection<Enzyme> enzymeCollection){
        this.enzymeSet = new HashSet(enzymeCollection);
        this.enzymeSet.add(this.getNeutralEnzyme());
    }
    
    /**
     * adds all given enzymes to the collection
     * @param enzymeCollection | collection of enzymes
     */
    public void addEnzymes(Collection<Enzyme> enzymeCollection){
        this.enzymeSet.addAll(enzymeCollection);
    }
    
    /**
     * replaces the enzymes from this collection with the new collection
     * @param enzymeCollection | collection of enzymes
     */
    public void replaceEnzymes(Collection<Enzyme> enzymeCollection){
        this.enzymeSet = new HashSet(enzymeCollection);
    }
    
    /**
     * 
     * @return a Collection of all known enzymes
     */
    public Collection<Enzyme> getAllEnzymes(){
        return this.enzymeSet;
    }
    
    /**
     * 
     * @param name String | the name of the enzyme
     * @return Enzyme | the enzyme with the given name
     */
    public Enzyme getEnzyme(String name){
        for (Enzyme enzyme : this.enzymeSet){
            if (enzyme.getName().toLowerCase().equals(name.toLowerCase())){
                return enzyme;
            }
        }
        return EnzymeEnum.NAN;
    }
    
    /**
     * returns the neutral NAN enzyme
     * @return Enzyme | a neutral nan enzyme
     * @see EnzymeEnum#NAN
     */
    public Enzyme getNeutralEnzyme(){
        return EnzymeEnum.NAN;
    }
    
    
    
}
