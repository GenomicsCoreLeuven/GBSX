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

import be.uzleuven.gc.logistics.GBSX.utils.sampleBarcodeEnzyme.model.BasePair;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class EnzymeClass implements Enzyme{
    
    
    private String name;
    private String[] initialCutSiteRemnants;
    
    /**
     * creates a new Enzyme from the name and the cut sites (as String array)
     * @param name String | the name of the Enzyme
     * @param initialCutSiteRemnants String[] | all the cut sites
     */
    public EnzymeClass(String name, String[] initialCutSiteRemnants){
        this.name = name;
        for (int i = 0; i < initialCutSiteRemnants.length; i++){
            initialCutSiteRemnants[i] = initialCutSiteRemnants[i].toUpperCase();
        }
        this.initialCutSiteRemnants = initialCutSiteRemnants;
    }
    
    /**
     * creates a new Enzyme from the name and the cut sites (as a collection)
     * @param name String | the name of the Enzyme
     * @param initialCutSiteRemnants collection of String | all the cut sites
     */
    public EnzymeClass(String name, Collection<String> initialCutSiteRemnants){
        this.name = name;
        this.initialCutSiteRemnants = (String[]) initialCutSiteRemnants.toArray();
        for (int i = 0; i < this.initialCutSiteRemnants.length; i++){
            this.initialCutSiteRemnants[i] = this.initialCutSiteRemnants[i].toUpperCase();
        }
    }

    /**
     * 
     * @return String | the name of the Enzyme
     */
    public String getName() {
        return this.name;
    }
    
    /**
     * returns all cutsites
     * @return Collection of String | all the cutsites of the enzyme in a collection
     */
    public Collection<String> getInitialCutSiteRemnant(){
        HashSet<String> cutsitesSet = new HashSet();
        cutsitesSet.addAll(Arrays.asList(this.initialCutSiteRemnants));
        return cutsitesSet;
    }
    
    /**
     * returns all complement cutsites
     * @return Collection of String | all the complement cutsites of the enzyme in a collection
     */
    public Collection<String> getComplementCutSiteRemnant() {
        HashSet<String> cutsitesSet = new HashSet();
        for (String cutsite : this.initialCutSiteRemnants){
            cutsitesSet.add(BasePair.getComplementSequence(cutsite));
        }
        return cutsitesSet;
    }
    
}
