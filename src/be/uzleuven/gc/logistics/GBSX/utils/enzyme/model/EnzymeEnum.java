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
public enum EnzymeEnum implements Enzyme{
    
    /*
     * TO ADD ENZYMES:
     * first name of the enzyme (caps) then the name, and the cutsites
     * copy:
     * NAME ("Name", new String[] {"first cutsite remains", "second cutsite remains", "third cutsite remains", ...}),
     */
    
    APEKI ("ApeKI", new String[] {"CAGC", "CTGC"}),
    
    PSTI ("PstI", new String[] {"TGCAG"}),
    
    ECOT22I ("EcoT22I", new String[] {"TGCAT"}),
    
    PASI ("PasI", new String[] {"CAGGG", "CTGGG"}),
    
    HPAII ("HpaII", new String[] {"CGG"}),
    
    MSPI ("MspI", new String[]{"CGG"}),
    
    PST1_ECOT22I ("PstI-EcoT22I", new String[]{"TGCAG", "TGCAT"}),
    
    PSTI_MSPI ("PstI-MspI", new String[]{"TGCAG"}),
    
    PSTI_TAQI ("PstI-TaqI", new String[]{"TGCAG"}),
    
    SBFI_MSPI ("SbfI-MspI", new String[]{"TGCAGG"}),
    
    ASISI_MSPI ("AsiSI-MspI", new String[]{"ATCGC"}),
    
    BSSHII_MSPI ("BssHII-MspI", new String[]{"CGCGC"}),
    
    FSEI_MSPI ("FseI-MspI", new String[]{"CCGGCC"}),
    
    SALI_MSPI ("SalI-MspI", new String[]{"TCGAC"}),
    
    APOI ("ApoI", new String[]{"AATTC", "AATTT"}),
    
    BAMHI ("BamHI", new String[]{"GATCC"}),
    
    MSEI ("MseI", new String[]{"TAA"}),
    
    SAU3AI ("Sau3AI", new String[]{"GATC"}),
    
    RBSTA ("RBSTA", new String[]{"TA"}),
    
    RBSCG ("RBSCG", new String[]{"CG"}),
    
    NSPI ("NspI", new String[]{"CATGT", "CATGC"}),
    
    NAN ("NA", new String[]{""});
    
    private final String name;
    private final String[] initialCutSiteRemnant;
    
    /**
     * creates the enzyme
     * @param name String | the name of the enzyme
     * @param initialCutSite String[] | the cutsites
     */
    private EnzymeEnum(String name, String[] initialCutSite){
        this.name = name;
        this.initialCutSiteRemnant = initialCutSite;
    }
    
    /**
     * 
     * @return String the name of the enzyme
     */
    public String getName(){
        return this.name;
    }
        
    /**
     * returns all cutsites
     * @return Collection of String | all the cutsites of the enzyme in a collection
     */
    public Collection<String> getInitialCutSiteRemnant(){
        HashSet<String> cutsitesSet = new HashSet();
        cutsitesSet.addAll(Arrays.asList(this.initialCutSiteRemnant));
        return cutsitesSet;
    }
    
    /**
     * returns all complement cutsites
     * @return Collection of String | all the complement cutsites of the enzyme in a collection
     */
    public Collection<String> getComplementCutSiteRemnant() {
        HashSet<String> cutsitesSet = new HashSet();
        for (String cutsite : this.initialCutSiteRemnant){
            cutsitesSet.add(BasePair.getComplementSequence(cutsite));
        }
        return cutsitesSet;
    }
    
    
    /**
     * 
     * @param enzymeName | String the name of the enzyme
     * @return EnzymeEnum the corresponding EnzymeEnum or null if no enzyme is known
     */
    public static EnzymeEnum getEnzyme(String enzymeName){
        for (EnzymeEnum e : EnzymeEnum.values()){
            if (enzymeName.toLowerCase().equals(e.getName().toLowerCase())){
                return e;
            }
        }
        return null;
    }

    
}
