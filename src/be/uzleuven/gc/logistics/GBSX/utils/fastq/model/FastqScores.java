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
package be.uzleuven.gc.logistics.GBSX.utils.fastq.model;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public enum FastqScores {
    /**
     * The quality of Illumina1.8
     * <br> Type = Illumina1.8
     * <br> min Score = 0
     * <br> max Score = 41
     * <br> start score = 33
     */
    ILLUMINA_1_8 ("Illumina1.8", 0, 41, 33),
    /**
     * The quality of Illumina1.3
     * <br> Type = Illumina1.3
     * <br> min Score = 3
     * <br> max Score = 40
     * <br> start score = 67
     */
    ILLUMINA_1_5 ("Illumina1.5", 3, 40, 67),
    /**
     * The quality of Illumina1.3
     * <br> Type = Illumina1.3
     * <br> min Score = 0
     * <br> max Score = 40
     * <br> start score = 64
     */
    ILLUMINA_1_3 ("Illumina1.3", 0, 40, 64),
    /**
     * The quality of Sanger
     * <br> Type = Sanger
     * <br> min Score = -5
     * <br> max Score = 40
     * <br> start score = 33
     */
    SANGER ("Sanger", -5, 40, 33),
    /**
     * The quality of Solid
     * <br> Type = Solid
     * <br> min Score = 0
     * <br> max Score = 40
     * <br> start score = 59
     */
    SOLID ("Solid", 0, 40, 59);
    
    private final int minScore;
    private final int maxScore;
    private final int startScore;
    private final String type;
    
    /**
     * returns the FastqScores of the given type, or the standard if the type doesn't exists
     * @param type String | the name of the type
     * @return the FastqScore of the type, or the standard
     * @see FastqScores#getStandard() 
     */
    public static FastqScores getFastqScores(String type){
        for (FastqScores score : FastqScores.values()){
            if (type.equals(score.getType())){
                return score;
            }
        }
        return FastqScores.getStandard();
    }
    
    /**
     * returns the standard FastqScores: Illumina1.8
     * @return FastqScores.ILLUMINA_1_8
     * @see FastqScores#ILLUMINA_1_8
     */
    public static FastqScores getStandard(){
        return FastqScores.ILLUMINA_1_8;
    }
    
    private FastqScores(String type, int minScore, int maxScore, int startScore){
        this.type = type;
        this.minScore = minScore;
        this.maxScore = maxScore;
        this.startScore = startScore;
    }
    
    /**
     * returns the type of the quality:
     * Illumina1.8, Illumina1.5, Illumina1.3, Sanger or Solid
     * @return 
     */
    public String getType(){
        return this.type;
    }
    
    /**
     * returns the minimum value of the score
     * @return int | minimum score
     */
    public int getMinScore(){
        return this.minScore;
    }
    
    /**
     * returns the maximum value of the score
     * @return int | maximum score
     */
    public int getMaxScore(){
        return this.maxScore;
    }
    
    /**
     * return the value of the startscore. If this value is substracted from the ASCII code of the quality, 
     * <br> the real score is get (also add the minscore if this is not 0)
     * 
     * @return int | the start score
     */
    public int getStartScore(){
        return this.startScore;
    }
    
    /**
     * 
     * @param baseQuality char | the basequality in ASCII format
     * @return int | the actual phred score
     */
    public int getPhredScoreOfBase(char baseQuality){
        return ((int) baseQuality - this.startScore + this.minScore);
    }
    
    /**
     * 
     * @param baseQuality int | the phred score of a base
     * @return char | the ASCII code of that phred score
     */
    public char getCodedScoreOfBase(int baseQuality){
        return Character.toChars(baseQuality - this.minScore + this.startScore)[0];
    }
    
}
