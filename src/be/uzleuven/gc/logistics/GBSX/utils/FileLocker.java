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
package be.uzleuven.gc.logistics.GBSX.utils;

import be.uzleuven.gc.logistics.GBSX.utils.fastq.infrastructure.FastqBufferedReader;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Koen Herten for the KU Leuven
 */
public class FileLocker {
    /**
     * LOCK PROCEDURES
     */
    
    private boolean isLocked = false;
    private boolean stopLock = false;
    private int waiting = 0;
    
    private int locked = 0;
    
    /**
     * will wait until the writer is unlocked, then lock it
     * also reserves space to lock (isCompleteUnlocked == false untill unlock)
     * @return true when locked, false when lock is impossible (the reader is closed)
     */
    public boolean lock(){
        locked++;
        this.waiting++;
        while (this.isLocked){
            if (this.stopLock){
                return false;
            }
            try {
                Thread.sleep(1);
            } catch (InterruptedException ex) {
                Logger.getLogger(FastqBufferedReader.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        this.isLocked = true;
        return this.isLocked;
    }
    
    /**
     * will unlock the writer
     * also give the space back
     * @return false
     */
    public boolean unlock(){
        this.waiting--;
        this.isLocked = false;
        return this.isLocked;
    }
    
    /**
     * checks if there are more processes waiting
     * @return true if no process is waiting for a lock
     */
    public boolean isCompleteUnlocked(){
        if (this.waiting == 0){
            return true;
        }else{
            return false;
        }
    }
    
    /**
     * waits till the writer is completely unlocked, then lock the writer (no unlock possible)
     * @return true
     */
    public boolean waitTillCompleteUnlock(){
        while (! this.isCompleteUnlocked()){
            try {
                Thread.sleep(10);
            } catch (InterruptedException ex) {
                Logger.getLogger(FastqBufferedReader.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        this.isLocked = true;
        this.stopLock = true;
        return true;
    }
    
    /**
     * END LOCK PROCEDURES
     */
}
