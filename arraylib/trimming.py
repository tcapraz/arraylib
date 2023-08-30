#!/usr/bin/env python3

import os
import numpy as np
from numba import njit
from pathlib import Path
import gzip 

from arraylib.io import txt_writer


@njit
def binary_subtract(array1,array2,mismatch):
    
    """ Used for matching 2 sequences based on the allowed mismatches.
    Requires the sequences to be in numerical form"""
    
    miss=0
    for arr1,arr2 in zip(array1,array2):
        if arr1-arr2 != 0:
            miss += 1
        if miss>mismatch:
            return 0
    return 1

@njit
def border_finder(seq,read,mismatch): 
    
    """ Matches 2 sequences (after converting to int8 format)
    based on the allowed mismatches. Used for sequencing searching
    a start/end place in a read"""
    
    s=seq.size
    r=read.size
    fall_over_index = r-s-1
    for i,bp in enumerate(read): #range doesnt exist in njit
        comparison = read[i:s+i]
        finder = binary_subtract(seq,comparison,mismatch)
        if i > fall_over_index:
            return
        if finder != 0:
            return i

    
def seq2bin(sequence):
    
    """ Converts a string to binary, and then to 
    a numpy array in int8 format"""
    
    byte_list = bytearray(sequence,'utf8')
    return np.array((byte_list), dtype=np.int8)

def barcode_extractor(border_up_bin,border_down_bin,read_bin):
    match_up=border_finder(border_up_bin,read_bin,1)
    match_down=border_finder(border_down_bin,read_bin,1)
    return match_up,match_down
    

def read_data(file, experiment):
    
    """ reads the fastq files. Searches for transposon border, and returns
    the subsequent trimmed sequence, read ID, and pool ID"""

    tn_seq_bin = seq2bin(experiment.tn_seq)
    
    if experiment.use_barcodes:
        border_up_bin = seq2bin(experiment.bar_upstream)
        border_down_bin = seq2bin(experiment.bar_downstream)
    
    reading,trimmed_seq = [],[]
    count = 0

    filename = Path(file).stem.split(".")[0]
    poolname = experiment.file2pool[filename]
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI"
    quality_set = set(quality_list[:int(experiment.seq_quality)-1])
    
    path_trimmed = os.path.join("temp", filename + "_trimmed_seq.fastq")
    with open(path_trimmed, "w") as output: # create empty file
        pass
    
    _, ext = os.path.splitext(file)
    if ext == ".gz":
        f = gzip.open(file, "rt") 
    else:
        f = open(file, "r") 
    for line in f:
        reading.append(line[:-1])

        if len(reading) == 4: #a read always has 4 lines
            count+=1
            ID = f"{count}:{poolname}:"
            read = reading[1]
            quality = reading[3]
            reading = []

            read_bin = seq2bin(read)
            match=border_finder(tn_seq_bin,read_bin,experiment.tn_mismatches)
            
            if match is not None:
                match+=len(tn_seq_bin)
                seq = read[match:]
                qual = quality[match:]
                # check if all bases have higher quality than param["phred"]
                if (len(quality_set.intersection(qual)) == 0):
                
                    if experiment.use_barcodes:
                        match_up,match_down=barcode_extractor(border_up_bin,border_down_bin,read_bin)
                        if (match_up is not None) & (match_down is not None):
                            barcode = read[match_up+len(border_up_bin):match_down]
                            ID = ID + f"{barcode}"
                        
                    trimmed_seq.append(f"@{ID}\n{seq}\n+\n{qual}\n")        
            
            #if psutil.virtual_memory().percent>=80:
            if len(trimmed_seq)>250000:
                txt_writer(path_trimmed, trimmed_seq)
                trimmed_seq = []
    f.close()
    txt_writer(path_trimmed, trimmed_seq)
    
    return path_trimmed


def detect_barcodes(file, experiment):
    
    """ reads the fastq files. Searches for transposon border, and returns
    barcodes read ID and pool ID"""

    
    if experiment.use_barcodes:
        border_up_bin = seq2bin(experiment.bar_upstream)
        border_down_bin = seq2bin(experiment.bar_downstream)
    
    reading,barcodes_id = [],[]
    count = 0

    filename = Path(file).stem.split(".")[0]
    poolname = experiment.file2pool[filename]
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI"
    quality_set = set(quality_list[:int(experiment.seq_quality)-1])
    
    path_temp = os.path.join("temp", filename + "_barcode_temp.csv")
    with open(path_temp, "w") as output: # create empty file
        pass
    
    _, ext = os.path.splitext(file)
    if ext == ".gz":
        f = gzip.open(file, "rt") 
    else:
        f = open(file, "r") 
    for line in f:
        reading.append(line[:-1])

        if len(reading) == 4: #a read always has 4 lines
            read = reading[1]
            quality = reading[3]
            reading = []

            read_bin = seq2bin(read)

            match_up,match_down=barcode_extractor(border_up_bin,border_down_bin,read_bin)
            if (match_up is not None) & (match_down is not None):
                count+=1

                barcode = read[match_up+len(border_up_bin):match_down]
                
                barcodes_id.append(f"{count},{poolname},0,nan,{barcode},nan\n")        
    
            #if psutil.virtual_memory().percent>=80:
            if len(barcodes_id)>250000:
                txt_writer(path_temp, barcodes_id)
                barcodes_id = []
    f.close()
    txt_writer(path_temp, barcodes_id)
    
    return path_temp
