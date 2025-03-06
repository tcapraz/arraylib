#!/usr/bin/env python3

import os
from pathlib import Path
import gzip 
from fast2q.fast2q import border_finder,seq2bin,sequence_tinder

from arraylib.io import txt_writer

def file_prep(file, experiment):

    """ Converts the barcode search sequences into its numeric form before 
    barcode_finder uses them for barcode extraction.    
    Saved as a list in a dictionary value"""

    file_name_extended = "_barcode_temp.csv"

    if experiment.use_barcodes:
        file_name_extended = "_trimmed_seq.fastq"

    filename = Path(file).stem.split(".")[0]
    poolname = experiment.file2pool[filename]

    path_trimmed = os.path.join("temp", filename + file_name_extended)
    with open(path_trimmed, "w") as output: # create empty file
        pass

    _, ext = os.path.splitext(file)
    if ext == ".gz":
        f = gzip.open(file, "rt") 
    else:
        f = open(file, "r") 

    barcode_info = {'upstream':experiment.bar_upstream,
                    'downstream':experiment.bar_downstream,
                    'upstream_bin':[seq2bin(experiment.bar_upstream)],
                    'downstream_bin':[seq2bin(experiment.bar_downstream)],
                    'miss_search_up':1,
                    'miss_search_down':1,
                    'quality_set_up':set(""),
                    'quality_set_down':set("")}

    return barcode_info,path_trimmed,poolname,f

def read_data(file, experiment):
    
    """ reads the fastq files. Searches for transposon border, and returns
    the subsequent trimmed sequence, read ID, and pool ID"""

    tn_seq_bin = seq2bin(experiment.tn_seq)
    barcode_info,path_trimmed,poolname,f = file_prep(file, experiment)

    reading,trimmed_seq = [],[]
    count = 0

    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI"
    quality_set = set(quality_list[:int(experiment.seq_quality)-1])

    for line in f:
        reading.append(line.rstrip())

        if len(reading) == 4:
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

                if (len(quality_set.intersection(qual)) == 0):
                
                    if experiment.use_barcodes:
                        start,end=sequence_tinder(read_bin,
                                                    quality.encode("utf-8"),
                                                    barcode_info)
                    
                        if (start is not None) & (end is not None):
                            barcode = read[start:end]
                            ID = ID + barcode
                        
                    trimmed_seq.append(f"@{ID}\n{seq}\n+\n{qual}\n")        
            
            if len(trimmed_seq)>250000:
                txt_writer(path_trimmed, trimmed_seq)
                trimmed_seq = []
    f.close()
    txt_writer(path_trimmed, trimmed_seq)
    
    return path_trimmed


def detect_barcodes(file, experiment):
    
    """ reads the fastq files. Searches for transposon border, and returns
    barcodes read ID and pool ID"""
    
    barcode_info,path_temp,poolname,f = file_prep(file, experiment)
    
    reading,barcodes_id = [],[]
    count = 0

    for line in f:
        reading.append(line.rstrip())

        if len(reading) == 4:
            read = reading[1]
            quality = reading[3]
            reading = []

            start,end=sequence_tinder(seq2bin(read),
                                        quality.encode("utf-8"),
                                        barcode_info)

            if (start is not None) & (end is not None):
                barcode = read[start:end]
                count+=1
                barcodes_id.append(f"{count},{poolname},0,nan,{barcode},nan\n")        
    
            if len(barcodes_id)>250000:
                txt_writer(path_temp, barcodes_id)
                barcodes_id = []
    f.close()
    txt_writer(path_temp, barcodes_id)
    
    return path_temp
