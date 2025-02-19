#!/usr/bin/env python3

import os
from pathlib import Path
import gzip 
from fast2q.fast2q import border_finder,seq2bin
from tnseeker.reads_trimer import barcode_finder

from arraylib.io import txt_writer

def barcode_converter(file, experiment):

    """ Converts the barcode search sequences into its numeric form before 
    barcode_finder uses them for barcode extraction.    
    Saved as a list in a dictionary value"""

    variables = {}
    file_name_extended = "_barcode_temp.csv"

    if experiment.use_barcodes:
        file_name_extended = "_trimmed_seq.fastq"
        variables["borders_bin"] = [seq2bin(experiment.bar_upstream),seq2bin(experiment.bar_downstream)]

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

    return variables,path_trimmed,poolname,f

def read_data(file, experiment):
    
    """ reads the fastq files. Searches for transposon border, and returns
    the subsequent trimmed sequence, read ID, and pool ID"""

    tn_seq_bin = seq2bin(experiment.tn_seq)
    variables,path_trimmed,poolname,f = barcode_converter(file, experiment)

    reading,trimmed_seq = [],[]
    count = 0

    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI"
    quality_set = set(quality_list[:int(experiment.seq_quality)-1])

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

                if (len(quality_set.intersection(qual)) == 0):
                
                    if experiment.use_barcodes:
                        barcode=barcode_finder(read_bin,variables)
                        if barcode is not None:
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
    
    variables,path_temp,poolname,f = barcode_converter(file, experiment)
    
    reading,barcodes_id = [],[]
    count = 0

    for line in f:
        reading.append(line[:-1])

        if len(reading) == 4: #a read always has 4 lines
            read = reading[1]
            reading = []

            barcode=barcode_finder(seq2bin(read),variables)
            if barcode is not None:
                count+=1
                barcodes_id.append(f"{count},{poolname},0,nan,{barcode},nan\n")        
    
            if len(barcodes_id)>250000:
                txt_writer(path_temp, barcodes_id)
                barcodes_id = []
    f.close()
    txt_writer(path_temp, barcodes_id)
    
    return path_temp
