import pytest
from fast2q.fast2q import seq2bin, binary_subtract, border_finder, sequence_tinder
from arraylib.trimming import read_data
from arraylib.libraryexperiment import LibraryExperiment
from arraylib.config import set_up_tmpdir
import numpy as np
import os
def test_seq2bin():
    
    
    seq = "ACTGC"
    binseq = seq2bin(seq)
    
    assert len(seq) == len(binseq)
    assert type(binseq) == np.ndarray

def test_binary_substract():
    
    arr1 = seq2bin("AAAAAA")
    arr2 = seq2bin("AAAAAA")    
    arr3 = seq2bin("AAAAAT")
    
    assert binary_subtract(arr1, arr2, 0) == 1
    assert binary_subtract(arr1, arr2, 2) == 1
    assert binary_subtract(arr1, arr3, 0) == 0
    assert binary_subtract(arr1, arr3, 1) == 1

def test_border_finder():
    
    seq = seq2bin("AAACTGCAAAAAAAAA")
    bar1 = seq2bin("ACTGC")
    bar2 = seq2bin("AAA")
    bar3 = seq2bin("ACCGC")
    
    assert border_finder(bar1, seq, 0) == 2
    assert border_finder(bar2, seq, 0) == 0
    assert border_finder(bar3, seq, 0) == None
    assert border_finder(bar3, seq, 1) == 2


def test_sequence_tinder():
    
    seq = "AAACTGCAAAACCC\n"
    up = "AAA"
    down = "ACCC"

    quality = "IIIIIIIIIIIIIIIIIIIIIIIIIIIII".encode("utf-8")

    barcode_info = {'upstream':up,
                    'downstream':down,
                    'upstream_bin':[seq2bin(up)],
                    'downstream_bin':[seq2bin(down)],
                    'miss_search_up':0,
                    'miss_search_down':0,
                    'quality_set_up':set(""),
                    'quality_set_down':set("")}
    
    barcode = "CTGCAAA"

    start,end = sequence_tinder(seq2bin(seq), quality, barcode_info)
    extracted = seq[start:end]

    assert extracted == barcode


def test_read_data(ddir):
    experiment = LibraryExperiment(8,30,10,os.path.join(str(ddir),"gb_ref/"), os.path.join(str(ddir),"bowtie_ref/UTI89"), 
                                    "AGATGTGTATAAGAGACAG", 1, str(ddir), os.path.join(str(ddir),"exp_design.csv"),
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr = 0.05,
                                    global_filter_thr = 5, min_counts=5)
    truecontents = open(os.path.join(str(ddir),"trimmed_expected.fastq")).read()
    set_up_tmpdir()
    try:
        outfile_path = read_data(os.path.join(str(ddir),"test.fastq"), experiment)
        contents = open(outfile_path).read()
    finally:
            os.remove(outfile_path)
    assert outfile_path  == os.path.join("temp","test_trimmed_seq.fastq")
    assert contents == truecontents

    
