import pytest
from pysudoku_package.trimming import seq2bin, binary_subtract, border_finder, barcode_extractor, read_data
from pysudoku_package.libraryexperiment import LibraryExperiment
from pysudoku_package.config import set_up_tmpdir
from pysudoku_package.generate_count_matrix import get_hash_ids, table_lookup, get_pool_table_line, get_count_matrix, filter_barcodes
import numpy as np
import os
import pandas as pd


def test_get_hash_ids():
    
    ids = ["the", "big", "brown", "fox", "jumps"]
    encoded_ids = get_hash_ids(ids)
    
    assert len(np.unique(ids)) == len(np.unique(encoded_ids))
    
    assert type(encoded_ids[0]) == np.int64
    assert type(encoded_ids) == np.ndarray


def test_table_lookup():
    ids = get_hash_ids(["the", "big", "the", "brown", "fox", "jumps"])
    
    assert all(table_lookup(ids[0], ids)    == np.array([0,2]))
    assert all(table_lookup(ids[1], ids)    == np.array(1))


def test_get_pool_table_line(ddir):
    bowtie_res =   pd.read_csv(os.path.join(str(ddir),"alignment_result_expected.csv"),
                               names=["id_in_pool","pool","coord","orientation","barcode", "reference"])

    pools =[3]
    
    # make sure pools are always on same position
    poolidx_dict = {p : i for i,p in enumerate(pools)}
    
    bowtie_res["barcode"] = bowtie_res["barcode"].astype(str)

    # barcodes that were not found are treated the same (can't distinguish between them)
    # concat orientation, coordinate and barcode to create unique ids
    bowtie_res["unique_id"] = bowtie_res["orientation"].astype(str) + ";"+ \
        bowtie_res["coord"].astype(str) + ";"  + bowtie_res["barcode"].astype(str)+ ";" + bowtie_res["reference"]
    
    all_ids = get_hash_ids(bowtie_res["unique_id"])
    count, id_ = get_pool_table_line(all_ids[0], bowtie_res.values, all_ids, pools, poolidx_dict)
    
    assert count == 6
    assert id_ == bowtie_res["unique_id"][0]


def test_get_count_matrix(ddir):
    experiment = LibraryExperiment(8,30,10,"gb_ref/", "bowtie_ref/UTI89", 
                                    "AGATGTGTATAAGAGACAG", 1, "test_data", os.path.join(str(ddir),"exp_design.csv"),
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr=0.05, 
                                    global_filter_thr = 5, min_counts=5)
    bowtie_res =   pd.read_csv(os.path.join(str(ddir),"alignment_result_expected.csv"),
                               names=["id_in_pool","pool","coord","orientation","barcode", "reference"])

    experiment.pools =[3]
    
    # make sure pools are always on same position
    
    bowtie_res["barcode"] = bowtie_res["barcode"].astype(str)

    # barcodes that were not found are treated the same (can't distinguish between them)
    # concat orientation, coordinate and barcode to create unique ids
    bowtie_res["unique_id"] = bowtie_res["orientation"].astype(str) + ";"+ \
        bowtie_res["coord"].astype(str) + ";"  + bowtie_res["barcode"].astype(str)+ ";" + bowtie_res["reference"]
    all_ids = get_hash_ids(bowtie_res["unique_id"])
    
    count_mat = get_count_matrix(all_ids, bowtie_res, experiment).iloc[0,:]
    expected = pd.Series([ 61, "+", "ATGACTTGCGGTAAAGAAGA","NC_007946.1", "6"], index = ["Feature", "Orientation" ,"Barcode" ,"Reference", 3])
    assert all(count_mat == expected)

def test_filter_barcodes(ddir):
    experiment = LibraryExperiment(8,30,10,"gb_ref/", "bowtie_ref/UTI89", 
                                    "AGATGTGTATAAGAGACAG", 1, "test_data", os.path.join(str(ddir),"exp_design.csv"),
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr=0.05, 
                                    global_filter_thr = 5, min_counts=5)
    
    mat = np.zeros((10,10))
    add1 = np.zeros(10)
    add1[0] = 10

    mat = mat + add1
    mat = mat.T
    mat[0,0] = mat[0,0] -9
    
    ref = np.array(["NC_007946.1"for i in range(10) ])
    feature = np.array(["100" for i in range(10)])
    orient = np.array(["+" for i in range(10)])
    cols = [str(i) for i in range(10)]
    experiment.pools = cols
    
    data = pd.DataFrame(np.vstack((ref, feature, orient, mat.T.astype(int))).T, columns = ["Reference", "Feature", "Orientation"]+cols)
    expected = pd.DataFrame(np.vstack((ref, feature, orient, mat.T.astype(int))).T, columns = ["Reference", "Feature", "Orientation"]+cols).iloc[0,:]

    filtered = filter_barcodes(data, experiment)
    
    assert all(expected == filtered)
    
    mat = np.zeros((10,10))
    add1 = np.zeros(10)
    add1[0] = 10

    mat = mat + add1
    mat = mat.T
    data = pd.DataFrame(np.vstack((ref, feature, orient, mat.astype(int))).T, columns = ["Reference", "Feature", "Orientation"]+cols)

    filtered = filter_barcodes(data, experiment)
    assert all(data == filtered)
