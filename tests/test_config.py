import pytest
from pysudoku_package.libraryexperiment import LibraryExperiment
from pysudoku_package.config import set_up_tmpdir, get_pool_dicts, get_input_files
import numpy as np
import os

def test_get_pool_dicts(ddir):
    
    
    file2pool_expected = {"test":"A", "test2":"B"}
    pool2dim_expected ={"row": ["A", "B"]}
    
    file2pool_dict, pool2dim_dict = get_pool_dicts(os.path.join(str(ddir),"exp_design.csv"))
    
    assert file2pool_expected == file2pool_dict
    assert pool2dim_expected == pool2dim_dict
    


def test_get_input_files(ddir):
    
    experiment = LibraryExperiment(8,30,10,os.path.join(str(ddir),"gb_ref/"), os.path.join(str(ddir),"bowtie_ref/UTI89"), 
                                    "AGATGTGTATAAGAGACAG", 1, str(ddir), os.path.join(str(ddir),"exp_design.csv"),
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr=0.05, 
                                    global_filter_thr = 5, min_counts=5)
    
    expected_paths = [os.path.join(str(ddir),'test.fastq'), os.path.join(str(ddir),'test2.fastq')]
    expected_names = ['test', 'test2']
    paths, names = get_input_files(experiment)
    assert expected_paths == paths
    assert expected_names == names
