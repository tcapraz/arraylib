import pytest
from pysudoku_package.libraryexperiment import LibraryExperiment
import numpy as np
import os
import string
import pandas as pd

@pytest.fixture(scope="session", name="experiment")
def get_example_data():
    experiment = LibraryExperiment(8,30,10,"test_data/gb_ref/", "test_data/bowtie_ref/UTI89", 
                                    "AGATGTGTATAAGAGACAG", 1, "test_data", "test_data/exp_design.csv",
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr=0.05,global_filter_thr=5, min_counts=5)
    
    numpr = 4
    numpc = 4
    rows = list(string.ascii_uppercase)[0:8] 
    cols = [str(i) for i in range(1,13)]
    prs = ["PR0"+str(i) if i < 10 else "PR" +str(i) for i in range(1,numpr+1)]
    pcs = ["PC0"+str(i) if i < 10 else "PC" +str(i) for i in range(1,numpc+1)]
    poolaxis=[rows,cols,prs, pcs]
    experiment.pools = rows + cols + prs+pcs
    experiment.pool_dims = {"row": rows, "col" :cols, "PR": prs, "PC": pcs}
    mat = np.zeros((10, len(experiment.pools)))
    
    mat[0:5,:] = 10
    mat[5:9, [0,10,23, 27]] = 10
    ref = np.array(["NC_007946.1" for i in range(10)])
    coord = np.array([i for i in range(230,240)])
    orient = np.array(["-" for i in range(10)])
    bar = np.array(["".join(np.array(np.random.choice( ["A", "C", "G", "T"],8, replace=True)).tolist()) for i in range(10)])
    data = pd.DataFrame(np.hstack((np.vstack((ref, coord, orient, bar)).T, mat.astype(int))), columns = ["Reference", "Feature", "Orientation", "Barcode"] + experiment.pools)
    
    experiment.filtered_count_mat = data
    experiment.count_mat = data

    return experiment