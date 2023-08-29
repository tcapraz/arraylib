import pytest
import subprocess
import filecmp
import pandas as pd
import numpy as np

@pytest.mark.integtest
def test_integration():
    
    # run program from command line
    subprocess.call(['arraylib-run', 
                     'test_data/input',
                     'test_data/full_exp_design.csv', 
                     "-c",
                     "8", 
                     "-t",
                     " AGATGTGTATAAGAGACAG",
                     "-bu",
                     "CGAGGTCTCT",
                     "-bd",
                     "CGTACGCTGC",
                     "-br",
                     "test_data/bowtie_ref/UTI89",
                     "-gb",
                     "test_data/gb_ref/"])
    
    assert filecmp.cmp('count_matrix.csv', 'test_data/expected_count_matrix.csv')
    assert filecmp.cmp('filtered_count_matrix.csv', 'test_data/expected_filtered_count_matrix.csv')
    test = pd.read_csv('mutant_location_summary.csv')
    expected = pd.read_csv('test_data/expected_mutant_location_summary.csv')
    cols_to_test = ['Reference', 'Coordinate', 'Gene', 'Locus_tag', 'Predicted_wells',
           'Gene_orientation', 'Transposon_orientation', 'Percent_from_start',
           'Gene_hit', 'Ambiguity', 'Possible_wells', 'Barcode']
    for i in cols_to_test:
        assert np.all(test[i].astype(str) == expected[i].astype(str))
    
    #assert filecmp.cmp('mutant_location_summary.csv', 'test_data/expected_mutant_location_summary.csv')
    
    test = pd.read_csv('well_location_summary.csv')
    expected = pd.read_csv('test_data/expected_well_location_summary.csv')
    for i in test.columns:
        assert np.all(test[i].astype(str) == expected[i].astype(str))
    #assert filecmp.cmp('well_location_summary.csv', 'test_data/expected_well_location_summary.csv')

    
