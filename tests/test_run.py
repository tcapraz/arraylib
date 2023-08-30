import pytest
import subprocess
import filecmp
import pandas as pd
import numpy as np
import os

@pytest.mark.integtest
def test_integration(ddir):
    
    # run program from command line
    subprocess.call(['arraylib-run', 
                     'tests/test_data/input',
                     'tests/test_data/full_exp_design.csv', 
                     "-c",
                     "8", 
                     "-t",
                     " AGATGTGTATAAGAGACAG",
                     "-bu",
                     "CGAGGTCTCT",
                     "-bd",
                     "CGTACGCTGC",
                     "-br",
                     "tests/test_data/bowtie_ref/UTI89",
                     "-gb",
                     "tests/test_data/gb_ref/"])
    
    #assert filecmp.cmp('count_matrix.csv', os.path.join(str(ddir),'expected_count_matrix.csv'))
    assert filecmp.cmp('count_matrix.csv', os.path.join(str(ddir),'expected_filtered_count_matrix.csv'))
    test = pd.read_csv('mutant_location_summary.csv')
    expected = pd.read_csv(os.path.join(str(ddir),'expected_mutant_location_summary.csv'))
    cols_to_test = ['Reference', 'Coordinate', 'Gene', 'Locus_tag', 'Predicted_wells',
           'Gene_orientation', 'Transposon_orientation', 'Percent_from_start',
           'Gene_hit', 'Ambiguity', 'Possible_wells', 'Barcode']
    for i in cols_to_test:
        assert np.all(test[i].astype(str) == expected[i].astype(str))
    
    #assert filecmp.cmp('mutant_location_summary.csv', 'test_data/expected_mutant_location_summary.csv')
    
    test = pd.read_csv('well_location_summary.csv')
    expected = pd.read_csv(os.path.join(str(ddir),'expected_well_location_summary.csv'))
    for i in test.columns:
        assert np.all(test[i].astype(str) == expected[i].astype(str))
    #assert filecmp.cmp('well_location_summary.csv', 'test_data/expected_well_location_summary.csv')

    
