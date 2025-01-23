import pytest
import pandas as pd
import numpy as np
import os
from arraylib.simulations import simulate_required_arraysize, simulate_unique_genes, plot_arraysize_vs_unique_genes


def test_simulate_required_arraysize(ddir):
    
    data = os.path.join(str(ddir),"tnseeker_test_output.csv")
    
    res, plot = simulate_required_arraysize(data, 
                                           minsize=0,
                                           maxsize=500,
                                           number_of_simulations=5, 
                                           number_of_repeats=3, 
                                           gene_start = 0.1, 
                                           gene_end = 0.9)
    assert type(res) == pd.DataFrame
    assert res.shape[0] == 5
    assert ~np.any(res.isna())
    assert plot is not None
    
def test_simulate_unique_genes(ddir):
    
    data = pd.read_csv(os.path.join(str(ddir),"tnseeker_test_output.csv"))
    arrsize=np.linspace(0,500, 5)

    
    res = simulate_unique_genes(data,
                            arrsize,
                            number_of_repeats=3, 
                            gene_start = 0.1, 
                            gene_end = 0.9)
    assert type(res) == pd.DataFrame
    assert res.shape[0] == 5
    assert ~np.any(res.isna())

def test_plot_arraysize_vs_unique_genes(ddir):
    data = pd.read_csv(os.path.join(str(ddir),"tnseeker_test_output.csv"))
    arrsize=np.linspace(0,500, 5)

    
    res = simulate_unique_genes(data,
                            arrsize,
                            number_of_repeats=3, 
                            gene_start = 0.1, 
                            gene_end = 0.9)
    plot = plot_arraysize_vs_unique_genes(res)
    assert plot is not None