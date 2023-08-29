import pytest
from pysudoku_package.libraryexperiment import LibraryExperiment
import numpy as np
import os
from pysudoku_package.deconvolution import get_ambiguity, get_unambiguous_locations, coords2genes,transpose_location_summary
import string
import pandas as pd
from pandas.testing import assert_frame_equal
from itertools import product

def test_get_ambiguity(experiment):
    amb_expected = experiment.count_mat.iloc[0:5,:]
    unamb_expected = experiment.count_mat.iloc[5:9,:]
    unamb,amb = get_ambiguity(experiment)
    
    assert all(amb == amb_expected)
    assert all(unamb == unamb_expected)


def test_get_unambiguous_locations(experiment):
    unamb,amb = get_ambiguity(experiment)
    
    ref = unamb["Reference"]

    feature = unamb["Feature"]
    orient = unamb["Orientation"]
    bar = unamb["Barcode"]
    
    ambiguity = np.array(["Unambiguous" for i in range(len(bar))])
    pred = np.array(["A_3_PR04_PC04" for i in range(len(bar))])
    possible = pred
    expected = pd.DataFrame(np.vstack((ref, feature, orient, bar, pred, possible, ambiguity)).T, 
                 columns = ["Reference", "Feature","Orientation", "Barcode", "Predicted_wells", "Possible_wells" ,  "Ambiguity"])
    
    locs = get_unambiguous_locations(unamb, experiment)
    assert np.all(locs == expected)

def test_coords2genes(experiment):
    unamb,amb = get_ambiguity(experiment)
    locs = get_unambiguous_locations(unamb, experiment)
    
    data = experiment.count_mat
    outcols = ["Reference", "Coordinate", "Gene", "Locus_tag", 
                                               "Predicted_wells","Gene_orientation",
                                               "Transposon_orientation",
                                               "Percent_from_start", "Gene_hit", 
                                               "Ambiguity", "Possible_wells", "Pool_counts", "Barcode"]
    out = pd.DataFrame(index=range(data.shape[0]), columns = outcols)
    
    ref = data["Reference"]
    coords = data["Feature"]
    orientation = data["Orientation"]
    out["Reference"]  = ref.values
    out["Coordinate"]  = coords.values
    out["Transposon_orientation"] = orientation.values
    out.loc[unamb.index,"Gene"] = "thrL"
    out.loc[unamb.index,"Locus_tag"] ="UTI89_RS00005"
    out.loc[unamb.index,"Gene_orientation"] = "+"
    out.loc[unamb.index,"Predicted_wells"] = locs["Predicted_wells"].values
    out.loc[unamb.index,"Possible_wells"] = locs["Possible_wells"].values
    out.loc[unamb.index,"Gene_hit"] = True
    out.loc[unamb.index,"Ambiguity"] = locs["Ambiguity"].values
    out["Barcode"] = data["Barcode"]
    pool_reads = data[experiment.pools].astype(float)
    for i in range(pool_reads.shape[0]):
        line = pool_reads.iloc[i,:]
        line = line[line > 0]
        pool_counts = [p+":" + str(i) for p,i in zip(line.index, line)]   
        pool_counts = ";".join(pool_counts      )    
        out.loc[i,"Pool_counts"] = pool_counts
    
    start = 189
    end = 255
    length = end-start
    percent_from_start = [(c-start)/length for c in locs["Feature"].astype(int)]
    out.loc[unamb.index,"Percent_from_start"] = percent_from_start
    
    summary = coords2genes(experiment, locs)
    assert_frame_equal(summary, out)

def test_transpose_location_summary(experiment):
    unamb,amb = get_ambiguity(experiment)
    locs = get_unambiguous_locations(unamb, experiment)
    summary = coords2genes(experiment, locs)
    start = 189
    end = 255
    length = end-start
    percent_from_start = [str((c-start)/length) for c in locs["Feature"].astype(int)]
    transposed = transpose_location_summary(experiment,summary)
    
    
    
    dims = experiment.pool_dims
    
    well_pools = list(product(*dims.values()))
    
    wells = np.array(["_".join(w) for w in well_pools])
    wells = wells.reshape(len(well_pools),1)
    well_pools = np.array(well_pools)
    
    incols = ["Well"] +list(dims.keys()) + ["Reference", "Coordinate", "Gene", "Locus_tag", 
                                              "Gene_orientation", "Transposon_orientation",
                                               "Percent_from_start", "Gene_hit", 
                                               "Ambiguity", "Barcode"]
    outcols = ["Reference", "Coordinate", "Gene", "Locus_tag", "Gene_orientation", 
               "Transposon_orientation","Percent_from_start", "Gene_hit", "Ambiguity", "Barcode"]
    
    dummy = np.empty((well_pools.shape[0],10))
    dummy[:] = np.nan
    grid = np.hstack([wells, well_pools, dummy])

    out = pd.DataFrame(grid, columns= incols)
    out.index = out["Well"]
    
    idx = locs["Predicted_wells"][0]
    out.loc[idx, "Reference"] = "NC_007946.1;NC_007946.1;NC_007946.1;NC_007946.1"
    out.loc[idx, "Coordinate"] = "235;236;237;238"
    out.loc[idx, "Locus_tag"] = "UTI89_RS00005;UTI89_RS00005;UTI89_RS00005;UTI89_RS00005"
    out.loc[idx, "Gene"] = "thrL;thrL;thrL;thrL"
    out.loc[idx, "Gene_orientation"] = "+;+;+;+"
    out.loc[idx, "Transposon_orientation"] = "-;-;-;-"
    out.loc[idx, "Percent_from_start"] = ";".join(percent_from_start)
    out.loc[idx, "Ambiguity"] = "Unambiguous;Unambiguous;Unambiguous;Unambiguous"
    out.loc[idx, "Barcode"] = ";".join(locs["Barcode"].tolist())

    out.loc[idx, "Gene_hit"] = "True;True;True;True"
    
    assert_frame_equal(transposed, out)
