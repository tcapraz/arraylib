#!/usr/bin/env python3

import pandas as pd
import numpy as np
from itertools import product
from Bio import SeqIO
import math
import functools
import os

def get_ambiguity(experiment):
    """
    Finds ambiguous and unambiguous mutants and returns them as separate arrays

    Parameters
    ----------
    experiment : LibraryExperiment
        Object holding input parameters.


    Returns
    -------
    unambiguous_data : numpy array
    ambiguous_data : numpy array

    """    
    

    # get only columns with read counts, so should also work with the pool_presence_table containing reads.
  
    dims = experiment.pool_dims
    data = experiment.filtered_count_mat
    
    pool_reads = data[experiment.pools].astype(int)
    
    nonzero_dimcounts = []
    for i in dims:
        nonzero_dimcounts.append(np.sum(pool_reads[dims[i]]>0, axis=1))
    
    dim_prod = functools.reduce(lambda x,y: x*y, nonzero_dimcounts)
    
    unpredictable  = dim_prod == 0
    
    nonzero_dimcounts = np.vstack(nonzero_dimcounts).T
    
    ambiguous =np.sum(nonzero_dimcounts > 1, axis=1) >=2
    ambiguous[unpredictable] = False
    ambiguous_data = data.loc[ambiguous,:]

    # we only allow one axis with more than 1 nonzero entry for unambiguous mutants
    unambiguous =np.sum(nonzero_dimcounts > 1, axis=1) <= 1
    unambiguous[unpredictable] = False
    unambiguous_data = data.loc[unambiguous,:]
    
    return unambiguous_data, ambiguous_data
    
def get_unambiguous_locations(data, experiment):
    """
    

    Parameters
    ----------
    data : pandas dataframe
        count matrix of unambiguous mutants
    experiment :  LibraryExperiment
        Object holding input parameters.


    Returns
    -------
    locations : pandas dataframe
        holding the well locations of unambiguous mutants

    """
    dims = experiment.pool_dims
    
    pool_reads = data[experiment.pools].astype(int)
    
    ref = data["Reference"]
    orientation = data["Orientation"]
    coords = data["Feature"]
    # get locations of unambiguous mutants in same format as predicted
    
    
    if "Barcode" in data.columns:
        outcols = ["Reference", "Feature", "Orientation", "Barcode",
                  "Predicted_wells", "Possible_wells", "Ambiguity"]
        barcode = True
        barcodes = data["Barcode"]
    else:
        barcode = False
        outcols = ["Reference", "Feature", "Orientation",
                  "Predicted_wells", "Possible_wells", "Ambiguity"]

    locations = pd.DataFrame(index=range(pool_reads.shape[0]), 
                                         columns=outcols)
    for i in range(pool_reads.shape[0]):
        line = pool_reads.iloc[i,:]
        nonzero_pools = []
        for j in dims:
            nonzero_pools.append(line[dims[j]][line[dims[j]]> 0].index.tolist())
        adds =  product(*nonzero_pools)

        wells = []
        for a in adds:
            wells.append("_".join(a))
        locations.loc[i,"Predicted_wells"] = ";".join(wells)
        locations.loc[i,"Possible_wells"] = ";".join(wells)
        locations.loc[i,"Orientation"] = orientation.values[i]
        locations.loc[i, "Feature"] =  coords.values[i]
        locations.loc[i, "Ambiguity"] =  "Unambiguous"
        locations.loc[i, "Reference"] =  ref.values[i]
        if barcode == True:
            locations.loc[i, "Barcode"] =  barcodes.values[i]
    
    return locations


def coords2genes(experiment, locations):
    """
    
    Infer genes from genomic coordinates using input genbank file.

    Parameters
    ----------
    experiment :  LibraryExperiment
        Object holding input parameters.
    locations : pandas dataframe
        holding the well locations the tn mutants


    Returns
    -------
    out : pandas dataframe
        output table holding mutant information and well locations

    """
    data = experiment.count_mat
    # reference genbank file
    references = np.unique(locations["Reference"].astype(str))
    
    recs = {i:[rec for rec in SeqIO.parse(os.path.join(experiment.gb_ref,i+ ".gb"), "genbank")] for i in references}
    # # - means complement strand -> reverse start end
    feats = {j:([feat for feat in recs[j][0].features if feat.type == "CDS"]) for j in recs}
    
    pools = experiment.pools
    
    outcols = ["Reference", "Coordinate", "Gene", "Locus_tag", 
                                               "Predicted_wells","Gene_orientation",
                                               "Transposon_orientation",
                                               "Percent_from_start", "Gene_hit", 
                                               "Ambiguity", "Possible_wells", "Pool_counts"]
    out = pd.DataFrame(index=range(data.shape[0]), columns = outcols)
    
    barcode = False
    if "Barcode" in locations.columns:
        barcode = True
        out["Barcode"] = data["Barcode"].reset_index(drop=True)


    # store nonzero pool counts in final summary table
    pool_reads = data[pools].astype(float)
    for i in range(pool_reads.shape[0]):
        line = pool_reads.iloc[i,:]
        line = line[line > 0]
        pool_counts = [p+":" + str(i) for p,i in zip(line.index, line)]   
        pool_counts = ";".join(pool_counts)    
        out.loc[i,"Pool_counts"] = pool_counts
        
    ref = data["Reference"]
    coords = data["Feature"]
    orientation = data["Orientation"]
    # populate with coords and orientation from pool presence table
    out["Coordinate"]  = coords.values
    out["Transposon_orientation"] = orientation.values
    out["Reference"] = ref.values

    # fill output dataframe, for now only with unambiguous mutants
    # this dataframe is equivalent to the progenitor_table of kosudoku, but ordered by mutants not wells
    for  i in range(locations.shape[0]):
        found = False
        r = locations["Reference"].values[i]
        c  = locations["Feature"].values[i]
        orientation_ = locations["Orientation"].values[i]
        # get idx of correct coordinate orientation pair
        if barcode == True:
            bar = locations["Barcode"].values[i]
            idx = np.where(np.logical_and(np.logical_and( out["Coordinate"] ==c , out["Transposon_orientation"] ==orientation_), out["Barcode"].astype(str)==str(bar)))[0]
        else:
            idx = np.where(np.logical_and( out["Coordinate"] ==c, out["Transposon_orientation"] ==orientation_))[0]
        c = int(c)
        for f in feats[r]:
            start = int(f.location.start)
            end = int(f.location.end)
            if (start < c < end):
                found = True
            if found == True:
                qual = f.qualifiers
                if "gene" in qual:
                    out.loc[idx,"Gene"] = ",".join(qual["gene"] )
                else:
                    out.loc[idx,"Gene"] =  "Not_found"
    
                if "locus_tag" in qual:
                    out.loc[idx, "Locus_tag"] = ",".join(qual["locus_tag"])
                else:
                    out.loc[idx,"Locus_tag"] =  "Not_found"
    
                if f.strand == 1:
                    length = end-start
                    length_from_start = c-start
                    percent_from_start = length_from_start/length
                    out.loc[idx, "Percent_from_start"] = percent_from_start
                    out.loc[idx, "Gene_orientation"] = "+"
    
                if f.strand == -1:
                    length = end-start
                    length_from_start = end - c
                    percent_from_start = length_from_start/length
                    out.loc[idx, "Percent_from_start"] = percent_from_start
                    out.loc[idx, "Gene_orientation"] = "-"
    
                break
        out.loc[idx, "Predicted_wells"] = locations["Predicted_wells"].values[i]
        out.loc[idx, "Possible_wells"] = locations["Possible_wells"].values[i]
        out.loc[idx, "Ambiguity"] = locations["Ambiguity"].values[i]
        out.loc[idx, "Transposon_orientation"] = orientation_
        if found: 
           out.loc[idx, "Gene_hit"] = True
        else:
           out.loc[idx, "Gene_hit"] = False 

        
    
    return out

def transpose_location_summary(experiment, summary):
    """
    Transpose mutant summary  to well summary where each row corresponds to 
    a different well.

    Parameters
    ----------
    experiment : LibraryExperiment
        Object holding input parameters.
    summary : pandas dataframe
        output table holding mutant information and well locations

    Returns
    -------
    out : pandas dataframe
       well summary

    """
    def append2cell(cell, val):
        if cell == "nan":
            out = val
        else:
            out = ";".join([str(cell),str(val)])
        return out

    
    dims = experiment.pool_dims
    
    well_pools = list(product(*dims.values()))
    
    wells = np.array(["_".join(w) for w in well_pools])
    wells = wells.reshape(len(well_pools),1)
    well_pools = np.array(well_pools)
    
    incols = ["Well"] +list(dims.keys()) + ["Reference", "Coordinate", "Gene", "Locus_tag", 
                                              "Gene_orientation", "Transposon_orientation",
                                               "Percent_from_start", "Gene_hit", 
                                               "Ambiguity"]
    outcols = ["Reference", "Coordinate", "Gene", "Locus_tag", "Gene_orientation", 
               "Transposon_orientation","Percent_from_start", "Gene_hit", "Ambiguity"]
    barcode = False
    if "Barcode" in summary.columns:
        barcode = True
        dummy = np.empty((well_pools.shape[0],10))
        dummy[:] = np.nan
        incols = incols + ["Barcode"]
        outcols = outcols + ["Barcode"]
    else:
        dummy = np.empty((well_pools.shape[0],9))
        dummy[:] = np.nan
    grid = np.hstack([wells, well_pools, dummy])

    out = pd.DataFrame(grid, columns= incols)
    out.index = out["Well"]

    for row  in summary.itertuples(index =False):
        w = row[4]
        if type(w) == float and  math.isnan(w):
            continue
        else:
            w = w.split(";")
            r = row[0]
            c = row[1]
            g = row[2]
            l = row[3]
            gorient = row[5]
            torient = row[6]
            pfromstart = row[7]
            ghit = row[8]
            amb = row[9]
            if barcode==True:
                bar = row[12]
                outvals = [r,c,g,l,gorient,torient,pfromstart,ghit,amb,bar]
            else:
                outvals = [r,c,g,l,gorient,torient,pfromstart,ghit,amb]
            for i in w:
                for j,v in zip(outcols,outvals):
                    out.loc[i,j] = append2cell(out.loc[i,j], v)
    
    return out
    
def transpose_barcode_summary(experiment, summary):
    """
    Transpose barcode summary  to well summary where each row corresponds to 
    a different well.

    Parameters
    ----------
    experiment : LibraryExperiment
        Object holding input parameters.
    summary : pandas dataframe
        output table holding barcode information and well locations

    Returns
    -------
    out : pandas dataframe
       well summary

    """
    def append2cell(cell, val):
        if cell == "nan":
            out = val
        else:
            out = ";".join([str(cell),str(val)])
        return out

    
    dims = experiment.pool_dims
    
    well_pools = list(product(*dims.values()))
    
    wells = np.array(["_".join(w) for w in well_pools])
    wells = wells.reshape(len(well_pools),1)
    well_pools = np.array(well_pools)
    
    incols = ["Well"] +list(dims.keys()) + ["Ambiguity", "Barcode"]
    outcols = ["Ambiguity", "Barcode"]

    dummy = np.empty((well_pools.shape[0],2))
    dummy[:] = np.nan

    grid = np.hstack([wells, well_pools, dummy])

    out = pd.DataFrame(grid, columns= incols)
    out.index = out["Well"]

    for row  in summary.itertuples(index =False):
        w = row[1]
        if type(w) == float and  math.isnan(w):
            continue
        else:
            w = w.split(";")
           
            amb = row[3]
            
            bar = row[0]
            outvals = [amb,bar]
           
            for i in w:
                for j,v in zip(outcols,outvals):
                    out.loc[i,j] = append2cell(out.loc[i,j], v)
    
    return out
    
