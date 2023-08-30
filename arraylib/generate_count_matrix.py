#!/usr/bin/env python3

import multiprocessing as mp
import pandas as pd
import numpy as np
from itertools import repeat
from numba import njit
import hashlib
import re
import scipy.stats as stats

@njit
def table_lookup(coord,data):
    """
    njit compiled np.where to speed up computation

    Parameters
    ----------
    coord : encoding genomic coordinate, barcode and 
    transposon orientation
    data : numpy array
        all ids

    Returns
    -------
    numpy array
        indices of coord in data

    """
    return np.where(data==coord)[0]



def get_filter_mask(id_, all_ids, data):
    """
    Get filter mask to remove all barcodes that have count sums below 
    10 % of the max barcode.

    Parameters
    ----------
    id_ : encoding genomic coordinate, barcode and 
    transposon orientation
    all_ids : numpy array
        all hashed ids encoding genomic coordinate, barcode and 
        transposon orientation
    data : numpy array
        count matrix

    Returns
    -------
    boolean array
        holding the filter mask

    """
    sub = data[table_lookup(id_,all_ids),:]
    
    total = np.sum(sub)
    rowsums = np.sum(sub, axis=1)
    row_fracs = rowsums/total

    # cutoff = stats.norm.ppf(0.0001, loc=np.mean(rowsums), scale=np.mean(rowsums)/4)
    # return rowsums > cutoff
    return row_fracs >= 0.1

def filter_barcodes(data,experiment):
    """
    Filter out barcodes with little read counts.

    Parameters
    ----------
    data : pandas dataframe
        count matrix
    experiment : LibraryExperiment
        Object holding input parameters.

    Returns
    -------
    out : pandas dataframe
        Filtered count matrix

    """
    out= data.copy()
    pools = experiment.pools
    ids = out["Feature"].astype(str) + out["Orientation"] + out["Reference"]
    hashed_ids = get_hash_ids(ids) 
    out["id"] = hashed_ids
    assert len(np.unique(ids)) == len(np.unique(out["id"])), "Hashed ids not unique!"

    unique_ids = out["id"].unique()
        
    pool = mp.Pool(experiment.cores)
    mask = pool.starmap_async(get_filter_mask, zip(unique_ids, repeat(hashed_ids),
                              repeat(out[pools].values.astype(int)))).get()
    pool.close()
    
    mask = np.hstack(mask)
    out = out[mask]
    out = out.drop("id", axis=1)
    return out

def local_filter_counts(counts, thr=0.05):
    """
    Apply local relative filter. Counts are set to zero if their ratio to the 
    max count is < thr.

    Parameters
    ----------
    counts : numpy array
        count matrix
    thr : int, optional
        relative filter threshold. The default is 0.05.

    Returns
    -------
    counts : numpy array
        filtered count matrix

    """
    m = np.max(counts, axis=1)
    mask = (counts.T/m).T < thr
    counts[mask] = 0
    
    return counts

def global_filter_counts(counts, thr=5):
    """
    Apply global absolut filter. Counts are set to zero if they are < thr

    Parameters
    ----------
    counts : numpy array
        count matrix
    thr : int, optional
        relative filter threshold. The default is 5.

    Returns
    -------
    counts : numpy array
        filtered count matrix

    """

    counts[counts < thr] = 0
    
    return counts

def get_pool_table_line(id_, bowtie_res, all_ids, pools, poolidx_dict):
    """
    

    Parameters
    ----------
    id_ : int
        hashed ids encoding genomic coordinate, barcode and 
        transposon orientation
    bowtie_res : pandas dataframe
        parsed output of bowtie2, each line corresponds to a read
    all_ids : numpy array
        all unique hashed ids encoding genomic coordinate, barcode and 
        transposon orientation
    pools : list of strings
       all pools
    poolidx_dict : dict
        of poolnames to indices

    Returns
    -------
    line : numpy array
       line of count matrix for mutant id
    line_id : int
        hashed id encoding genomic coordinate, barcode and 
        transposon orientation

    """
    coord_subset = bowtie_res[table_lookup(id_, all_ids)]
    
    # pools are in column 1
    poolcounts = np.unique(coord_subset[:,1], return_counts=True)
    line_id = np.unique(coord_subset[:,6])[0]
    line = np.zeros((1,len(pools)), dtype=int)
    for j in range(len(poolcounts[0])):
        line[0,poolidx_dict[poolcounts[0][j]]] = poolcounts[1][j]

    return line, line_id


def get_count_matrix(ids, bowtie_res, experiment):
    """
    Assembles count matrix from bowtie2 output.

    Parameters
    ----------
    ids : int
        hashed ids encoding genomic coordinate, barcode and 
        transposon orientation
    bowtie_res : pandas dataframe
        parsed output of bowtie2, each line corresponds to a read
    experiment : LibraryExperiment
        Object holding input parameters.

    Returns
    -------
    out : pandas dataframe
        Mutant count matrix

    """
    # get all pools
    pools = experiment.pools
    
    # make sure pools are always on same position
    poolidx_dict = {p : i for i,p in enumerate(pools)}
    
    unique_ids, counts = np.unique(ids, return_counts=True)
    # exclude barcodes which occur only < experiment.min_counts
    unique_ids = unique_ids[counts>=experiment.min_counts]
    print("Generating count matrix with", len(unique_ids), "mutants")

    pool = mp.Pool(experiment.cores)
    pool_table_lines,line_ids = zip(*pool.starmap_async(get_pool_table_line, 
                                          zip(unique_ids, 
                                              repeat(bowtie_res.values), 
                                              repeat(ids),
                                              repeat(pools), 
                                              repeat(poolidx_dict))
                                          ).get())

    pool.close()
    pool_presence_table = np.vstack(pool_table_lines)
    table_ids  = np.hstack(line_ids)
    #orient = [i[0] for i in table_ids]
    # coord_match  = re.compile(r'(\d+)(\D+)')
    orient = []
    coord = []
    bar = []
    ref =  []
    for i in table_ids:
        o, c, b, r = i.split(";")
        orient.append(o)
        coord.append(c)
        bar.append(b)
        ref.append(r)
        # m = coord_match.search(i)
        # coord.append(m.group(1))
        # bar.append(m.group(2))
    
    coord = np.array(coord, dtype=int)
    orient = np.array(orient)
    bar = np.array(bar)
    ref = np.array(ref)
    out_ids = np.vstack((coord, orient, bar, ref)).T
    out = pd.DataFrame(np.hstack((out_ids,pool_presence_table)), columns = ["Feature","Orientation", "Barcode", "Reference"] + list(pools) )
    out["Feature"] = out["Feature"].astype(int)
    # important to sort for barcode filtering
    out = out.sort_values(["Feature", "Orientation"], axis=0)
    out = out.reset_index(drop=True)
    print("Generated count matrix!")
 
    return out


def get_hash_ids(ids):
    """
    Concatenate orientation, coordinate and barcode to create unique ids
    Return unique int identifiers along with original str ids

    Parameters
    ----------
    ids : list of str
        Ids indicating genomic coordinate, barcode and Tn orientation of mutant

    Returns
    -------
    encoded_id : list of int
        hashed id

    """

    # generate unique int id by creating a hash from the string and convert to int
    # numba only accepts int or bools
    encoded_id = []
    for i in ids:
        m = hashlib.md5()
        m.update(i.encode("utf8"))
        encoded_id.append(str(int(m.hexdigest(), 16))[0:16])
    encoded_id = np.array(encoded_id, dtype=int)
    
    return encoded_id

    
    


    
    