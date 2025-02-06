import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import string
from arraylib.predict_locations import getX, fit_nnomp
from itertools import product
import seaborn as sns

def simulate_required_arraysize(data, 
                                minsize,
                                maxsize,
                                number_of_simulations, 
                                number_of_repeats=30, 
                                gene_start = 0.1, 
                                gene_end = 0.9,
                                seed=None):
    """
    Wrapper for simulation of how many mutants need to be picked for the arrayed 
    library to reach a certain number of unique genes. Simulations are based 
    on the mutant distribution of a pooled library. 
    Genes can be filtered whether only transposon hits are considered between 
    gene_start and gene_end.
    
    Parameters
    ----------
    data : str
        filepath to tnseeker output file all_insertions.csv. 
        It should contain the columns "Read Counts", "Gene Name" and "Relative Position in Gene (0-1)".
    minsize : int
        Minimum array size to simulate.
    maxsize : int
       Maximum array size to simulate.
    number_of_simulations : int
        Number of simulations to perform between minsize and maxsize. 
        I.e. if the number_of_simulations is 2, 
        only simulations of minsize and maxsize are performed. 
    number_of_repeats : int, optional
        Number of times the simulations are repeated, by default 30
    gene_start : float, optional
        Minimum distance to the start of a gene to be counted as a transposon hit, 
        by default 0.1
    gene_end : float, optional
        Maximum distance to the start of a gene to be counted as a transposon hit, 
        by default 0.9
    
    Returns
    -------
    pd.DataFrame
        DataFrame with mean and standard deviation of the number of unique genes and array size
    matplotlib.Axis
        Scatter plot of unique genes vs size of arrayed library 

    """
    arrsize=np.linspace(minsize,maxsize, number_of_simulations, dtype=int)
    
    
    sim_result = simulate_unique_genes(data, 
                                       arrsize,
                                       number_of_repeats, 
                                       gene_start, 
                                       gene_end,
                                       seed=seed)
    plot = plot_arraysize_vs_unique_genes(sim_result)
    
    return sim_result, plot

def simulate_unique_genes(data,
                        arrsize,
                        number_of_repeats=30, 
                        gene_start = 0.1, 
                        gene_end = 0.9, 
                        seed=None):
    """Simulate how many unique genes are picked for a given array size 
    and mutant distribution of a pooled library. Genes can be filtered whether 
    only transposon hits are considered between gene_start and gene_end.

    Parameters
    ----------
    data : str
        filepath to tnseeker output file all_insertions.csv. 
        It should contain the columns "Read Counts", "Gene Name" and "Relative Position in Gene (0-1)".
    arrsize : np.array
        array containing arraysizes to simulate
    number_of_repeats : int, optional
        Number of times the simulations are repeated, by default 30
    gene_start : float, optional
        Minimum distance to the start of a gene to be counted as a transposon hit, 
        by default 0.1
    gene_end : float, optional
        Maximum distance to the start of a gene to be counted as a transposon hit, 
        by default 0.9
    seed : int, optional
        random seed for simulations

    Returns
    -------
    pd.DataFrame
        DataFrame with mean and standard deviation of the number of unique genes and array size

    """
    data = pd.read_csv(data)    
    dist = data["Read Counts"]/np.sum(data["Read Counts"])

    result_df = pd.DataFrame(columns = ["Arraysize", "Mean", "Std"], index = arrsize)
    if seed is not None:
        np.random.seed(seed)
    for i in arrsize:
        result=[]
        for j in range(number_of_repeats):
            picks = np.random.choice(np.arange(dist.shape[0]), replace=True, size=int(i),p=dist)
            genes = data.loc[picks, "Gene Name"]
            gene_pos =  data.loc[picks, "Relative Position in Gene (0-1)"]
            gene_pos_mask = np.logical_and(gene_pos >gene_start, gene_pos< gene_end)
            genes_filtered = genes.loc[gene_pos_mask]
            result.append(len(genes_filtered.unique()))
        result_df.loc[i, "Arraysize"] = i
        result_df.loc[i, "Mean"] = np.mean(result)
        result_df.loc[i, "Std"] = np.std(result)
    return result_df

def plot_arraysize_vs_unique_genes(result_df):
    """Plot the number of unique genes picked for a given array size.

    Parameters
    ----------
    result_df : pd.DataFrame
        DataFrame output of simulate_unique_genes_for_arraysize
    Returns
    matplotlib.Figure
    """
    result_df = result_df.astype(float)
    fig, ax = plt.subplots()
    ax.plot(result_df["Arraysize"], result_df["Mean"], c="black")
    ax.fill_between(result_df["Arraysize"], result_df["Mean"]-result_df["Std"], result_df["Mean"]+result_df["Std"], alpha=0.5)
    ax.set_xlabel("Number of mutants in arrayed library")
    ax.set_ylabel("Number of unique genes in arrayed library")
    ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()

    return fig

def recall(pred, true):
    correct = 0
    for i in pred:
        p = pred[i]
        t = true[i]
        correct += len(np.intersect1d(p,t))
    return correct/len([x for v in true.values() for x in v])

def precision(pred,true):
    total = 0
    correct = 0
    for i in pred:
        p = pred[i]
        t = true[i]
        correct += len(np.intersect1d(p,t))
        total += len(p)
    return correct/total

    
def simulate_deconvolution(data,
                           minsize,
                           maxsize,
                           number_of_simulations,
                           seed=None
                           ):
    """
    Simulate deconvolution runs and calculate precision and recall. Simulations are based 
    on the mutant distribution of a pooled library. 
    
    Parameters
    ----------
    data : str
        filepath to tnseeker output file all_insertions.csv. 
        It should contain the columns "Read Counts", "Gene Name" and "Relative Position in Gene (0-1)".
    minsize : int
        Minimum array size to simulate.
    maxsize : int
       Maximum array size to simulate.
    number_of_simulations : int
        Number of simulations to perform between minsize and maxsize. 
        I.e. if the number_of_simulations is 2, 
        only simulations of minsize and maxsize are performed. 
    seed : int, optional
        random seed for simulations
    Returns
    -------
    pd.DataFrame
        DataFrame with precision and recall for 3D and 4D deconvolution
    
    """
     
    gridsize=np.linspace(minsize,maxsize, number_of_simulations, dtype=int)
    data = pd.read_csv(data)
    poolsize = data.shape[0]
    pool = data["Read Counts"]/np.sum(data["Read Counts"])
    
    result_df = pd.DataFrame(columns = ["Arraysize", "Precision_4D", "Recall_4D", "Precision_3D", "Recall_3D"], index = [s*s*96 for s in gridsize])
    if seed is not None:
        np.random.seed(seed)
    for  s in gridsize: 
        row=8
        col=12
        arrsize=s*s*96
    
        n_plates_3d = s*s
    
        rows = list(string.ascii_uppercase)[0:row] 
        cols = [str(i) for i in range(1,col+1)]
        prs = ["PR0"+str(i) if i < 10 else "PR" +str(i) for i in range(1,s+1)]
        pcs = ["PC0"+str(i) if i < 10 else "PC" +str(i) for i in range(1,s+1)]
        plates_3d = ["plate"+str(i+1) for i in range(n_plates_3d)]
        
        dims= {"Row": rows, "Col": cols, "PR": prs, "PC":pcs}
        dims_3d= {"Row": rows, "Col": cols, "Plates": plates_3d}
    
        X = getX(dims, rows+cols+prs+pcs)
        X_3d = getX(dims_3d, rows+cols+plates_3d)
    
        well_pools = list(product(*dims.values()))
        well_pools_3d = list(product(*dims_3d.values()))
    
        wells = np.array(["_".join(w) for w in well_pools])
        wells_3d = np.array(["_".join(w) for w in well_pools_3d])
    
    
        picks =np.random.choice(np.arange(len(pool)), replace=True, size=arrsize,p=pool)
    
        count_mat = pd.DataFrame(0, index=np.unique(picks), columns = rows+cols+prs+pcs)
        count_mat_3d = pd.DataFrame(0, index=np.unique(picks), columns = rows+cols+plates_3d)
    
        picked_wells = {m:[] for m in np.unique(picks)}
        picked_wells_3d = {m:[] for m in np.unique(picks)}
    
    
    
        for i in range(arrsize):
            mutant = picks[i]
            w =wells[i]
            picked_wells[mutant].append(str(w))
            p = w.split("_")
            poisson_mean = np.random.gamma(100,1)
            count_mat.loc[mutant,p] = count_mat.loc[mutant,p] + np.random.poisson(poisson_mean,size=4)
            #3D
            w = wells_3d[i]
            picked_wells_3d[mutant].append(str(w))
            p = w.split("_")
            poisson_mean = np.random.gamma(100,1)
            count_mat_3d.loc[mutant,p] = count_mat_3d.loc[mutant,p] + np.random.poisson(poisson_mean,size=3)
            
            loc = {}
            loc_3d = {}
        for i in range(count_mat.shape[0]):
            res = fit_nnomp(count_mat.iloc[i], X, dims)
            loc[count_mat.index[i]]=res[0]
            res=fit_nnomp(count_mat_3d.iloc[i], X_3d, dims_3d)
            loc_3d[count_mat_3d.index[i]]=res[0]

        result_df.loc[arrsize, "Recall_4D"] = recall(loc,picked_wells)
        result_df.loc[arrsize, "Precision_4D"] = precision(loc,picked_wells)
        result_df.loc[arrsize, "Recall_3D"] = recall(loc_3d,picked_wells_3d)
        result_df.loc[arrsize, "Precision_3D"] = precision(loc_3d,picked_wells_3d)
        result_df.loc[arrsize, "Arraysize"] = arrsize
    return result_df
    
def plot_precision_recall(result_df):
    """Plot the number of unique genes picked for a given array size.

    Parameters
    ----------
    result_df : pd.DataFrame
        DataFrame output of simulate_deconvolution
    Returns
    matplotlib.Figure
    """
    
    result_df = result_df.astype(float)
    cblind_colors = sns.color_palette("colorblind", as_cmap=True)[0:2]

    prec_fig, prec_ax = plt.subplots()
    prec_ax.plot(result_df["Arraysize"], result_df["Precision_4D"], c= cblind_colors[0], label="4D")
    prec_ax.scatter(result_df["Arraysize"], result_df["Precision_4D"], c= cblind_colors[0], s=10)
    prec_ax.plot(result_df["Arraysize"], result_df["Precision_3D"], c= cblind_colors[1], label="3D")
    prec_ax.scatter(result_df["Arraysize"], result_df["Precision_3D"], c= cblind_colors[1], s=10)

    prec_ax.set_xlabel("Number of mutants in arrayed library")
    prec_ax.set_ylabel("Precision")
    prec_ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()


    
    rec_fig, rec_ax = plt.subplots()
    rec_ax.plot(result_df["Arraysize"], result_df["Recall_4D"], c= cblind_colors[0], label="4D")
    rec_ax.scatter(result_df["Arraysize"], result_df["Recall_4D"], c= cblind_colors[0], s=10)
    rec_ax.plot(result_df["Arraysize"], result_df["Recall_3D"], c= cblind_colors[1], label="3D")
    rec_ax.scatter(result_df["Arraysize"], result_df["Recall_3D"], c= cblind_colors[1], s=10)

    rec_ax.set_xlabel("Number of mutants in arrayed library")
    rec_ax.set_ylabel("Recall")
    rec_ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()

    return prec_fig, rec_fig
