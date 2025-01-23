import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def simulate_required_arraysize(data, 
                                minsize,
                                maxsize,
                                number_of_simulations, 
                                number_of_repeats=30, 
                                gene_start = 0.1, 
                                gene_end = 0.9):
    """
    Run simulation of how many mutants need to be picked for the arrayed 
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
    arrsize=np.linspace(minsize,maxsize, number_of_simulations)
    data = pd.read_csv(data)
    
    sim_result = simulate_unique_genes(data, 
                                       arrsize,
                                       number_of_repeats, 
                                       gene_start, 
                                       gene_end)
    plot = plot_arraysize_vs_unique_genes(sim_result)
    
    return sim_result, plot

def simulate_unique_genes(data,
                        arrsize,
                        number_of_repeats=30, 
                        gene_start = 0.1, 
                        gene_end = 0.9):
    """Simulate how many unique genes are picked for a given array size 
    and mutant distribution of a pooled library. Genes can be filtered whether 
    only transposon hits are considered between gene_start and gene_end.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame of tnseeker output
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

    Returns
    -------
    pd.DataFrame
        DataFrame with mean and standard deviation of the number of unique genes and array size

    """
        
    dist = data["Read Counts"]/np.sum(data["Read Counts"])

    result_df = pd.DataFrame(columns = ["Arraysize", "Mean", "Std"], index = arrsize)

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
    matplotlib.Axis
    """
    result_df = result_df.astype(float)
    fig, ax = plt.subplots()
    ax.plot(result_df["Arraysize"], result_df["Mean"], c="black")
    ax.fill_between(result_df["Arraysize"], result_df["Mean"]-result_df["Std"], result_df["Mean"]+result_df["Std"], alpha=0.5)
    ax.set_xlabel("Number of mutants in arrayed library")
    ax.set_ylabel("Number of unique genes in arrayed library")
    ax.spines[['right', 'top']].set_visible(False)
    plt.tight_layout()

    return ax
