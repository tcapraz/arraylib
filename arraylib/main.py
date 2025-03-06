from arraylib.libraryexperiment import LibraryExperiment
import os
import click
import pandas as pd
import numpy as np
import matplotlib
from arraylib.simulations import simulate_required_arraysize, simulate_deconvolution, plot_precision_recall
from arraylib.config import run_tnseeker

@click.command()
@click.argument("input_dir", type=str)
@click.argument("exp_design", type=str)
@click.option("--cores", '-c',
              metavar='<str>',
              type=int,
              default=1,
              help='Number of cpu cores to use')
@click.option("--map_quality", '-mq',
              metavar='<int>',
              type=int,
              default=30,
              help='Minimum bowtie2 alignment quality score to include read')
@click.option("--seq_quality", '-sq',
              metavar='<int>',
              type=int,
              default=10,
              help='Minimum phred score for each base to include read')
@click.option("--gb_ref", '-gb',
              metavar='<str>',
              type=str,
              help='Path to genbank reference file')
@click.option("--bowtie_ref", '-br',
            metavar='<str>',
            type=str,
            help='Path to bowtie reference files')
@click.option("--tn_seq", '-t',
              metavar='<str>',
              type=str,
              help='Transposon sequence')
@click.option("--tn_mismatch", '-tm',
              metavar='<int>',
              type=int,
              default=1,
              help='Number of transposon mismatches allowed')
@click.option("--use_barcodes", '-b',
              type=bool,
              metavar='<bool>',
              default=True,
              help='Whether to use barcodes for deconvolution')
@click.option("--bar_upstream", '-bu',
              metavar='<str>',
              type=str,
              help='Upstream sequence of barcode')
@click.option("--bar_downstream", '-bd',
              metavar='<str>',
              type=str,
              help='Downstream sequence of barcode')
@click.option("--filter_thr", '-thr',
              metavar='<int>',
              type=float,
              default=0.05,
              help='Count filter threshold')
@click.option("--global_filter_thr", '-g_thr',
              metavar='<int>',
              type=float,
              default=5,
              help='Count filter threshold')
@click.option("--min_counts", '-mc',
              metavar='<int>',
              type=int,
              default=5,
              help='Minimum read counts for barcode to be inclucded')
def run(input_dir, exp_design, cores,map_quality,seq_quality,gb_ref, bowtie_ref, 
                                tn_seq, tn_mismatch, 
                                use_barcodes, bar_upstream, bar_downstream, filter_thr,
                                global_filter_thr, min_counts):
    """ 
    Run pysudoku to process sequencing data from combnatorially pooled 
    arrayed libraries and infer the most likely locations for each mutant.
    INPUT_DIR should be a path to the directory containing the input fastq
    files and EXP_DESIGN a csv file with columns Filename, Poolname and 
    Pooldimension indicating which files belong to which pool and pooling 
    dimension.
    """
    if tn_seq == None:
        raise ValueError("Please provide the transposon sequence!")
    if use_barcodes == True:
        if bar_upstream == None:
            raise ValueError("Please provide the sequence upstream to the transposon! If you don't want to consider Barcodes please set -b to False!")
        if bar_downstream == None:
            raise ValueError("Please provide the sequence downstream to the transposon! If you don't want to consider Barcodes please set -b to False!")
            
    
    experiment = LibraryExperiment(cores,map_quality,seq_quality,gb_ref, bowtie_ref, 
                                    tn_seq, tn_mismatch, input_dir, exp_design,
                                    use_barcodes, bar_upstream, bar_downstream, filter_thr,
                                    global_filter_thr,min_counts)
    experiment.pools
    print("Found", str(len(experiment.file2pool)), "input files!")
    print("Running pysudoku with", str(len(experiment.pools)), "pools and ",
          str(len(experiment.pool_dims)), "pooling dimensions!")
    print("Finding transposon sequence and trimming reads to genomic sequence!")
    experiment.get_genomic_seq()    
    print("Trimming done!")
    print("Starting alignment to reference!")    
    experiment.align_genomic_seq()       
    print("Alignment done!")    
    experiment.write_count_matrix()
    print("Inferring most likely locations and writing location summary!")
    experiment.deconvolve()
    print("Done!")

@click.command()
@click.argument("input_dir", type=str)
@click.argument("exp_design", type=str)
@click.option("--cores", '-c',
              metavar='<str>',
              type=int,
              default=1,
              help='Number of cpu cores to use')
@click.option("--map_quality", '-mq',
              metavar='<int>',
              type=int,
              default=30,
              help='Minimum bowtie2 alignment quality score to include read')
@click.option("--seq_quality", '-sq',
              metavar='<int>',
              type=int,
              default=10,
              help='Minimum phred score for each base to include read')
@click.option("--tn_mismatch", '-tm',
              metavar='<int>',
              type=int,
              default=1,
              help='Number of transposon mismatches allowed')
@click.option("--bar_upstream", '-bu',
              metavar='<str>',
              type=str,
              help='Upstream sequence of barcode')
@click.option("--bar_downstream", '-bd',
              metavar='<str>',
              type=str,
              help='Downstream sequence of barcode')
@click.option("--filter_thr", '-thr',
              metavar='<int>',
              type=float,
              default=0.05,
              help='Count filter threshold')
@click.option("--global_filter_thr", '-g_thr',
              metavar='<int>',
              type=float,
              default=5,
              help='Count filter threshold')
@click.option("--min_counts", '-mc',
              metavar='<int>',
              type=int,
              default=20,
              help='Minimum read counts for barcode to be inclucded')
def run_on_barcodes(input_dir, exp_design, cores,map_quality,seq_quality, 
                                 tn_mismatch, 
                                 bar_upstream, bar_downstream, filter_thr,
                                 global_filter_thr, min_counts):
    """ 
    Run pysudoku in barcode only mode, where no alignment to the reference 
    genome is performed.
    INPUT_DIR should be a path to the directory containing the input fastq
    files and EXP_DESIGN a csv file with columns Filename, Poolname and 
    Pooldimension indicating which files belong to which pool and pooling 
    dimension.
    """
   
    use_barcodes = True
    bowtie_ref = None
    gb_ref = None
    tn_seq = None
    experiment = LibraryExperiment(cores,map_quality,seq_quality,gb_ref, bowtie_ref, 
                                    tn_seq, tn_mismatch, input_dir, exp_design,
                                    use_barcodes, bar_upstream, bar_downstream, filter_thr, 
                                    global_filter_thr, min_counts)
    experiment.pools
    print("Found", str(len(experiment.file2pool)), "input files!")
    print("Running pysudoku with", str(len(experiment.pools)), "pools and ",
          str(len(experiment.pool_dims)), "pooling dimensions!")
    print("Finding transposon sequence and trimming reads to genomic sequence!")
    experiment.get_genomic_seq(barcode_only=True)    
    print("Barcode detection done!")     
    experiment.write_count_matrix(barcode_only=True)
    print("Inferring most likely locations and writing location summary!")
    experiment.deconvolve(barcode_only=True)
    print("Done!")

@click.command()
@click.argument("count_matrix", type=str)
@click.argument("exp_design", type=str)
@click.option("--cores", '-c',
              metavar='<str>',
              type=int,
              default=1,
              help='Number of cpu cores to use')
@click.option("--gb_ref", '-gb',
              metavar='<str>',
              type=str,
              help='Path to genbank reference file')
@click.option("--filter_thr", '-thr',
              metavar='<int>',
              type=float,
              default=0.05,
              help='Count filter threshold')
@click.option("--global_filter_thr", '-g_thr',
              metavar='<int>',
              type=float,
              default=5,
              help='Count filter threshold')
def deconvolve(count_matrix, exp_design, cores, gb_ref, filter_thr, global_filter_thr):
    """ 
    Run pysudoku to infer the most likely locations for each mutant from a 
    precomputed count matrix.
    COUNT_MATRIX should be a path to a precomputed count matrix and EXP_DESIGN 
    a csv file with columns Filename, Poolname and 
    Pooldimension indicating which files belong to which pool and pooling 
    dimension.
    """
    if count_matrix == None:
        raise ValueError("Please provide path to count matrix!")
   
    experiment = LibraryExperiment(cores,map_quality =None,seq_quality=None,gb_ref=gb_ref, bowtie_ref=None, 
                                    tn_seq=None, tn_mismatches=None, input_dir=None, exp_design=exp_design,
                                    use_barcodes=True, bar_upstream=None, bar_downstream=None, filter_thr=filter_thr, 
                                    global_filter_thr=global_filter_thr, min_counts=5)
    c = pd.read_csv(count_matrix)
    experiment.count_mat = c
    print("Inferring most likely locations and writing location summary!")
    experiment.deconvolve()
    print("Done!")

@click.command()
@click.argument("count_matrix", type=str)
@click.argument("exp_design", type=str)
@click.option("--cores", '-c',
              metavar='<str>',
              type=int,
              default=1,
              help='Number of cpu cores to use')
@click.option("--gb_ref", '-gb',
              metavar='<str>',
              type=str,
              help='Path to genbank reference file')
@click.option("--filter_thr", '-thr',
              metavar='<int>',
              type=float,
              default=0.05,
              help='Count filter threshold')
@click.option("--global_filter_thr", '-g_thr',
              metavar='<int>',
              type=float,
              default=5,
              help='Count filter threshold')
def deconvolve_validation(count_matrix, exp_design, cores, gb_ref, filter_thr, global_filter_thr):
    """ 
    Run pysudoku to infer the most likely locations for each mutant from a 
    precomputed count matrix.
    COUNT_MATRIX should be a path to a precomputed count matrix and EXP_DESIGN 
    a csv file with columns Filename, Poolname and 
    Pooldimension indicating which files belong to which pool and pooling 
    dimension.
    """
    if count_matrix == None:
        raise ValueError("Please provide path to count matrix!")
   
    experiment = LibraryExperiment(cores,map_quality =None,seq_quality=None,gb_ref=gb_ref, bowtie_ref=None, 
                                    tn_seq=None, tn_mismatches=None, input_dir=None, exp_design=exp_design,
                                    use_barcodes=True, bar_upstream=None, bar_downstream=None, filter_thr=filter_thr, 
                                    global_filter_thr=global_filter_thr, min_counts=5)
    c = pd.read_csv(count_matrix)
    experiment.count_mat = c
    print("Inferring most likely locations and writing location summary!")
    experiment.deconvolve_validation()
    print("Done!")

@click.command()
@click.option("--tnseeker_output_path", 
                type=str,
                default=None,
                help='Path to tnseeker mapped insertions file. If none is provided, tnseeker will be automatically ran and create one.')
@click.option("--minsize", '-min',
              metavar='<int>',
              type=int,
              default=1000,
              help='Minimum array size to simulate')
@click.option("--maxsize", '-max',
              metavar='<int>',
              type=int,
              default=5000,
              help='Maximum array size to simulate')
@click.option("--number_of_simulations", '-nsim',
              metavar='<int>',
              type=int,
              default=100,
              help='Number of simulations of sizes interpolated between minsize and maxsize')
@click.option("--number_of_repeats", '-niter',
              metavar='<int>',
              type=int,
              default=10,
              help='Number of repeats of simulations to perform')
@click.option("--gene_start", '-gstart',
              metavar='<float>',
              type=float,
              default=0.1,
              help='Gene start')
@click.option("--gene_end", '-gend',
              metavar='<float>',
              type=float,
              default=0.9,
              help='Gene end')
@click.option("--output_filename", '-o',
              metavar='<str>',
              type=str,
              default="arraysize_simulation.csv",
              help='Path for simulation output file')
@click.option("--output_plot", '-p',
              metavar='<str>',
              type=str,
              default="arraysize_simulation.pdf",
              help='Path for simulation output plot')
@click.option("--seed", '-s',
              metavar='<int>',
              type=int,
              default=None,
              help='Seed for random number generator')
@click.option("--input_dir", 
                type=str)
@click.option("--cores", '-c',
              metavar='<str>',
              type=int,
              default=1,
              help='Number of cpu cores to use')
@click.option("--map_quality", '-mq',
              metavar='<int>',
              type=int,
              default=30,
              help='Minimum bowtie2 alignment quality score to include read')
@click.option("--seq_quality", '-sq',
              metavar='<int>',
              type=int,
              default=10,
              help='Minimum phred score for each base to include read')
@click.option("--gb_ref", '-gb',
              metavar='<str>',
              type=str,
              help='Path to the reference file + FASTA (both are required)')
@click.option("--tn_seq", '-t',
              metavar='<str>',
              type=str,
              help='Transposon sequence')
@click.option("--tn_mismatch", '-tm',
              metavar='<int>',
              type=int,
              default=1,
              help='Number of transposon mismatches allowed')
@click.option("--annotation_type", '-at',
              metavar='<str>',
              type=str,
              default="GB",
              help='Reference file type: GB or GFF')
@click.option("--seq_type", '-st',
              metavar='<str>',
              type=str,
              default="SE",
              help='Single-end (SE) or pair-end (PE)')
@click.option("--name", '-na',
              metavar='<str>',
              type=str,
              default="tnseeker_out",
              help='Name of the reference files (both reference and FASTA need to match)')
def plot_required_arraysize(tnseeker_output_path, minsize, maxsize, number_of_simulations, 
                                number_of_repeats, gene_start, 
                                gene_end, output_filename, output_plot, seed, input_dir, cores,
                                map_quality, seq_quality, gb_ref, tn_seq, tn_mismatch, 
                                annotation_type, seq_type, name):
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
    output_filename : str, optional
        Filepath for simulation output.
        by default : arraysize_simulation.csv
    output_plot : str, optional
        Filepath for simulation plot.
        by default : arraysize_simulation.pdf
        
    Returns
    -------
    pd.DataFrame
        DataFrame with mean and standard deviation of the number of unique genes and array size
    matplotlib.Axis
        Scatter plot of unique genes vs size of arrayed library 
    """

    if tnseeker_output_path is None:
        run_tnseeker(input_dir,
                    cores,
                    map_quality,
                    seq_quality,
                    gb_ref, 
                    tn_seq, 
                    tn_mismatch, 
                    annotation_type,
                    seq_type,
                    name)
        tnseeker_output_path = os.path.join(os.getcwd(), name, f"all_insertions_{name}.csv")

    matplotlib.use('Agg')
    print("Simulating unique genes for arraysizes between", str(minsize), 
          "and", str(maxsize), "!")

    res, plot = simulate_required_arraysize(data = tnseeker_output_path, 
                                minsize = minsize,
                                maxsize = maxsize,
                                number_of_simulations = number_of_simulations, 
                                number_of_repeats = number_of_repeats, 
                                gene_start = gene_start,
                                gene_end = gene_end,
                                seed = seed)
    
    print("Plotting results to", str(output_plot), "!")
    plot.savefig(output_plot)
    print("Saving results to", str(output_filename), "!")
    res.to_csv(output_filename, index=False)
    print("Done!")
    
@click.command()
@click.option("--input_dir",
              type=str)
@click.option("--tnseeker_output_path", 
                type=str,
                default=None,
                help='Path to tnseeker mapped insertions file. If none is provided, tnseeker will be automatically ran and create one.')
@click.option("--minsize", '-min',
              metavar='<int>',
              type=int,
              default=10,
              help='Minimum grid size to simulate')
@click.option("--maxsize", '-max',
              metavar='<int>',
              type=int,
              default=20,
              help='Maximum grid size to simulate')
@click.option("--number_of_simulations", '-nsim',
              metavar='<int>',
              type=int,
              default=2,
              help='Number of simulations of sizes interpolated between minsize and maxsize')
@click.option("--output_filename", '-o',
              metavar='<str>',
              type=str,
              default="precision_recall_simulation.csv",
              help='Path for simulation output file')
@click.option("--output_plot", '-p',
              metavar='<str>',
              type=str,
              default="simulation.pdf",
              help='Base path for simulation output plot')
@click.option("--seed", '-s',
              metavar='<int>',
              type=int,
              default=None,
              help='Seed for random number generator')
@click.option("--cores", '-c',
              metavar='<str>',
              type=int,
              default=1,
              help='Number of cpu cores to use')
@click.option("--map_quality", '-mq',
              metavar='<int>',
              type=int,
              default=30,
              help='Minimum bowtie2 alignment quality score to include read')
@click.option("--seq_quality", '-sq',
              metavar='<int>',
              type=int,
              default=10,
              help='Minimum phred score for each base to include read')
@click.option("--gb_ref", '-gb',
              metavar='<str>',
              type=str,
              help='Path to the reference file + FASTA (both are required)')
@click.option("--tn_seq", '-t',
              metavar='<str>',
              type=str,
              help='Transposon sequence')
@click.option("--tn_mismatch", '-tm',
              metavar='<int>',
              type=int,
              default=1,
              help='Number of transposon mismatches allowed')
@click.option("--annotation_type", '-at',
              metavar='<str>',
              type=str,
              default="GB",
              help='Reference file type: GB or GFF')
@click.option("--seq_type", '-st',
              metavar='<str>',
              type=str,
              default="SE",
              help='Single-end (SE) or pair-end (PE)')
@click.option("--name", '-na',
              metavar='<str>',
              type=str,
              default="SE",
              help='Name of the reference files (both reference and FASTA need to match)')
def plot_expected_deconvolution_accuracy(tnseeker_output_path, 
                                         minsize, 
                                         maxsize, 
                                         number_of_simulations, 
                                         output_filename, 
                                         output_plot,
                                         seed, input_dir, cores,
                                         map_quality, seq_quality, gb_ref, tn_seq, tn_mismatch, 
                                         annotation_type, seq_type, name):
    """
    Run simulations of count matrices followed by deconvolution. For each simulation
    precision and recall are calculated for a 3D and 4D grid of well plates. 
    Simulations are based 
    on the mutant distribution of a pooled library. 
    
    Parameters
    ----------
    data : str
        filepath to tnseeker output file all_insertions.csv. 
        It should contain the columns "Read Counts", "Gene Name" and "Relative Position in Gene (0-1)".
    minsize : int
        Minimum grid size to simulate.
    maxsize : int
       Maximum grid size to simulate.
    number_of_simulations : int
        Number of simulations to perform between minsize and maxsize. 
        I.e. if the number_of_simulations is 2, 
        only simulations of minsize and maxsize are performed. 
    output_filename : str, optional
        Filepath for simulation output.
        by default : accuracy_simulation.csv
    output_plot : str, optional
        Filepath for simulation plot.
        by default : accuracy_simulation.pdf
        
    Returns
    -------
    pd.DataFrame
        DataFrame with precision and recall for a simulated deconvolution for 
        a 3D and 4D grid of well plates.
    matplotlib.Axis
        Scatter plot of number of mutants vs recall (or precision)
    """
    if tnseeker_output_path is None:
        run_tnseeker(input_dir,
                    cores,
                    map_quality,
                    seq_quality,
                    gb_ref, 
                    tn_seq, 
                    tn_mismatch, 
                    annotation_type,
                    seq_type,
                    name)
        tnseeker_output_path = os.path.join(os.getcwd(), name, f"all_insertions_{name}.csv")
    
    matplotlib.use('Agg')
    gridsize=np.linspace(minsize,maxsize, number_of_simulations, dtype=int)

    print("Simulating deconvolution runs for grids of sizes", gridsize.astype(str), "!")

    res = simulate_deconvolution(data = tnseeker_output_path, 
                                minsize = minsize,
                                maxsize = maxsize,
                                number_of_simulations = number_of_simulations,
                                seed = seed) 
    
    prec_plot, rec_plot = plot_precision_recall(res)

    print("Plotting results to precision_"+str(output_plot), "!")
    prec_plot.savefig("precision_" +output_plot)
    print("Plotting results to recall_"+str(output_plot), "!")
    rec_plot.savefig("recall_" + output_plot)
    print("Saving results to", str(output_filename), "!")
    res.to_csv(output_filename, index=False)
    print("Done!")