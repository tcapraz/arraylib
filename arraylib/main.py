from pysudoku_package.libraryexperiment import LibraryExperiment
import sys
import click
import pandas as pd


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
def deconvolve(count_matrix, exp_design, cores, gb_ref, filter_thr, global_filter_thr) :
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
def deconvolve_validation(count_matrix, exp_design, cores, gb_ref, filter_thr, global_filter_thr) :
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