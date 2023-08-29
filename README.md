# arraylib-solve

# Introduction

`arraylib-solve` is a tool to deconvolve combinatorially pooled arrayed random mutagenesis libraries (e.g. by transposon mutagenesis). In a typical experiment generating arrayed mutagenesis libraries, first a pooled version of the library is created and arrayed on a grid of well plates. To infer the identities of each mutant on the well plate, wells are pooled in combinatorial manner such that each mutant appears in a unique combination of pools. The pools are then sequenced using NGS and sequenced reads are stored in individual fastq files per pool. `arraylib-solve` deconvolves the pools and returns summaries stating the identity and location of each mutant on the original well grid. The package is based on the approach described in [[1]](#1).

# Installation

To install `arraylib-solve` first create `Python 3.8` environment e.g. by

```
conda create --name arraylib-env python=3.8
conda activate arraylib-env
```

and install the package using 

```
pip install arraylib-solve
```

`arraylib-solve` uses bowtie2 [[2]](#2) to align reads to the reference genome. Please ensure that bowtie2 is installed in your environment by running:

```
conda install -c bioconda bowtie2
```


# How to run `arraylib-solve`

To run `arraylib-solve` on a library deconvolution experiment with default parameters run:

```
arraylib-run <input_directory> <experimental_design.csv> -c <number_of_cpu_cores_to_use> -gb <path_to_genbank_reference_directory> -br <path_to_bowtie2_indices> -t <transposon_sequence> -bu <upstream_sequence_of_barcodes> -bd <downstream_sequence_of_barcodes>
```

## Input parameters

Required parameters:
* input_dir: path to directory holding the input fastq files
* exp_design: path to csv file indicating experimental design (values should be separated by a comma). The experimental design file 
       should have columns, Filename, Poolname and Pooldimension. (see example in tests/test_data/full_exp_design.csv)
  * Filename should contain all the unqiue input fastq filenames.
  * Poolname should indicate to which pool a given file belongs. Multiple files per poolname are allowed.
  * Pooldimension indicates the pooling dimension a pool belongs to. All pools sharing the same pooling dimension should have the same string in the Pooldimension column.
  

An example of how an exp_design file could look like:

| Filename          | Poolname        | Pooldimension  |
| :---------------: | :-------------: | :------------: |
| column1.fastq     | column1         | columns        |
| column2.fastq     | column2         | columns        |
| row1.fastq        | row1            | rows           |
| row2.fastq        | row2            | rows           |
| platerow1.fastq   | platerow1       | platerows      |
| platerow2.fastq   | platerow2       | platerows      |
| platecol1.fastq   | platecol1       | platecols      |
| platecol2.fastq   | platecol2       | platecols      |

* -gb path to genbank reference file
* -br path to bowtie index files, ending with the basename of your index (if the basename of your index is UTI89 and you store your bowtie2 references in bowtie_ref it should be bowtie_ref/UTI89). Please visit https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer for a manual how to create bowtie2 indices.
* -t transposon sequence (e.g. AGATGTGTATAAGAGACAG)
* -bu upstream sequence of barcode (e.g. CGAGGTCTCT)
* -bd downstream sequence of barcode (e.g. CGTACGCTGC)

Optional parameters:
* -mq minimum bowtie2 alignment quality score for each base to include read
* -sq minimum phred score for each base to include read
* -tm number of transposon mismatches allowed
* -thr threshold for local filter (e.g. a threshold of 0.05 would filter out all reads < 0.05 of the maximum read count for a given mutant)
* -g\_thr threshold for global filter (all reads below g_thr will be set to 0) 

## Run only on barcodes
If you want to run arraylib-solve only on barcodes without alignment to the reference genome use the following command:

```
arraylib-run_on_barcodes <input_directory> <experimental_design.csv> -c <number_of_cpu_cores_to_use>  -bu <upstream_sequence_of_barcodes> -bd <downstream_sequence_of_barcodes>
```

Optional parameters:

* -thr threshold for local filter (e.g. a threshold of 0.05 would filter out all reads < 0.05 of the maximum read count for a given mutant)
* -g\_thr threshold for global filter (all reads below g_thr will be set to 0) 

## Output

`arraylib-solve` outputs 4 files: 
* count_matrix.csv: Read counts per pool for each mutant.
* filtered_matrix.csv: Read counts per pool for each mutant, but mutants with barcodes with low read counts for a given genomic location are filtered out.
* mutant_location_summary.csv: A summary of mutants found in the well plate grid, where each row corresponds to a different mutant.
* well_location_summary.csv: A summary of the deconvolved well plate grid, where each row corresponds to a different well.



# References
<a id="1">[1]</a> 
Baym, M., Shaket, L., Anzai, I.A., Adesina, O. and Barstow, B., 2016. Rapid construction of a whole-genome transposon insertion collection for Shewanella oneidensis by Knockout Sudoku. Nature communications, 7(1), p.13270.\
<a id="2">[2]</a> 
Langmead, B. and Salzberg, S.L., 2012. Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), pp.357-359.

