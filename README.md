# arraylib-solve

[![PyPI version](https://badge.fury.io/py/arraylib-solve.svg)](https://badge.fury.io/py/arraylib-solve)

# Introduction

`arraylib-solve` is a tool to deconvolve combinatorially pooled arrayed random mutagenesis libraries (e.g. by transposon mutagenesis). In a typical experiment generating arrayed mutagenesis libraries, first a pooled version of the library is created and arrayed on a grid of well plates. To infer the identities of each mutant on the well plate, wells are pooled in combinatorial manner such that each mutant appears in a unique combination of pools. The pools are then sequenced using NGS and sequenced reads are stored in individual fastq files per pool. `arraylib-solve` deconvolves the pools and returns summaries stating the identity and location of each mutant on the original well grid. The package is based on the approach described in [[1]](#1).

# Installation as python package

To install `arraylib` first create `Python 3.8` environment e.g. by

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

# Installation with singularity

Install singularity on your system 

```
conda create -n singularity -c conda-forge singularity -y
```

Active the conda environment
```
conda activate singularity
```

Download the docker image from dockerhub. This will write a singularity container named `arraylib_latest.sif` into your current work directory.
```
singularity pull docker://tcapraz/arraylib:latest

```

Start an interactive session of the container. Importantly, you need to `--bind` the directory with the input files. 
Note: the `:rw` at the end of the path is crucial for singularity to obtain read/write permission and hence be able to compute.
```
singularity shell --bind /path/to/folder/containing/input:/data:rw \
                  arraylib_latest.sif
```

# How to run `arraylib`

A detailed manual how to run `arraylib` interactively and from the command line can be found here https://tcapraz.github.io/arraylib/index.html.

# References
<a id="1">[1]</a> 
Baym, M., Shaket, L., Anzai, I.A., Adesina, O. and Barstow, B., 2016. Rapid construction of a whole-genome transposon insertion collection for Shewanella oneidensis by Knockout Sudoku. Nature communications, 7(1), p.13270.\
<a id="2">[2]</a> 
Langmead, B. and Salzberg, S.L., 2012. Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), pp.357-359.

