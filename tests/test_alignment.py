import pytest
from pysudoku_package.libraryexperiment import LibraryExperiment
import os
from pysudoku_package.alignment import  parse_bowtie2_output
from pysudoku_package.config import set_up_tmpdir


def test_parse_bowtie2_output():
    experiment = LibraryExperiment(8,30,10,"test_data/gb_ref/", "test_data/bowtie_ref/UTI89", 
                                    "AGATGTGTATAAGAGACAG", 1, "test_data", "test_data/exp_design.csv",
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr=0.05, 
                                    global_filter_thr = 5, min_counts=5)
    truecontents = open(os.path.join("test_data","alignment_result_expected.csv")).read()
    experiment.alignment = os.path.join("test_data", "test_alignment.sam")
    set_up_tmpdir()

    try:
        parse_bowtie2_output(experiment)
        contents = open(os.path.join("temp","alignment_result.csv")).read()
    finally:
        os.remove(os.path.join("temp","alignment_result.csv"))
    assert contents == truecontents
