import pytest
from pysudoku_package.libraryexperiment import LibraryExperiment
import os
from pysudoku_package.alignment import  parse_bowtie2_output
from pysudoku_package.config import set_up_tmpdir


def test_parse_bowtie2_output(cdir):
    experiment = LibraryExperiment(8,30,10,os.path.join(str(cdir),"gb_ref/"), 
                                   os.path.join(str(cdir),"bowtie_ref/UTI89"), 
                                    "AGATGTGTATAAGAGACAG", 1, str(cdir), 
                                    os.path.join(str(cdir),"exp_design.csv"),
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr=0.05, 
                                    global_filter_thr = 5, min_counts=5)
    truecontents = open(os.path.join(cdir,"alignment_result_expected.csv")).read()
    experiment.alignment = os.path.join(cdir, "test_alignment.sam")
    set_up_tmpdir()

    try:
        parse_bowtie2_output(experiment)
        contents = open(os.path.join(str(cdir),"temp","alignment_result.csv")).read()
    finally:
        os.remove(os.path.join(str(cdir),"temp","alignment_result.csv"))
    assert contents == truecontents
