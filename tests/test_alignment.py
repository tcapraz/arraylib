import pytest
from arraylib.libraryexperiment import LibraryExperiment
import os
from arraylib.alignment import  parse_bowtie2_output
from arraylib.config import set_up_tmpdir


def test_parse_bowtie2_output(ddir,cdir):
    experiment = LibraryExperiment(8,30,10,os.path.join(str(ddir),"gb_ref/"), 
                                   os.path.join(str(ddir),"bowtie_ref/UTI89"), 
                                    "AGATGTGTATAAGAGACAG", 1, str(ddir), 
                                    os.path.join(str(ddir),"exp_design.csv"),
                                    True, "CGAGGTCTCT", "CGTACGCTGC", filter_thr=0.05, 
                                    global_filter_thr = 5, min_counts=5)
    truecontents = open(os.path.join(ddir,"alignment_result_expected.csv")).read()
    experiment.alignment = os.path.join(ddir, "test_alignment.sam")
    set_up_tmpdir()


    parse_bowtie2_output(experiment)

    contents = open(os.path.join("temp","alignment_result.csv")).read()
 
    os.remove(os.path.join("temp","alignment_result.csv"))
    assert contents == truecontents
