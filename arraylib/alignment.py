#!/usr/bin/env python3

import subprocess
import os
from tnseeker.extras.helper_functions import bowtie2parser

from arraylib.io import txt_writer


def parse_bowtie2_output(experiment):
    """
    Parse output of bowtie2 and write to alignment_result.csv in 
    temp directory.

    Parameters
    ----------
    experiment : LibraryExperiment
        Object holding input parameters.

    Returns
    -------
    None.

    """
    flag_list = [0, 16]
    # if param['seq_type']=="PE": #paired ended
    #     flag_list = [83, 99] 

    out=[]
    outfile = os.path.join("temp", "alignment_result.csv")
    align_res = experiment.alignment
    with open(align_res) as current:
        for line in current:
            line = line.split('\t')
            sam_output = bowtie2parser(line,experiment.map_quality,flag_list)

            if "aligned_valid_reads" in sam_output:
                read_id,pool,barcode = line[0].split(":")
                out.append(f"{read_id},{pool},{sam_output['local']},{sam_output['orientation']},{barcode},{sam_output['contig']}\n")

            if len(out)>500000:
                txt_writer(outfile, out)
                out=[]
                        
    txt_writer(outfile, out)


def run_bowtie2(experiment):
    """
    Calls bowtie2 alignment from subprocess and runs alignment on 
    trimmed_sequences.fastq. Bowtie2 writes alignments to temp/aligment.sam.

    Parameters
    ----------
    experiment : LibraryExperiment
        Object holding input parameters.


    -------
    None.

    """

    out = experiment.alignment
    seq = os.path.join("temp", "trimmed_sequences.fastq")
    subprocess.call(['bowtie2', 
                     '-x', 
                     experiment.bowtie_ref, 
                     '-q', 
                     seq, 
                     "--end-to-end", 
                     f"-S {out}",
                     '--no-unal',
                     '--threads', 
                     str(experiment.cores)])

