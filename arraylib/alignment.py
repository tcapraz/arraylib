#!/usr/bin/env python3

import subprocess
import os

from pysudoku_package.io import txt_writer


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
            if (line[0][0] != "@") and (line[0][0] != ">"): #ignores headers and unaligned contigs
                id_str = line[0]
                flag_sum = int(line[1])
                ref = line[2]
                coord = int(line[3])
                quality = int(line[4])
                cigar = line[5]
                multi = "XS:i:" in line
                orientation = '+'
                if flag_sum == flag_list[1]:
                    orientation = '-'
                    # considers CIGAR scores. S means that the read is clipped, so the position will actually be different
                    # do we need that in end to end mode?
                    # matches = findall(r'(\d+)([A-Z]{1})', cigar)
                    # clipped = 0
                    # for match in matches:
                    #     if match[1] == "S":
                    #         clipped=int(match[0])
                    #         break #only consideres the first one at the start
                    # coord = coord + len(line[9]) - clipped
                    # subtract 2 if tn is in reverse direction
                    coord = coord + len(line[9]) - 2
                if (flag_sum in flag_list) & (multi==False) & (quality >= experiment.map_quality):
                    read_id,pool,barcode = id_str.split(":")
                    out.append(f"{read_id},{pool},{coord},{orientation},{barcode},{ref}\n")
            #if psutil.virtual_memory().percent>=80:
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

