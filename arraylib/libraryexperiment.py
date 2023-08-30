#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from arraylib.config import get_input_files, get_pool_dicts, set_up_tmpdir
from arraylib.io import txt_writer, file_compiler
from arraylib.alignment import run_bowtie2, parse_bowtie2_output
from arraylib.trimming import read_data, detect_barcodes
from arraylib.generate_count_matrix import  get_count_matrix, get_hash_ids, filter_barcodes, local_filter_counts, global_filter_counts
from arraylib.deconvolution import get_ambiguity, get_unambiguous_locations, coords2genes, transpose_location_summary, transpose_barcode_summary
from arraylib.predict_locations import predict_ambiguous_locations
import multiprocessing as mp
from itertools import repeat
import os
from pandas import concat, read_csv
import numpy as np

class LibraryExperiment(object): 
        """
        LibraryExperiment class used to perform individual analysis steps and
        store intermediate results.
        
        Parameters
        ----------
        cores : int
            number of cores to use
        map_quality : int
            minimum bowtie2 alignment quality score for each base to include read
        seq_quality : int
            minimum phred score for each base to include read
        gb_ref : str
            path to genbank reference file
        bowtie_ref : str
            path to bowtie reference file
        tn_seq : str
            transposon sequence (e.g. ATTGCCTA)
        tn_mismatches : int
            number of transposon mismatches allowed
        input_dir : str
            path to directory holding the input fastq files
        exp_design : str
            path to file indicating experimental design. The experimental design file 
            should have columns, Filename, Poolname and Pooldimension
        use_barcodes : bool
            whether to perform deconvolution only on barcodes without genomic alignment
        bar_upstream : str
            upstream sequence of barcode
        bar_downstream : str
            downstream sequence of barcode
        filter_thr : float
            threshold for local filter. read counts of a given mutant whose 
            percentage of the max read count is lower than filter_thr are set 
            to 0
        global_filter_thr : int
            threshold for global filter. read counts lower than 
            global_filter_thr are set to 0



        """
        def __init__(self, 
                     cores, 
                     map_quality, 
                     seq_quality, 
                     gb_ref, 
                     bowtie_ref, 
                     tn_seq,
                     tn_mismatches,
                     input_dir, 
                     exp_design,
                     use_barcodes, 
                     bar_upstream,
                     bar_downstream,
                     filter_thr,
                     global_filter_thr,
                     min_counts):
        
            # store input parameters
            self.cores = cores
            self.map_quality = map_quality
            self.seq_quality = seq_quality
            self.gb_ref = gb_ref
            self.bowtie_ref = bowtie_ref
            self.tn_seq = tn_seq
            self.tn_mismatches = tn_mismatches
            self.input_dir = input_dir
            self.exp_design = exp_design
            self.use_barcodes = use_barcodes
            self.filter_thr = filter_thr
            self.global_filter_thr = global_filter_thr
            self.min_counts = min_counts
            if use_barcodes:
                self.bar_upstream = bar_upstream
                self.bar_downstream = bar_downstream
            
            # hardcode paths for intermediate results
            self.bowtie_res = os.path.join("temp", "alignment_result.csv")
            self.alignment = os.path.join("temp", "alignment.sam")
            
            
            self.file2pool, self.pool_dims = get_pool_dicts(self.exp_design)
            self.pools = np.sort(np.unique([i for j in list(self.pool_dims.values()) for i in j])).tolist()
        



        def get_genomic_seq(self, barcode_only=False):
            
            """ 
            Trims all the reads based on the presence of the transposon recognizing
            sequence. Only keeps the downstream genomic sequences after the 
            transposon border site.
            
            Trimmed genomic sequences are written to temp/trimmed_sequences.fastq
            
            """
            
            self.input_file_paths, self.input_file_names = get_input_files(self)
    
            set_up_tmpdir()
            
            # if barcode only mode is run no need to trim sequences for alignment
            if barcode_only:
                pool = mp.Pool(self.cores)
        
                paths = pool.starmap_async(detect_barcodes, zip(self.input_file_paths, repeat(self))).get()
        
                pool.close()
                pool.join()
                
                barcode_paths = [item for item in paths]
        
                file_compiler(os.path.join("temp", "alignment_result.csv"), barcode_paths)
    
            else:
                pool = mp.Pool(self.cores)
        
                paths = pool.starmap_async(read_data, zip(self.input_file_paths, repeat(self))).get()
        
                pool.close()
                pool.join()
                
                trimmed_paths = [item for item in paths]
        
                file_compiler(os.path.join("temp", "trimmed_sequences.fastq"), trimmed_paths)
    
        
        def align_genomic_seq(self):
            """
            Aligns trimmed reads to reference using bowtie2. The output of bowtie2 
            is parsed and stored in alignment_result.csv.
            """
            set_up_tmpdir()
    
            # clean previous alignments
            align_res = os.path.join("temp", "alignment_result.csv")
            if os.path.exists(align_res):
                os.remove(align_res)
                
            run_bowtie2(self)
            parse_bowtie2_output(self)
        
        def write_count_matrix(self, barcode_only=False):
            """
            Assembles count matrix from bowtie2 alignment. Filters out spurious 
            barcodes with very little read counts for a given coordinate.
    
            """
            bowtie_res = read_csv(self.bowtie_res, 
                                   names=["id_in_pool","pool","coord","orientation","barcode", "ref"],
                                   dtype= {"id_in_pool": str, "pool": str, "coord":str, "orientation":str, "barcode":str, "ref":str})
                
            bowtie_res["barcode"] = bowtie_res["barcode"].astype(str)
        
            # barcodes that were not found are treated the same (can't distinguish between them)
            # concat orientation, coordinate and barcode to create unique ids
            bowtie_res["unique_id"] = bowtie_res["orientation"].astype(str) +";" +\
                bowtie_res["coord"].astype(str)+ ";"   + bowtie_res["barcode"].astype(str)+ ";"  + bowtie_res["ref"].astype(str)
            
            encoded_ids = get_hash_ids(bowtie_res["unique_id"])
            
            assert len(np.unique(encoded_ids)) == len(np.unique(bowtie_res["unique_id"])), "Hashed ids not unique!"
    
            count_mat = get_count_matrix(encoded_ids, bowtie_res, self)
            self.raw_count_mat = count_mat
            if barcode_only:
                self.count_mat = count_mat
                self.filtered_count_mat = count_mat
                count_mat.to_csv("count_matrix.csv", index=False, float_format="%.0f")
            else:
                filtered_count_mat = filter_barcodes(count_mat, self)
                filtered_count_mat.to_csv("count_matrix.csv", index=False, float_format="%.0f")
                
                self.count_mat = filtered_count_mat
    
            # normalize count matrix to counts per million
            normalized_count_mat = self.count_mat.copy()
            normalized_count_mat_ = normalized_count_mat.copy()
            counts = normalized_count_mat_[self.pools].values.astype(float)
            
            counts = (counts/np.sum(counts, axis=0))*1e6
            normalized_count_mat_[self.pools] = counts
            
            self.normalized_count_mat = normalized_count_mat_
            
            # apply count filters to reduce noise
            counts = global_filter_counts(counts, thr=self.global_filter_thr)
            counts = local_filter_counts(counts, thr=self.filter_thr)
            filtered_count_mat = normalized_count_mat_.copy()
            filtered_count_mat[self.pools] = counts
            self.filtered_count_mat = filtered_count_mat
            
        def deconvolve(self, barcode_only=False, count_mat="filtered"):
            """
            Deconvolve mutant count matrix and return summary output with 
            genes names.
            
            """
    
            unambiguous_data, ambiguous_data= get_ambiguity(self, count_mat)
            unambiguous_locations =  get_unambiguous_locations(unambiguous_data, self)
            ambiguous_locations = predict_ambiguous_locations(ambiguous_data, self)
            locations = concat([unambiguous_locations, ambiguous_locations])
            temppath = os.path.join("temp", "locations.csv")
            locations.to_csv(temppath)
           
            if barcode_only:
                locations.drop(["Reference", "Feature", "Orientation"], axis=1, inplace=True)
                locations.to_csv("barcode_location_summary.csv", index=False)
                
                transposed_barcode_summary = transpose_barcode_summary(self, locations)
                self.transposed_barcode_summary = transposed_barcode_summary
                transposed_barcode_summary.to_csv("well_barcode_summary.csv", index=False)
    
            else:
                location_summary = coords2genes(self, locations)
                self.location_summary = location_summary
                location_summary.to_csv("mutant_location_summary.csv", index=False)
                
                transposed_location_summary = transpose_location_summary(self, location_summary)
                self.transposed_location_summary = transposed_location_summary
                transposed_location_summary.to_csv("well_location_summary.csv", index=False)
            
        def deconvolve_validation(self, barcode_only=False):
            """
            Deconvolve mutant count matrix and return summary output with 
            genes names.
            
            """
            # normalize to library size and filter low reads
    
            
        
            unambiguous_data, ambiguous_data= get_ambiguity(self)
            unambiguous_locations =  get_unambiguous_locations(unambiguous_data, self)
            ambiguous_locations = predict_ambiguous_locations(ambiguous_data, self)
            locations = concat([unambiguous_locations, ambiguous_locations])
            temppath = os.path.join("temp", "locations.csv")
            locations.to_csv(temppath)
      

# experiment = LibraryExperiment(8,30,10,"gb_ref/UTI89.gb", "bowtie_ref/UTI89", 
#                                 "AGATGTGTATAAGAGACAG", 1, "input", "barcodes.csv",
#                                 True, "CGAGGTCTCT", "CGTACGCTGC")
# experiment.get_genomic_seq()        
# experiment.align_genomic_seq()        
# experiment.write_count_matrix()
# experiment.deconvolve()
