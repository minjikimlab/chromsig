import numpy as np
import random
import statsmodels.api as sm
import os
from statsmodels.sandbox.stats.multicomp import multipletests
import pybedtools
from pyBedGraph import BedGraph
from pybedtools import BedTool
from sys import argv
import pandas as pd
#import multiprocessing
#from memory_profiler import profile
import multiprocess as mp
#import defs_v2

import matplotlib.pyplot as plt
import time


# FDR: 0.1
# sample_size = 5000
# For 'chr1':
#   # Pass
#   # Fail
#   # Runtime/memory

# Read in samples from bed/bedgraph file
def read_bed(directory, file_name, chr_name):
    """
    Read a bed file.
    Args:
       directory (str): directory of the file location (ex: '/chromsig/Data_Directory/')
       file_name (str): name of the file (ex: 'GM12878_ChIP-seq_RAD51_ENCFF431RHO_bam2bed.std.bed')
       chr_name  (str): name of the current chromosome being processed

    Returns:
       subset_reads (dict): dictionary of reads in bed file, sorted by read ID
       cnt_all       (int): total number of fragments in the bed file
       cnt_subset    (int): number of fragments in the current chromosome
    """
    subset_reads = {}
    cnt_all = 0
    cnt_subset = 0
    search_str = chr_name + "\t"
    with open(directory + file_name) as f:
        for line in f:
            if search_str in line:
                tmp = line.strip().split("\t")
                subset_reads.setdefault(tmp[3].split("/")[0], []).append(tmp)
                cnt_subset += 1
            elif (search_str not in line) and (cnt_subset > 1):
                break
        cnt_all = len(f.readlines())
    read_span = int(subset_reads[list(subset_reads.keys())[0]][0][2]) - int(subset_reads[list(subset_reads.keys())[0]][0][1])
    return subset_reads, cnt_all, cnt_subset, read_span

    # Read the hg38.chrom.sizes file into a dictionary
def read_chroms(directory, genome_name):
    """
    Read a tab-delimited text file with list of chromosomes and their sizes.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       genome_name (str): name of reference genome (ex: 'dm3', 'mm10', 'hg38', 'hg19', etc.)
                        Note: must have a text file <genome_name>.chrom.sizes in the directory.

    Returns:
       chrom_dict (dictionary): dictionary with chromosomes as keys and sizes as values
    """
    chrom_dict = {}
    with open(directory + genome_name + '.chrom.sizes') as f:
        for line in f:
            tmp_list = line.strip().split("\t")
            chrom_dict[tmp_list[0]] = int(tmp_list[1])
    return chrom_dict


# Determine random start locations for samples to compare to observed sample
# Read span should be the same for all single-end fragments
def random_loc_se(chrom_size, read_span, chrom, sample_size, cov):
    """
    Determine all pseudo enrichment values for sample_size number of pseudo fragments.
    Args:
       chrom_size  (int): size of a given chromosome (ex: 23011544 if 'chr2L')
       read_span   (int): read from start to end of fragment (ex: 5000)
       sample_size (int): number of random locations to sample
       cov    (bedgraph): bedgraph file loaded in for the current chromosome by pyBedGraph module

    Returns:
       exp_enrich_lst (array): enrichment values for all pseudo reads
    """
    np.random.seed(12345)
    exp_enrich_lst = []
    i = 0
    while i < sample_size:
        startpos = np.random.randint(low = 0, high = chrom_size - read_span - 1)
        current_frag = [chrom, startpos, startpos + read_span]
        exp_enrich = cov.stats(stat = 'max', intervals = [current_frag])[0]
        if int(exp_enrich) == -1:
            i -= 1
        else:
            exp_enrich_lst.append(exp_enrich)
        i += 1
    return exp_enrich_lst


def random_loc_pe(subset_reads, chrom_size, chrom, sample_size, cov, read_span, out_directory, plot_hist):
    """
    Determine all pseudo enrichment values for sample_size number of pseudo fragments.
    Args:
       subset_reads (dict): dictionary of all reads for current chromosome
       chrom_size    (int): size of a the current chromosome being processed (ex: 23011544 if 'chr2L')
       chrom         (str): name of current chromosome (ex: 'chr10')
       sample_size   (int): number of random locations to sample
       cov      (bedgraph): bedgraph file loaded in for the current chromosome by pyBedGraph module
       read_span     (int): read from start to end of fragment (ex: 5000)

    Returns:
       exp_enrich_dct (dict): enrichment values for pseudo reads arranged by corresponding span
       spans_lst      (list): list of spans for all reads in subset_reads, arranged in order of reads
    """
    np.random.seed(12345)
    exp_enrich_dct = {}
    spans_lst = []
    all_tots = []
    print("Number of values total: " + str(len(subset_reads.values())))
    count = 0
    # Creating list of spans: spans for paired fragments are tuples
    # Each tuple contains four elements, with each element being the distance from the starting index of that read
    # Ex: [[‘chr1’, 10187, 10281], [‘chr1’, 10293, 10392]] --> (0, 94, 106, 205)
    for idx, value in enumerate(subset_reads.values()):
        total_span = int(value[-1][2]) - int(value[0][1])
        all_tots.append(total_span)
        all_dists = []
        start_idx = int(value[0][1])
        for elem in value:
            all_dists.append(int(elem[1]) - start_idx)
            all_dists.append(int(elem[2]) - start_idx)
        info_tup = tuple(all_dists)
        if len(info_tup) == 2:
            print("Current Value: " + str(value))
            print("Span: " + str(info_tup))
            count += 1
        spans_lst.append(info_tup)
    print("Number of Total Spans: " + str(len(spans_lst)))
    print("Number of Single End Values: " + str(count))
    # Creating a histogram of the spans
    if (plot_hist == "True"):
        hist_start = time.time()
        bins = np.arange(min(all_tots), max(all_tots) + 1, 1)
        print("Created Bins")
        plt.hist(all_tots, bins=bins, edgecolor='black', density=False)
        plt.xlabel("Distance between Paired Ends")
        plt.ylabel("Number of Paired Fragments with Corresponding Distance")
        plt.savefig(out_directory + '/plots/unique_dist_histogram_' + chrom + '.png')
        print("Made figure")
        hist_end = time.time()
        hist_time = hist_end - hist_start
        print("Time to create histogram: " + str(round(hist_time, 2)))
    # Creating a set of all spans to eliminate duplicates
    # This way, random locations are only generated for unique spans, preventing unnecessary repetition
    start_pos_dct = {}
    all_spans = set(spans_lst)
    for span in all_spans:
        start_pos_dct[span] = np.random.randint(low = 0, high = chrom_size - max(span) - 1, size = sample_size)
    for ind, span in enumerate(all_spans):
        current_frag = []
        print("Span: " + str(span))
        for pos in start_pos_dct[span]:
            current_frag.append([chrom, pos, pos + span[1]])
            current_frag.append([chrom, pos + span[2], pos + span[3]])
        exp_lst = list(cov.stats(stat = 'max', intervals = current_frag))
        # Replacing "faulty fragments": If a pseudo fragment happened to be located in an area of 
        # the bedgraph without any data, pyBedGraph returns a -1 value, which we want to avoid.
        # By sorting through the enrichments, we can determine which fragments resulted in -1, and
        # replace them with a fragment in an area with data.
        faulty_idxs = [index for index, enr_val in enumerate(exp_lst) if enr_val == -1]
        skip_next = False
        for badidx in faulty_idxs:
            if skip_next:
                skip_next = False
                continue
            int_divide = badidx // 2
            replace_idx = int_divide*2 
            new_exps = exp_lst[replace_idx:replace_idx + 2]
            while -1 in new_exps:
                startpos2 = np.random.randint(low = 0, high = chrom_size - max(span) - 1)
                better_frag = [[chrom, startpos2, startpos2 + span[1]], [chrom, startpos2 + span[2], startpos2 + span[3]]]
                new_exps = list(cov.stats(stat = 'max', intervals = better_frag))
            exp_lst[replace_idx] = new_exps[0]
            exp_lst[replace_idx + 1] = new_exps[1]
            if (badidx + 1) in faulty_idxs:
                skip_next = True
        # Arrange the enrichment lists into a dictionary organized by span
        exp_enrich_dct[span] = [sum(exp_lst[i:i + 2]) / 2 for i in range(0, len(exp_lst), 2)]
    return exp_enrich_dct, spans_lst



#def get_raw_pval(obs_frags, obs_span, obs_fraglen, obs_f2f, sample_size, cov, bin_size, chrom_size):
def get_raw_pval(exp_enrich, obs_frags, sample_size, cov, strand_type):
    """
    Compute raw p-value for a read.
    Args:
       exp_enrich (list): list of expected enrichments to compare the observed enrichment to
       obs_frags (list of list): [chrom,start,end] for each fragment
       sample_size (int): number of pseudo-reads to sample
       cov: bedgraph coverage track
       strand_type (str): Single End or Paired End (SE or PE)

    Returns:
       raw_pval (float): raw p-value computed
       obs_enrich (int): observed enrichment
    """
    # observed enrichment
    if str(strand_type) == "SE":
        obs_enrich = cov.stats(stat = 'max', intervals =  [[obs_frags[0][0], int(obs_frags[0][1]), int(obs_frags[0][2])]])[0]
    else:
        current_frag = []
        for elem in obs_frags:
            current_frag.append([elem[0], int(elem[1]), int(elem[2])])
        obs_lst = list(cov.stats(stat = 'max', intervals = current_frag))
        obs_enrich = sum(obs_lst)/len(obs_lst)
    tot = 0
    for i in exp_enrich:
        tot += (i >= obs_enrich)
    raw_pval = tot/sample_size
    return(raw_pval, obs_enrich)

def get_adj_pval(raw_pval_list, fdr_thresh, method):
    """ 
    Adjust raw pvalues by Benjamini Hochberg multiple testing adjustment. 
    Args:
       raw_pval_list (list): list of raw p-values (ex: [0.1, 0.04, 0.1])
       fdr_thresh (float): false discovery rate (ex: 0.05)
       method (string): adjustment method (ex: 'fdr_bh')

    Returns:
       adj_pval_list (array): array of booleans and adjusted p-values (ex: [0.1, 0.1, 0.1])
    """
    adj_pval_list = multipletests(raw_pval_list, alpha = fdr_thresh, method = method)
    return(adj_pval_list)

def write_master_result(out_read_list, out_name):
    """ 
    Write out significance test results.
    Args: 
       out_read_list (list): list with 11 items including read_id, frag_str, p-values, etc.
       out_name (string): output file name

    Returns:
       None
    """
    with open(out_name, 'a') as file1:
        header = ['READ_ID', 'Start Coord', 'End Coord', 'Obs', 'rawpval1', 'adjpval1', 'decis1']
        file1.write('\t'.join(map(str, header)) + '\n')
        for i in range(len(out_read_list)):
            file1.write('\t'.join(map(str, out_read_list[i])) + '\n')

def write_output_beds(lst, out_name):
    with open(out_name, 'a') as file:
        for line in lst:
            file.write('\t'.join(line) + '\n')



if __name__ == '__main__':
    total_start = time.time()
    argv = argv[1].split(" ")
    ### Set parameters ###
    library_name = argv[0] ## Library name of our data ##
    genome_name = argv[1] ## Name of the reference genome ##
    fdr_thresh = float(argv[2])  # should be argument; Benjamini-Hochberg FDR; p-value cutoff ##
    chr_name = argv[3] # should be argument
    samp_size = int(argv[4]) ## Number of pseudo-reads ##
    bg_name = argv[5] # bedgraph file name
    directory = argv[6]
    file_name = argv[7]
    strand_type = argv[8] ## Paired End or Single End (PE or SE) ##
    plot_hist = argv[9] ## Paired End or Single End (PE or SE) ##
    print("Read in Args\n")

    out_directory = directory + library_name + "_EnrichTest_FDR_" + str(fdr_thresh) + "_pseudoRead_" + str(samp_size) + '/'

    prefix = library_name + "_" + chr_name + "_FDR_" + str(fdr_thresh) + "_pseudoRead_" + str(samp_size) + str(strand_type)

    ### Set directory and input file name ###

    #### Log file ####
    out = open(out_directory + prefix + "_enrichTest_logFile.txt", "a")

    out.write("Software version: v0.0.1 (2025-07-01, Gupta)" + "\n")
    out.write("Input directory: " + directory + "\n")
    out.write("Input file name: " + file_name + "\n")
    out.write("Input bedgraph file: " + bg_name + "\n")
    out.write("Output directory: " + out_directory + "\n")
    out.write("Library name: " + library_name + "\n")
    out.write("Reference genome: " + genome_name + "\n")
    out.write("FDR threshold: " + str(fdr_thresh) + "\n")
    out.write("Number of pseudo-Reads: " + str(samp_size) + "\n")
    out.write("Single or Paired End: " + str(strand_type) + "\n")
    out.write("Started processing frequency-based enrichment test. \n")
    out.write("================================= \n")
    out.write("===== Chromosome is: " + chr_name + ". ===== \n")
    out.write("================================= \n")

    read_bed_start = time.time()
    subset_reads, cnt_all, cnt_subset, read_span = read_bed(directory, file_name, chr_name)
    out.write("Finished reading the input Read file. \n")
    out.write(str(cnt_all-1) + " total Reads. \n")
    out.write("================================= \n")
    read_bed_end = time.time()

    print("Read Bed File\n")

    ### Read chrom sizes file ###
    read_sizes_start = time.time()
    chrom_size = read_chroms(directory, genome_name)
    read_sizes_end = time.time()
    ### Subset input Read files by chromosome ###
    #subset_reads = [x for x in in_reads if chr_name == x[1].split(":")[0]]
    #out.write("Finished subsetting Reads. \n")
    out.write(str(cnt_subset) + " Reads in " + chr_name + " (of " + str(cnt_all-1) + " total Reads). \n")
    out.write("================================= \n")
    print("Read Chroms Sizes File\n")
    
    #del in_reads
    
    ### Read bedgraph coverage file ###
    #bg_name = library_name + '_binned_' + str(bin_size) + 'bp_' + chr_name + '.bedgraph'
    if not os.path.isfile(directory+bg_name):
        out.write("Error: no bedgraph coverage file." + "\n")
    #cov = BedGraph(directory + genome_name+".chrom.sizes", directory + bg_name, str(chr_name))
    pyBedGraph_start = time.time()
    cov = BedGraph(directory + genome_name+".chrom.sizes", directory + bg_name, [chr_name])
    cov.load_chrom_data(chr_name)
    #coverage = read_bedgraph(directory, library_name, bin_size, chr_name)
    out.write("Finished reading the coverage file " + bg_name + ". \n")
    out.write("================================= \n")
    pyBedGraph_end = time.time()
    print("Loaded PyBedGraph\n")

    ### Initialize output Read list ###
    out_reads = []
    
    ### Calculate raw p-values for subset_reads ###
    pvals = []

    random_loc_start = time.time()
    if str(strand_type) == "SE":
        exp_enrich = random_loc_se(chrom_size[chr_name], read_span, chr_name, samp_size, cov)
    else:
        #print("Starting random location")
        exp_enrich, spans_lst = random_loc_pe(subset_reads, chrom_size[chr_name], chr_name, samp_size, cov, read_span, out_directory, plot_hist)
    random_loc_end = time.time()
    #out.write(str(exp_enrich[list(subset_reads.keys())[0]]))
    #exp_enrich = []
    #pseudo = []
    
    #out.write("Pseudo: " + str(pseudo[i]) + "\n")
    #out.write("Exp_enrich: " + str(exp_enrich) + "\n")
        
    #exp_enrich.append(sum(exp)/len(exp))
    out.write("Finished generating pseudo samples.\n")
    out.write("================================= \n")
    raw_pvals_start = time.time()
    read_ids = list(subset_reads.keys())
    for k in range(len(subset_reads)):
        read_id = read_ids[k]
        #read_id, chrom, frags, span, fraglen, f2fdist, fragnum, frag_str, coord = list(extract_info(subset_reads[k]))
        if str(strand_type) == "SE":
            #out.write(str(subset_reads[read_id]) + "\n")
            pval, enr = get_raw_pval(exp_enrich, subset_reads[read_id], samp_size, cov, strand_type)
        else:
            pval, enr = get_raw_pval(exp_enrich[spans_lst[k]], subset_reads[read_id], samp_size, cov, strand_type)
        pvals.append(pval)
        #out.write("Pval: " + str(pval) + "\n")
        out_reads.append([read_id, subset_reads[read_id][0][1], subset_reads[read_id][-1][2], round(enr,1), round(pval,3)])
    out.write("Finished calculating raw p-values for " + str(k + 1) + " Reads. \n")
    #raw_pass = sum(i <= fdr_thresh for i in pvals)
    raw_pass = 0
    for i in pvals:
        raw_pass += (i <= fdr_thresh)
    out.write(str(raw_pass) + " Reads have raw p-val <= " + str(fdr_thresh) + "(i.e., " + "{0:.2f}".format(raw_pass/(k+1)*100) +" %). \n")
    out.write("================================= \n")
    raw_pvals_end = time.time()
    
    ### Adjust p-values ###
    adj_pvals_start = time.time()
    pass_pileup = []
    fail_pileup = []
    adj_pval = get_adj_pval(pvals, fdr_thresh, 'fdr_bh')
    out.write("Finished adjusting p-values for " + str(k + 1) + " Reads. \n")
    adj_pass = 0
    full_elem_lst = list(subset_reads.values())
    for i in range(len(adj_pval[1])):
        #print(adj_pval[0][i])
        out_reads[i].append(round(adj_pval[1][i], 3))
        elem_lst = full_elem_lst[i]
        if adj_pval[0][i]:
            out_reads[i].extend(['PASS'])
            pass_pileup.extend(elem_lst)
            adj_pass += 1
        else:
            out_reads[i].extend(['FAIL'])
            fail_pileup.extend(elem_lst)
    out.write("Pass Pileup: " + str(len(pass_pileup)) + "\n")
    out.write("Fail Pileup: " + str(len(fail_pileup)) + "\n")
    out.write(str(adj_pass) + " Reads have adjusted p-val <= " + str(fdr_thresh) + "(i.e.," + "{0:.2f}".format(adj_pass/(k+1)*100) +" %). \n")
    out.write("================================= \n")
    adj_pvals_end = time.time()
    
    ### Write results ###
    write_res_start = time.time()
    write_master_result(out_reads, out_directory+prefix+'_enrichTest_master.txt')
    write_output_beds(pass_pileup, out_directory+prefix+'_pass_pileup.bed')
    write_output_beds(fail_pileup, out_directory+prefix+'_fail_pileup.bed')
    write_res_end = time.time()
    out.write("Finished writing files. \n")
    out.write("DONE. \n")
    total_end = time.time()
    out.write("================================= \n")
    out.write("Timings: \n")
    out.write("Reading Bed File: " + str(round(read_bed_end - read_bed_start, 2)) + "\n")
    out.write("Reading Sizes: " + str(round(read_sizes_end - read_sizes_start, 2)) + "\n")
    out.write("pyBedGraph: " + str(round(pyBedGraph_end - pyBedGraph_start, 2)) + "\n")
    out.write("Random Locations: " + str(round(random_loc_end - random_loc_start, 2)) + "\n")
    out.write("Raw Pvals: " + str(round(raw_pvals_end - raw_pvals_start, 2)) + "\n")
    out.write("Adjusted Pvals: " + str(round(adj_pvals_end - adj_pvals_start, 2)) + "\n")
    out.write("Writing Results: " + str(round(write_res_end - write_res_start, 2)) + "\n")
    out.write("Total: " + str(round(total_end - total_start, 2)) + "\n")
    out.write("================================= \n")
    #out.write("Memory Usage: \n")
    #out.write(str(tracemalloc.get_traced_memory()))
    #tracemalloc.stop()
    out.close()
