import numpy as np
import random
import statsmodels.api as sm
import os
from statsmodels.sandbox.stats.multicomp import multipletests
import pybedtools
from pyBedGraph import BedGraph
from pybedtools import BedTool
from sys import argv
import multiprocessing
from memory_profiler import profile
import multiprocess as mp
#import defs_v2

import time



# instantiating the decorator
@profile

# FDR: 0.1
# sample_size = 5000
# For 'chr1':
#   # Pass
#   # Fail
#   # Runtime/memory

# Read in samples from bed/bedgraph file
def read_gems(directory, file_name, chr_name):
    """
    Read a GEM file of the form "PlinePgem".
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       file_name (str): name of the file (ex: 'SHG0008H.Fragnum_PlinePgem')
    Returns:
       in_gems (list): tab-separated lists of lists
    """
    subset_gems = {}
    cnt_all = 0
    cnt_subset = 0
    with open(directory + file_name) as f:
        for line in f:
            tmp = line.strip().split("\t")
            if tmp[0] == chr_name:
                subset_gems.setdefault(tmp[3].split("/")[0], []).append(tmp)
                cnt_subset += 1
            cnt_all += 1
        #in_gems = [line.strip().split("\t") for line in f]
    gem_span = int(subset_gems[list(subset_gems.keys())[0]][0][2]) - int(subset_gems[list(subset_gems.keys())[0]][0][1])
    return subset_gems, cnt_all, cnt_subset, gem_span

    # Read the hg38.chrom.sizes file into a dictionary
def read_chroms(directory, genome_name):
    """
    Read a tab-delimited text file with list of chromosomes and their sizes.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       genome_name (str): name of reference genome (ex: 'dm3', 'mm10', 'hg38', 'hg19', etc.)
                        Note: must have a text file <genome_name>.chrom.sizes in the directory.
    Returns:
       chrom_dict (dictionary): tab-separated lists of lists
    """
    chrom_dict = {}
    with open(directory + genome_name + '.chrom.sizes') as f:
        for line in f:
            tmp_list = line.strip().split("\t")
            chrom_dict[tmp_list[0]] = int(tmp_list[1])
    return chrom_dict


# Determine random start locations for samples to compare to observed sample
# Gem span should be the same for all single-end fragments
def random_loc_se(chrom_size, gem_span, chrom, sample_size, cov):
    """
    Randomly locate GEM with same span in a chromosome sample_size times.
    Args:
       chrom_size (int): size of a given chromosome (ex: 23011544 if 'chr2L')
       gem_span (int): GEM from start of first fragment to end of last fragment (ex: 5000)
       sample_size (int): number of random locations to sample
    Returns:
       startpos (array): random integers array of length 'sample_size'
    """
    np.random.seed(12345)
    exp_enrich_lst = []
    i = 0
    while i < sample_size:
        startpos = np.random.randint(low = 0, high = chrom_size - gem_span - 1)
        current_frag = [chrom, startpos, startpos + gem_span]
        exp_enrich = cov.stats(stat = 'max', intervals = [current_frag])[0]
        #outfile.write("CurrentIndex: " + str(i) + "\t")
        #outfile.write("Exp Enrichment: " + str(exp_enrich) + "\n")
        #outfile.write("Current Pseudo: " + str(pseudo) + "\n")
        if int(exp_enrich) == -1:
            i -= 1
        else:
            exp_enrich_lst.append(exp_enrich)
        #outfile.write("CurrentIndex SameCycle: " + str(i) + "\t")
        #outfile.write("Exp Enrichment SameCycle: " + str(exp_enrich) + "\n")
        #outfile.write("Current Pseudo SameCycle: " + str(pseudo) + "\n")
        i += 1
    return exp_enrich_lst
    
'''
def random_loc_pe(subset_gems, chrom_size, chrom, sample_size, cov):
    print("In old random_loc\n")
    np.random.seed(12345)
    exp_enrich_lst = []
    for value in subset_gems.values():
        i = 0
        start_idx = int(value[0][1])
        total_span = int(value[-1][2]) - start_idx
        one_samp_exp_enrich = []
        while i < sample_size:
            startpos = np.random.randint(low = 0, high = chrom_size - total_span - 1)
            current_frag = []
            for elem in value:
                current_frag.append([chrom, int(startpos + int(elem[1])-start_idx), int(startpos + int(elem[2])-start_idx)])
            exp_lst = list(cov.stats(stat = 'max', intervals = current_frag))
            exp_enrich = sum(exp_lst)/len(exp_lst)
            if int(exp_enrich) == -1:
                i -= 1
            else:
                #exp_enrich_dct.setdefault(key, []).append(exp_enrich)
                one_samp_exp_enrich.append(exp_enrich)
            i += 1
        exp_enrich_lst.append(one_samp_exp_enrich)
    return exp_enrich_lst
'''
'''
def random_loc_pe(subset_gems, chrom_size, chrom, sample_size, cov, out):
    print("In new random_loc_pe")
    out.write("In new random_loc_pe\n")
    np.random.seed(12345)
    exp_enrich_lst = []
    tot_pseudo_frag = []
    for value in subset_gems.values():
        start_idx = int(value[0][1])
        total_span = int(value[-1][2]) - start_idx
        startpos = np.random.randint(low = 0, high = chrom_size - total_span - 1, size = sample_size)
        for pos in startpos:
            pseudo_frag = []
            for elem in value:
                pseudo_frag.append([chrom, pos + int(elem[1]) - start_idx, pos + int(elem[2]) - start_idx])
            tot_pseudo_frag.append(pseudo_frag)
            #out.write(str(tot_pseudo_frag))
    print("Pseudos created")
    out.write("Pseudos created\n")
    one_d_frags = sum(tot_pseudo_frag, [])
    exp_lst = list(cov.stats(stat = 'max', intervals = one_d_frags))
    out.write("Enrichments found\n")
    exp_lst_arr = []
    # Rearrange exp_lst to keep maximums for the same value grouped together
    for validx, value in enumerate(subset_gems.values()):
        avg_len = len(value)
        exp_lst_arr.append([exp_lst[i:i + avg_len] for i in range(validx*sample_size, validx*sample_size + sample_size, avg_len)])
    out.write("Enrichments Arranged: " + str(exp_lst_arr[0]) + "\n")
    ## All occurences of -1 in exp_lst:
    faulty_idxs = [index for index, enr_val in enumerate(exp_lst_arr) if (-1 in enr_val)]
    for badidx in faulty_idxs:
        out.write("Bad: " + str(exp_lst_arr[badidx]) + "\t")
        while -1 in exp_lst_arr[badidx]:
            badfrag = tot_pseudo_frag[badidx]
            start_idx = int(badfrag[0][1])
            total_span = int(badfrag[-1][2]) - start_idx
            startpos2 = np.random.randint(low = 0, high = chrom_size - total_span - 1) - start_idx
            better_frag = []
            for elem in tot_pseudo_frag[badfrag]:
                better_frag.append([chrom, startpos2 + int(elem[1]), startpos2 + int(elem[2])])
            exp_lst_arr[badidx] = list(cov.stats(stat = 'max', intervals = better_frag))
        out.write("Better: " + str(exp_lst_arr[badidx]) + "\n")
    for val in range(0, len(exp_lst_arr), sample_size):
        exp_enrich_lst.append([sum(exp_lst_arr[i]) / len(exp_lst_arr[i]) for i in range(val, val + sample_size)])
    #exp_enrich = sum(exp_lst)/len(exp_lst)
    exp_enrich_lst = 0
    return exp_enrich_lst
'''


def random_loc_pe(subset_gems, chrom_size, chrom, sample_size, cov, gem_span):
    np.random.seed(12345)
    exp_enrich_dct = {}
    spans_lst = []
    for value in subset_gems.values():
        all_dists = []
        start_idx = int(value[0][1])
        for elem in value:
            all_dists.append(int(elem[1]) - start_idx)
        info_tup = (len(value), tuple(all_dists))
        spans_lst.append(info_tup)
    all_spans = list(set(spans_lst))
    start_pos_dct = {}
    for span in all_spans:
        start_pos_dct[span] = np.random.randint(low = 0, high = chrom_size - (span[1][-1] + gem_span) - 1, size = sample_size)
    for span in all_spans:
        current_frag = []
        for pos in start_pos_dct[span]:
            for dist_idx in range(span[0]):
                current_frag.append([chrom, pos + span[1][dist_idx], pos + span[1][dist_idx] + gem_span])
        exp_lst = list(cov.stats(stat = 'max', intervals = current_frag))
        ## All occurences of -1 in exp_lst:
        avg_len = span[0]
        faulty_idxs = [index for index, enr_val in enumerate(exp_lst) if enr_val == -1]
        for badidx in faulty_idxs:
            replace_idx = badidx // avg_len
            new_exps = exp_lst[replace_idx:replace_idx + avg_len]
            while -1 in new_exps:
                startpos2 = np.random.randint(low = 0, high = chrom_size - (span[1][-1] + gem_span) - 1)
                better_frag = []
                for dist_idx in range(span[0]):
                    better_frag.append([chrom, startpos2 + span[1][dist_idx], startpos2 + span[1][dist_idx] + gem_span])
                new_exps = list(cov.stats(stat = 'max', intervals = better_frag))
            for i in range(0, avg_len):
                exp_lst[replace_idx + i] = new_exps[i]
        exp_enrich_dct[span] = [sum(exp_lst[i:i + avg_len]) / avg_len for i in range(0, len(exp_lst), avg_len)]
        #exp_enrich = sum(exp_lst)/len(exp_lst)
    return exp_enrich_dct, spans_lst



#def get_raw_pval(obs_frags, obs_span, obs_fraglen, obs_f2f, sample_size, cov, bin_size, chrom_size):
def get_raw_pval(exp_enrich, obs_frags, sample_size, cov, strand_type):
    """
    Compute raw p-value for a GEM.
    Args:
       obs_frags (list of list): [chrom,start,end] for each fragment 
       obs_span (int): start of leftmost fragment to end of rightmost fragment
       obs_fraglen (list): length of each fragment
       obs_f2f (list): distance between neighboring fragments
       sample_size (int): number of pseudo-GEMs to sample
       cov: bedgraph coverage track
       bin_size (int): size of the bin used to generate cov
       chrom_size (int): size of the chromosome
    Returns:
       raw_pval (float): raw p-value computed
       obs_enrich (int): observed enrichment
    """
    # observed enrichment
#    frag_cov = []
#    for i in range(len(obs_frags)):
#        frag_cov.append(get_mean_cov(obs_frags[i][0], obs_frags[i][1], cov, bin_size))
    #frag_cov = list(cov.stats(stat = 'mean', intervals = [[obs_frags[0], int(obs_frags[1]), int(obs_frags[2])]]))
    #obs_enrich = cov.stats(stat = 'max', intervals = [[obs_frags[0], int(obs_frags[1]), int(obs_frags[2])]])[0]
    if str(strand_type) == "SE":
        obs_enrich = cov.stats(stat = 'max', intervals =  [[obs_frags[0][0], int(obs_frags[0][1]), int(obs_frags[0][2])]])[0]
    else:
        current_frag = []
        for elem in obs_frags:
            current_frag.append([elem[0], int(elem[1]), int(elem[2])])
        obs_lst = list(cov.stats(stat = 'max', intervals = current_frag))
        #new_lst = [obs_frags[:, 0].tolist(), obs_frags[:, 1].tolist(), obs_frags[:, 2].astype(int).tolist()]
        #new_lst = [obs_frags[:, 0].tolist(), obs_frags[:, 1].tolist(), obs_frags[:, 2].tolist()]
        #obs_lst = list(cov.stats(stat = 'max', intervals =  [[row[i] for row in new_lst] for i in range(len(new_lst[0]))]))
        obs_enrich = sum(obs_lst)/len(obs_lst)
        #obs_enrich = max(frag_cov)
    # expected enrichment (sampling background)
    #startpos = random_loc(chrom_size, obs_span, sample_size)
    '''
    exp_enrich = []
    for k in range(sample_size):
        pseudo = pseudo_gem(startpos[k], obs_span, chr_name)
        exp = list(cov.stats(stat = 'mean', intervals = pseudo))
        exp_enrich.append(sum(exp)/len(exp))
        #exp_enrich.append(max(exp))
    '''
    # write 1000 null (expected, sampled) values; 20191226
#    prefix_null = "_".join(prefix.split("_")[:-1]) + "_1000"
#    with open(out_directory + prefix_null + "_enrichTest_null.txt", "a") as f1:
#        f1.write(','.join(map(str, [round(x,1) for x in exp_enrich])) + '\n')
#    f1.close()
    # compute raw p-value
    #raw_pval = sum(i > obs_enrich for i in exp_enrich)/sample_size
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

def write_master_result(out_gem_list, out_name):
    """ 
    Write out significance test results.
    Args: 
       out_gem_list (list): list with 11 items including gem_id, frag_str, p-values, etc.
       out_name (string): output file name
    Returns:
       None
    """
    with open(out_name, 'a') as file1:
        header = ['GEM_ID', 'Start Coord', 'End Coord', 'Category', 'Obs', 'rawpval1', 'adjpval1', 'decis1', 'rawpval2', 'adjpval2', 'decis2']
        file1.write('\t'.join(map(str, header)) + '\n')
        for i in range(len(out_gem_list)):
            file1.write('\t'.join(map(str, out_gem_list[i])) + '\n')

def write_output_beds(lst, out_name):
    with open(out_name, 'a') as file:
        for line in lst:
            file.write('\t'.join(line) + '\n')





def full_program(library_name, genome_name, fdr_thresh, chr_name, samp_size, bg_name, directory, file_name, strand_type):
    
    prefix = library_name + "_" + chr_name + "_FDR_" + str(fdr_thresh) + "_pseudoGEM_" + str(samp_size) + str(strand_type)

    ### Set directory and input file name ###
    #gem_span = int(argv[8])

    #### Log file ####
    out = open(out_directory + prefix + "_enrichTest_logFile.txt", "a")

    #out.write("Software version: v1.0 (2020-04-01, Kim)" + "\n")
    out.write("Input directory: " + directory + "\n")
    out.write("Input file name: " + file_name + "\n")
    out.write("Input bedgraph file: " + bg_name + "\n")
    out.write("Output directory: " + out_directory + "\n")
    out.write("Library name: " + library_name + "\n")
    out.write("Reference genome: " + genome_name + "\n")
    out.write("FDR threshold: " + str(fdr_thresh) + "\n")
    out.write("Number of pseudo-GEMs: " + str(samp_size) + "\n")
    out.write("Single or Paired End: " + str(strand_type) + "\n")
    out.write("Started processing frequency-based enrichment test. \n")
    out.write("================================= \n")
    out.write("===== Chromosome is: " + chr_name + ". ===== \n")
    out.write("================================= \n")

    read_bed_start = time.time()
    subset_gems, cnt_all, cnt_subset, gem_span = read_gems(directory, file_name, chr_name)
    out.write("Finished reading the input GEM file. \n")
    out.write(str(cnt_all-1) + " total GEMs. \n")
    out.write("================================= \n")
    read_bed_end = time.time()

    ### Read chrom sizes file ###
    read_sizes_start = time.time()
    chrom_size = read_chroms(directory, genome_name)
    read_sizes_end = time.time()
    ### Subset input GEM files by chromosome ###
    #subset_gems = [x for x in in_gems if chr_name == x[1].split(":")[0]]
    #out.write("Finished subsetting GEMs. \n")
    out.write(str(cnt_subset) + " GEMs in " + chr_name + " (of " + str(cnt_all-1) + " total GEMs). \n")
    out.write("================================= \n")
    
    #del in_gems
    
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

    ### Initialize output GEM list ###
    out_gems = []
    
    ### Calculate raw p-values for subset_gems ###
    pvals = []

    random_loc_start = time.time()
    if str(strand_type) == "SE":
        exp_enrich = random_loc_se(chrom_size[chr_name], gem_span, chr_name, samp_size, cov)
    else:
        exp_enrich, spans_lst = random_loc_pe(subset_gems, chrom_size[chr_name], chr_name, samp_size, cov, gem_span)
    random_loc_end = time.time()
    #out.write(str(exp_enrich[list(subset_gems.keys())[0]]))
    #exp_enrich = []
    #pseudo = []
    
    #out.write("Pseudo: " + str(pseudo[i]) + "\n")
    #out.write("Exp_enrich: " + str(exp_enrich) + "\n")
        
    #exp_enrich.append(sum(exp)/len(exp))
    out.write("Finished generating pseudo samples.\n")
    out.write("================================= \n")
    raw_pvals_start = time.time()
    gem_ids = list(subset_gems.keys())
    for k in range(len(subset_gems)):
        gem_id = gem_ids[k]
        #gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, frag_str, coord = list(extract_info(subset_gems[k]))
        if str(strand_type) == "SE":
            #out.write(str(subset_gems[gem_id]) + "\n")
            pval, enr = get_raw_pval(exp_enrich, subset_gems[gem_id], samp_size, cov, strand_type)
        else:
            pval, enr = get_raw_pval(exp_enrich[spans_lst[k]], subset_gems[gem_id], samp_size, cov, strand_type)
        pvals.append(pval)
        #out.write("Pval: " + str(pval) + "\n")
        out_gems.append([gem_id, subset_gems[gem_id][0][1], subset_gems[gem_id][-1][2], 'Orig', round(enr,1), round(pval,3)])
    out.write("Finished calculating raw p-values for " + str(k + 1) + " GEMs. \n")
    #raw_pass = sum(i <= fdr_thresh for i in pvals)
    raw_pass = 0
    for i in pvals:
        raw_pass += (i <= fdr_thresh)
    out.write(str(raw_pass) + " GEMs have raw p-val <= " + str(fdr_thresh) + "(i.e., " + "{0:.2f}".format(raw_pass/(k+1)*100) +" %). \n")
    out.write("================================= \n")
    raw_pvals_end = time.time()
    
    ### Adjust p-values ###
    adj_pvals_start = time.time()
    out.write(str(adj_pvals_start) + "\n")
    pass_pileup = []
    fail_pileup = []
    adj_pval = get_adj_pval(pvals, fdr_thresh, 'fdr_bh')
    out.write("Finished adjusting p-values for " + str(k + 1) + " GEMs. \n")
    adj_pass = 0
    full_elem_lst = list(subset_gems.values())
    for i in range(len(adj_pval[1])):
        #print(adj_pval[0][i])
        out_gems[i].append(round(adj_pval[1][i], 3))
        elem_lst = full_elem_lst[i]
        if adj_pval[0][i]:
            out_gems[i].extend(['PASS','.','.','.'])
            pass_pileup.extend(elem_lst)
            adj_pass += 1
        else:
            out_gems[i].extend(['FAIL','.','.','.'])
            fail_pileup.extend(elem_lst)
    out.write("Pass Pileup: " + str(len(pass_pileup)) + "\n")
    out.write("Fail Pileup: " + str(len(fail_pileup)) + "\n")
    out.write(str(adj_pass) + " GEMs have adjusted p-val <= " + str(fdr_thresh) + "(i.e.," + "{0:.2f}".format(adj_pass/(k+1)*100) +" %). \n")
    out.write("================================= \n")
    adj_pvals_end = time.time()
    out.write(str(adj_pvals_end))
    
    ### Write results ###
    write_res_start = time.time()
    write_master_result(out_gems, out_directory+prefix+'_enrichTest_master.txt')
    write_output_beds(pass_pileup, out_directory+prefix+'_pass_pileup.bed')
    write_output_beds(fail_pileup, out_directory+prefix+'_fail_pileup.bed')
    write_res_end = time.time()
    out.write("Finished writing files. \n")
    out.write("DONE. \n")
    total_end = time.time()
    out.write("================================= \n")
    out.write("Timings: \n")
    out.write("Reading Bed File: " + str(read_bed_end - read_bed_start) + "\n")
    out.write("Reading Sizes: " + str(read_sizes_end - read_sizes_start) + "\n")
    out.write("pyBedGraph: " + str(pyBedGraph_end - pyBedGraph_start) + "\n")
    out.write("Random Locations: " + str(random_loc_end - random_loc_start) + "\n")
    out.write("Raw Pvals: " + str(raw_pvals_end - raw_pvals_start) + "\n")
    out.write("Adjusted Pvals: " + str(adj_pvals_end - adj_pvals_start) + "\n")
    out.write("Writing Results: " + str(write_res_end - write_res_start) + "\n")
    out.write("Total: " + str(total_end - total_start) + "\n")
    out.write("================================= \n")
    #out.write("Memory Usage: \n")
    #out.write(str(tracemalloc.get_traced_memory()))
    #tracemalloc.stop()
    out.close()


if __name__ == '__main__':
    #tracemalloc.start()
    total_start = time.time()
    #argv = os.environ.get('NB_ARGS')
    #print(str(argv))
    argv = argv[1].split(" ")
    print(argv)
    ### Set parameters ###
    library_name = argv[0] ## Library name of our data ##
    genome_name = argv[1] ## Name of the reference genome ##
    fdr_thresh = float(argv[2])  # should be argument; Benjamini-Hochberg FDR; p-value cutoff ##
    chr_lst = argv[3] # should be argument
    samp_size = int(argv[4]) ## Number of pseudo-GEMs ##
#    bin_size = int(argv[6])
    bg_name = argv[5] # bedgraph file name
    directory = argv[6]
    file_name = argv[7]
    strand_type = argv[8] ## Paired End or Single End (PE or SE) ##
    #out_directory = directory + library_name + "_EnrichTest_FDR_" + str(fdr_thresh) + '/'
    #if not os.path.exists(out_directory):
    #    os.mkdir(out_directory)

    out_directory = directory + library_name + "_EnrichTest_FDR_" + str(fdr_thresh) + "_pseudoGEM_" + str(samp_size) + '/'
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    #print(chr_lst)
    chr_lst = chr_lst.split("\n")
    #print(chr_lst)

    #for chr_name in chr_lst:
    #    prefix = library_name + "_" + chr_name + "_FDR_" + str(fdr_thresh) + "_pseudoGEM_" + str(samp_size) + str(strand_type)
    args_lst = []
    for chr_name in chr_lst:
        args_lst.append((library_name, genome_name, fdr_thresh, chr_name, samp_size, bg_name, directory, file_name, strand_type))
    max_parallel = mp.cpu_count()
        ### Set directory and input file name ###
        #gem_span = int(argv[8])
    
        #### Log file ####
        #print("Current args: " + str(elem))
    with mp.Pool(processes=1) as pool:
        pool.starmap(full_program, args_lst)