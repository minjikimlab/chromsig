def full_program(library_name, genome_name, fdr_thresh, chr_name, samp_size, bg_name, directory, file_name, strand_type):
    out_directory = directory + library_name + "_EnrichTest_FDR_" + str(fdr_thresh) + '/'
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    
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
        exp_enrich = random_loc_pe(subset_gems, chrom_size[chr_name], chr_name, samp_size, cov)
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
            pval, enr = get_raw_pval(exp_enrich[k], subset_gems[gem_id], samp_size, cov, strand_type)
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
