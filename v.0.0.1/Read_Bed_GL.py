import os
from sys import argv

if __name__ == '__main__':
    #argv = os.environ.get('NB_ARGS')
    argv = argv[1].split(" ")
    ### Set parameters ###
    directory = argv[0] ## Directory of bed + sizes files ##
    bed_file = argv[1] ## Name of the bed file ##
    size_file = argv[2] ## Name of the reference file ##

    std_bed_file = bed_file.replace(".bed", ".std.bed")

    if os.path.exists(directory + std_bed_file) == False:
        chrom_lst = []
        with open(directory + size_file, "r") as file:
            for line in file:
                chrom_lst.append(line.split("\t")[0])
        with open(directory + std_bed_file, "x") as dest, open(directory + bed_file, "r") as source:
            for line in source:
                if line.split("\t")[0] in chrom_lst:
                    dest.write(line)
