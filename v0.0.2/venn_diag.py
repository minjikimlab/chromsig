from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from sys import argv
import numpy as np

def getset(filename):
    return_lst = []
    with open(filename, "r") as file:
        file_size = sum(1 for _ in file)
        for line in file:
            line = line.split("\t")
            if line[0] == "chrM":
                continue
            return_lst.append(float(line[3]))
    return return_lst, file_size


if __name__ == '__main__':
    argv = argv[1].split(" ")
    onlya = argv[0]
    onlyb = argv[1]
    a_and_b = argv[2]
    output_dir = argv[3]
    lib_name = argv[4]

    seta, sizea = getset(onlya)
    setb, sizeb = getset(onlyb)
    setc, sizeab = getset(a_and_b)

    plt.figure(1)
    subsets = (sizea, sizeb, sizeab)
    venn2(subsets, set_labels=('Peaks in Total Bed', 'Peaks in Pass Pileup'))
    plt.title("Venn Diagram of Peaks Called by SICER")
    plt.savefig(output_dir + lib_name + '_venn_fig.png')


    plt.figure(2)
    print("Set A: " + str(seta))
    plt.boxplot([seta, setc, setb])
    plt.xticks([1, 2, 3], ['Total Peaks Values', 'Peaks Detected in Both Sets', 'Pass Pileup Peaks Values'])
    plt.ylabel('Confidence Values')
    plt.title('Box Plots of Confidence Values in Each Set of SICER Results')

    plt.savefig(output_dir + lib_name + '_box_fig.png')
