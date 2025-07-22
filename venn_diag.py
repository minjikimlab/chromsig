from matplotlib_venn import venn2
import matplotlib.pyplot as plt
from sys import argv
import numpy as np

def getset(filename):
    return_set = ""
    return_lst = []
    total_set = {}
    with open(filename, "r") as file:
        for line in file:
            line = line.split("\t")
            return_set = return_set + line[0] + "\t" + line[1] + "\t" + line[2] + "\t"
            return_lst.append(line[0] + "\t" + line[1] + "\t" + line[2])
            total_set[return_lst] = float(line[3])
    return_set = return_set.rstrip()
    return return_set, return_lst, total_set


def convert_to_set(start_str):
    lst = start_str.split("\t")

    two_d_list = []
    for i in range(0, len(lst), 3):
        two_d_list.append(tuple(lst[i:i + 3]))
    final_set = set(two_d_list)
    return final_set


if __name__ == '__main__':
    argv = argv[1].split(" ")
    filea = argv[0]
    fileb = argv[1]

    seta, return_lsta, total_seta = getset(filea)
    setb, return_lstb, total_setb = getset(fileb)

    set1 = convert_to_set(seta)
    set2 = convert_to_set(setb)

    expanded_lsta = []
    for line in return_lsta:
        line = line.split("\t")
        for i in range(int(line[1]), int(line[2])):
            expanded_lsta.append(line[0] + "\t" + str(int(line[1]) + i) + "\t" + str(int(line[1]) + i))
    expanded_lstb = []
    for line in return_lstb:
        line = line.split("\t")
        for i in range(int(line[1]), int(line[2])):
            expanded_lstb.append(line[0] + "\t" + str(int(line[1]) + i) + "\t" + str(int(line[1]) + i))


    fig, axs = plt.subplots(2, 1, figsize=(5, 10))

    venn2([set(expanded_lsta), set(expanded_lstb)], set_labels=('Total Bedgraph', 'Pass Pileup Bedgraph'), ax=axs[0])
    axs[0].set_title("SICER Peak Calling Results Pre and Post Denoising")

    return_lsta = set(expanded_lsta)
    return_lstb = set(expanded_lstb)

    onlya = return_lsta - return_lstb
    onlyb = return_lstb - return_lsta
    intersectionab = return_lsta.intersection(return_lstb)

    # Find values for only a
    avals = []
    curr_file = open("/home/njgupta/chromsig/list_a.txt", "a")
    count = 0
    for elem in onlya:
        if elem in total_seta.keys():
            avals.append(total_seta[elem])
            count += 1
    curr_file.write("Matches made: " + str(count) + " Total count: " + str(len(total_seta)))
    curr_file.close()

    # Find values for only b
    bvals = []
    curr_file = open("/home/njgupta/chromsig/list_b.txt", "a")
    count = 0
    for elem in onlyb:
        if elem in total_setb.keys():
            bvals.append(total_seta[elem])
            count += 1
    curr_file.write("Matches made: " + str(count) + " Total count: " + str(len(total_seta)))
    curr_file.close()
    
    # Find values for both a and b
    abvals = []
    curr_file = open("/home/njgupta/chromsig/list_ab.txt", "a")
    count = 0
    for elem in intersectionab:
        if elem in total_seta.keys():
            abvals.append(total_seta[elem])
            count += 1
    curr_file.write("Matches made: " + str(count) + " Total count: " + str(len(total_seta)))
    curr_file.close()

    box_plot_data = [avals, abvals, bvals]

    axs[1].boxplot(box_plot_data)
    axs[1].set_xticks([1, 2, 3], ['Peaks only in Total Bed', 'Peaks detected in both', 'Peaks only in Pass Pileup']) # Set x-axis labels
    axs[1].set_ylabel('Confidence Values')
    axs[1].set_title('Vertical Box Plots for Confidence Levels of Peaks')


    plt.savefig('venn_fig.png')
