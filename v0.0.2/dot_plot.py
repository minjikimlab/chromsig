import matplotlib.pyplot as plt
from sys import argv
import os
import numpy as np

def test_line_chart(runtime_lst, title, labels, runtypes):
    """
    Use Matplotlib to visualize the difference between replicate scores and non-replicate scores with 1D line chart.
    Returns: updated ax. Usage:
        [calls ...]
        plt.tight_layout()
        plt.show()
    """
    label_size=5
    label_stagger=.026
    stagger_order=20
    s=2

    _, ax = plt.subplots(figsize=(12, 3.5 + len(labels)*1.5))
    if len(labels) > 1:
        plt.subplots_adjust(right=1.4)
    for i, exp in enumerate(labels):
        x_run = np.array(runtime_lst[i])
        #x_mem = np.array(memory_lst[i])
        rep_lbls = labels[i]
        run_name = runtypes[i]
        #nonrep_lbls = nonreps_labels[i]
        l = i+1

        for idx in range(len(rep_lbls)):
            print("File: " + rep_lbls[idx] + " Runtime: " + str(x_run[idx]))

        ax.scatter(x_run, [-l+.05] * len(x_run), color='tab:blue', alpha=0.7, label='Runtimes (s)' if i == 0 else "", s=s)
        #ax.scatter(x_mem, [-l+.1] * len(x_nonrep), color='tab:red', alpha=0.7, label='Memory Usage' if i == 0 else "", s=s)
        ax.text(1.05, -l, run_name, verticalalignment='center')

        #verbose = CONFIG.cfg['visualization']['verbose_labels']
        verbose = True
        if verbose:
            for j, val in enumerate(x_run):
                offset = 0.03 + label_stagger * (j % stagger_order)
                ax.text(val, -l - offset, str(rep_lbls[j]), fontsize=label_size, rotation=0, verticalalignment='bottom', color='darkblue')
            prev_offset = offset
            #for j, val in enumerate(x_nonrep):
            #    offset = prev_offset + label_stagger * (j % stagger_order)
            #    ax.text(val, -l - offset, str(nonrep_lbls[j]), fontsize=label_size, rotation=0, verticalalignment='bottom', color='darkred')

    ax.set_yticks([])
    if len(labels) > 1:
        ax.set_xlim(0, 3000)
    ax.set_xlabel("Runtimes (s)")
    ax.set_title(title)
    ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False)
    plt.tight_layout()

    plt.savefig("runtime_dot_plot.png")
    #VIZ_CONFIG['line_chart'].save()
    #plt.show()

def getdirs(lst):
    all_dirs = []
    for item in lst:
        if "_EnrichTest_" in item and "_100" in item:
            all_dirs.append(item)
    return all_dirs

def get_max_runtime(run_type, current_dir):
    all_files = os.listdir(run_type + "/" + current_dir)
    all_maxes = []
    for file in all_files:
        if "_logFile.txt" in file:
            with open(run_type + "/" + current_dir + "/" + file, "r") as logfile:
                for line in logfile:
                    if "Total: " in line:
                        all_maxes.append(float(line.split(" ")[1]))
    if (len(all_maxes) == 0):
        return 0
    else:
        return max(all_maxes)

if __name__ == '__main__':
    argv = argv[1].split(" ")
    atac_seq = argv[0]
    chip_seq = argv[1]
    cutandrun = argv[2]
    snatac_seq = argv[3]

    atac_files = os.listdir(atac_seq)
    chip_files = os.listdir(chip_seq)
    candr_files = os.listdir(cutandrun)
    snatac_files = os.listdir(snatac_seq)

    # Labels
    atac_dirs1 = getdirs(atac_files)
    chip_dirs1 = getdirs(chip_files)
    candr_dirs1 = getdirs(candr_files)
    snatac_dirs1 = getdirs(snatac_files)

    atac_times = []
    atac_dirs = []
    for item in atac_dirs1:
        max_time = get_max_runtime(atac_seq, item)
        if max_time == 0:
            continue
        else:
            atac_dirs.append(item)
            atac_times.append(max_time)

    chip_times = []
    chip_dirs = []
    for item in chip_dirs1:
        max_time = get_max_runtime(chip_seq, item)
        if max_time == 0:
            continue
        else:
            chip_dirs.append(item)
            chip_times.append(max_time)

    candr_times = []
    candr_dirs = []
    for item in candr_dirs1:
        max_time = get_max_runtime(cutandrun, item)
        if max_time == 0:
            continue
        else:
            candr_dirs.append(item)
            candr_times.append(max_time)
    
    snatac_times = []
    snatac_dirs = []
    for item in snatac_dirs1:
        max_time = get_max_runtime(snatac_seq, item)
        if max_time == 0:
            continue
        else:
            snatac_dirs.append(item)
            snatac_times.append(max_time)

    test_line_chart([atac_times, chip_times, candr_times, snatac_times], "Runtimes for Datasets", [atac_dirs, chip_dirs, candr_dirs, snatac_dirs], [atac_seq, chip_seq, cutandrun, snatac_seq])






