import os
import matplotlib.pyplot as plt
from sys import argv
import numpy as np
import pandas as pd


if __name__ == '__main__':
    summ_folder = argv[1]

    files_in_folder = []
    for item in os.listdir(summ_folder):
        item_path = os.path.join(summ_folder, item)
        if os.path.isfile(item_path):
            item_name = item_path.split("/")[-1]
            files_in_folder.append(item_name)

    # categories: files_in_folder
    perc_pass = []
    perc_fail = []
    raw_pass = []
    raw_fail = []

    percpass_flag = False
    percfail_flag = False
    rawpass_flag = False
    rawfail_flag = False
    for file in files_in_folder:
        with open(summ_folder + "/" + file, "r") as f:
            for line in f:
                if percpass_flag:
                    perc_pass.append(float(line.strip()))
                    percpass_flag = False
                elif percfail_flag:
                    perc_fail.append(float(line.strip()))
                    percfail_flag = False
                elif rawpass_flag:
                    raw_pass.append(float(line.strip()))
                    rawpass_flag = False
                elif rawfail_flag:
                    raw_fail.append(float(line.strip()))
                    rawfail_flag = False

                if "Pass Percentage:" in line:
                    percpass_flag = True
                elif "Fail Percentage:" in line:
                    percfail_flag = True
                elif "Pass Reads:" in line:
                    rawpass_flag = True
                elif "Fail Reads:" in line:
                    rawfail_flag = True
    print(str(perc_pass))
    print(str(perc_fail))

    data = {
        'Category': files_in_folder,
        'Pass_Perc': perc_pass,
        'Fail_Perc': perc_fail,
        'Pass_Reads': raw_pass,
        'Fail_Reads': raw_fail
    }
    df = pd.DataFrame(data)

    # Set the width of the bars
    x_pos = np.arange(len(files_in_folder))
    bar_width = 0.2

    fig, ax1 = plt.subplots(figsize=(10, 10))
    
    ax1.bar(x_pos - 0.3, perc_pass, bar_width, label='Pass Percentage', color='skyblue')
    ax1.bar(x_pos - 0.1, perc_fail, bar_width, label='Fail Percentage', color='lightcoral')
    ax1.set_ylabel('Percent')

    ax2 = ax1.twinx()
    ax2.bar(x_pos + 0.1, raw_pass, bar_width, label='Pass Reads', color='lightgreen') # Example with bars
    ax2.bar(x_pos + 0.3, raw_fail, bar_width, label='Fail Reads', color='gold') # Example with bars
    # Alternatively, for a line plot:
    # ax2.plot(df['Category'], df['Value2'], color='lightcoral', marker='o', label='Value 2')
    ax2.set_ylabel('Number of Reads')

    ax1.set_xlabel('Files')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(files_in_folder)

    plt.tight_layout() # Adjust layout to prevent labels from overlapping
    plt.title('Pass/Fail Summary for Each File')
    fig.legend(loc="upper left", bbox_to_anchor=(0.1, 0.9))

    plt.savefig(summ_folder + '/bar_plot.png')
