#!/usr/bin/env python3
import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcol
import matplotlib.collections as mc
import matplotlib.transforms as mtrans
import matplotlib.ticker as mtick
from itertools import chain
from collections import deque
import numpy as np
import scipy.stats as spsta
import math
import traceback

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
default_colors_hardcoded = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
blue = np.array([0.12156863, 0.46666667, 0.70588235])
orange = np.array([1.        , 0.49803922, 0.05490196])
green = np.array([0.17254902, 0.62745098, 0.17254902])
red =  np.array(mcol.to_rgb(default_colors[3]))

blue       = np.array([0x1f, 0x77, 0xb4], dtype=float) / 0xff
orange     = np.array([0xff, 0x7f, 0x0e], dtype=float) / 0xff
green      = np.array([0x2c, 0xa0, 0x2c], dtype=float) / 0xff
red        = np.array([0xd6, 0x27, 0x28], dtype=float) / 0xff
purple     = np.array([0x94, 0x67, 0xbd], dtype=float) / 0xff
brown      = np.array([0x8c, 0x56, 0x4b], dtype=float) / 0xff
pink       = np.array([0xe3, 0x77, 0xc2], dtype=float) / 0xff
grey       = np.array([0x7f, 0x7f, 0x7f], dtype=float) / 0xff
ugly_green = np.array([0xbc, 0xbd, 0x22], dtype=float) / 0xff
teal       = np.array([0x17, 0xbe, 0xcf], dtype=float) / 0xff

def generate_plot_of_postcorrection_probabilities(raw_data, interactive=True, savefig_basename=None):
    histograms = {}
    for fname in raw_data:
        for n_weak_bits in raw_data[fname]:
            if n_weak_bits not in histograms:
                histograms[n_weak_bits] = []
            for data_point in raw_data[fname][n_weak_bits]:
                histograms[n_weak_bits] += [i for i in data_point['post_probabilities'] if i != 0.0]
    
    for n_weak_bits in histograms:
        print(n_weak_bits, len(histograms[n_weak_bits]))
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 2.25))
    x = sorted(histograms.keys())
    y = [j for i, j in sorted(histograms.items(), key=lambda entry: entry[0])]
    y_raw = [[0.5] for _ in y]
    y_log = [[d for d in row] for row in y]

    labels = []
    v1 = ax.violinplot(y_raw, np.array(x) * 1.5 + 0, points=100, bw_method=0.2, showmedians=True, widths=0.5)
    v2 = ax.violinplot(y_log, np.array(x) * 1.5 + 0.5, points=100, bw_method=0.2, showmedians=True, widths=0.5)
    labels.append(mpatches.Patch(color=v1["bodies"][0].get_facecolor().flatten(), label="Pre-Correction"))
    labels.append(mpatches.Patch(color=v2["bodies"][0].get_facecolor().flatten(), label="Post-Correction"))

    # Make all the violin statistics marks red:
    for vp in [v1, v2]:
        vp['cmedians'].set_color('black')
        vp['cmedians'].set_linewidth(1)
        
        vp['cmins'].set_linewidth(2)
        vp['cmaxes'].set_linewidth(2)

    
    ax.set_xticks(np.array(x) * 1.5 + 0.25)
    ax.set_xticklabels(x)        
    ax.set_xlabel('Number of Pre-Correction Errors Per ECC Word')
    ax.set_ylabel('Per-Bit Probability\nof Post-Correction Error', position=(0, 0.4))
    ax.legend(handles=labels, fontsize=7, loc='upper left')
    fig.tight_layout()
    fig.canvas.manager.set_window_title("Figure 4: Distribution of error probabilites")
    if interactive:
        plt.show()
    else:
        fname = "fig_4_postcorrection_probabilities.pdf" 
        if savefig_basename != None:
            if os.path.isdir(savefig_basename):
                savefig_basename = os.path.join(os.path.normpath(savefig_basename), '')
            fname = savefig_basename + fname
        print("Saving figure:", fname)
        fig.savefig(fname)

def extract_bracketlist(spt, start_pos):
    assert spt[start_pos] == '['
    epos = start_pos + 1
    while epos < len(spt) and spt[epos] != ']':
        epos += 1
    assert epos < len(spt), "mismatched brace list!"
    return spt[start_pos + 1 : epos], epos

def parse_file(fname):
    data = {}
    with open(fname, 'r') as f:
        line_num = 1
        try:
            for line in f:
                if line.startswith("[INFO]"):
                    if line[7:].startswith("generating Hamming code"):
                        spt = line.strip().split(' ')
                        k = int(spt[6])
                        p = int(spt[8])
                        assert k == 64, "ERROR - only did this analysis for K=64"
                    elif line[7:].startswith("n_weak_bits"):
                        spt = line.strip().split(' ')
                        n_weak_bits = int(spt[2])
                        data[n_weak_bits] = []
                elif line.startswith("[DATA] weak_bits"):
                    end_pos = line.find('post_probabilities:')
                    assert end_pos != -1, "malformed line: " + line
                    weak_bit_list = eval(line[18:end_pos])
                    post_probabilities = eval(line[end_pos + 19:])
                    data[n_weak_bits].append({
                          'weak_bits' : weak_bit_list
                        , 'post_probabilities' : post_probabilities
                    })

                line_num += 1
        except:
            # print("[DEBUG] data:", data)
            print("[ERROR] failed parsing at fname:", fname, "line:", line_num)
            raise

    # print("[DEBUG]", data)
    print("[DEBUG] fname:", fname, "parsed")
    return data

def parse_files(filenames):
    print("[INFO] parsing", len(filenames), "input files")

    all_data = {}
    for fname in filenames:
        try:
            all_data[fname] = parse_file(fname)
        except Exception as e:
            print(e)
            print(traceback.format_exc())

    return all_data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--interactive", help="draw interactive plot with matplotlib (else save plot to file)", action='store_true')
    parser.add_argument("-o", "--output-filename-base", help="custom filename base of plots to save - will be suffixed with -figN.pdf", type=str, default=None)
    parser.add_argument("input_files_and_dirs", metavar="input-files-and-dirs", help="files (or directory containing files) to parse", nargs='+')
    args = parser.parse_args()

    # canonicalize input filenames
    cleaned_fnames = []
    for f_or_d in args.input_files_and_dirs:
        if os.path.isfile(f_or_d):
            cleaned_fnames.append(f_or_d)
        elif os.path.isdir(f_or_d):
            for root, dirs, files in os.walk(f_or_d):
                cleaned_fnames += [os.path.join(root, f) for f in files]
        else:
            print("[ERROR] invalid file/dir:", f_or_d)
            sys.exit(-1)

    if not cleaned_fnames:
        print("[WARN] no input files provided to parse!")
    else:
        all_data = parse_files(cleaned_fnames)
        generate_plot_of_postcorrection_probabilities(all_data, args.interactive, args.output_filename_base)

if __name__ == "__main__":
    main()    
