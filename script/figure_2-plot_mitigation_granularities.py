#!/usr/bin/env python3
import os
import sys
import json
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

def generate_plot_of_expected_wasted_capacity(interactive=True, savefig_filename=None):
    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))

    for mitigation_granularity in [1024, 512, 64, 32, 1]:
        x = []
        y = []
        print("[INFO] computing probabilites for a mitigation granularity of", mitigation_granularity)
        for rber in 10.0**np.arange(-8, 0, 0.01):
            expectation_wasted_bits = 0
            for n in range(1, mitigation_granularity + 1):
                P_n_errors = spsta.binom.pmf(n, mitigation_granularity, rber)
                n_wasted_bits = mitigation_granularity - n
                expectation_wasted_bits += P_n_errors * n_wasted_bits

            yval = expectation_wasted_bits / mitigation_granularity
            x.append(rber)
            y.append(yval)
        ax.plot(x, y, label=str(mitigation_granularity))

    ax.set_xscale('log')
    ax.grid(axis='y', ls='--')
    
    # ax.set_yscale('log')
    ax.legend(loc='upper left', fontsize=7, title=" Repair Granularity (Bits)", title_fontsize=7)
    ax.set_ylabel('Expected Wasted Storage\n(Ratio of Original Capacity)')
    ax.set_xlabel('Raw Bit Error Rate (RBER)')
    fig.tight_layout()
    fig.canvas.manager.set_window_title('Figure 2: Expected Wasted Storage')
    
    if interactive:
        plt.show()
    else:
        fname = 'fig_2_expected_wasted_storage.pdf' if savefig_filename == None else savefig_filename
        print("Saving figure:", fname)
        fig.savefig(fname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--interactive", help="draw interactive plot with matplotlib (else save plot to file)", action='store_true')
    parser.add_argument("-o", "--output-filename", help="custom filename of plot to save", type=str, default=None)
    args = parser.parse_args()

    generate_plot_of_expected_wasted_capacity(args.interactive, args.output_filename)