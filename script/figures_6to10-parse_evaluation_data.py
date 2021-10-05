#!/usr/bin/env python3
import os
import sys
import json
import argparse
import matplotlib as mpl
display_valid = "DISPLAY" in os.environ
if not display_valid:
    mpl.use('Agg')
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

json_file_list = {}
json_database = {}

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
black      = np.array([0x00, 0x00, 0x00], dtype=float) / 0xff

# deep sys.getsizeof() function
# https://code.activestate.com/recipes/577504/ 
def total_size(o, handlers={}, verbose=False):
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = sys.getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = sys.getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)

def compute_hamming_code_n_parity_bits(K):
    NmK = 0;
    while (1 << NmK) < (NmK + K + 1):
        NmK += 1
    assert NmK == math.ceil(math.log2(K + NmK + 1)), "MISMATCH - INCORRECT LOGIC"
    return NmK

# Replaecement for: axis.violinplot(ydata, xdata)
def discrete_violinplot(axis, ydata, xdata):
    vparts = axis.violinplot(ydata, xdata)
    for pc, yvals, xval in zip(vparts['bodies'], ydata, xdata):
        nbins = max(yvals) - min(yvals)
        if nbins:
            hist, bins = np.histogram(yvals, bins=nbins, density=True)
            widths = np.diff(bins)
            hist *= 0.2 / max(hist)
            # print(xval, set(yvals), hist)
            axis.hist(bins[:-1], bins, weights=hist, align='mid', bottom=xval, orientation='horizontal', histtype='bar', color=pc.get_facecolor().flatten())
            axis.hist(bins[:-1], bins, weights=-hist, align='mid', bottom=xval, orientation='horizontal', histtype='bar', color=pc.get_facecolor().flatten())
        # axis.hist(yvals, density=False, bottom=xval, orientation='horizontal', histtype='stepfilled', color=pc.get_facecolor().flatten())
        pc.set_visible(False)
    return vparts
       
def plot_figure6_direct_error_coverage_and_figure7_bootstrapping(raw_data, K, interactive, savefig_basename):

    # step 1: compute empirical distribution of rber 
    all_pbem = {"PBEM_100", "PBEM_75", "PBEM_50", "PBEM_25"}
    all_dp = {"RANDOM", "ALL_ONES", "COLSTRIPE", "BEEP"}
    all_n_prec_errors = {2, 3, 4, 5}
    
    N = K + compute_hamming_code_n_parity_bits(K)
    print("[INFO] analyzing codes with K:", K)

    n_words = 0
    direct_cov = {}
    direct_ideal = {}
    n_bootstrap_rounds = {}
    n_words_tested = {}
    for fname in raw_data:
        for json_cfg_fname in raw_data[fname]:
            if json_database[json_cfg_fname]['k'] != K:
                continue
            
            for n_prec_errors in raw_data[fname][json_cfg_fname]:
                if n_prec_errors not in direct_cov:
                    direct_cov[n_prec_errors] = {}
                    direct_ideal[n_prec_errors] = {}
                    n_bootstrap_rounds[n_prec_errors] = {}
                    n_words_tested[n_prec_errors] = {}

                for word in raw_data[fname][json_cfg_fname][n_prec_errors]:
                    prec_set = set(word['prec'])
                    posc_set = set(word['posc'])
                    direct_set = {i for i in posc_set if i in prec_set}
                    # indirect_set = {i for i in posc_set if i not in direct_set} # posc_set \ type1_set (i.e., set difference; basically everything else)

                    n_words += 1
                    for pbem in word['profiles']:
                        if pbem not in direct_cov[n_prec_errors]:
                            direct_cov[n_prec_errors][pbem] = {}
                            direct_ideal[n_prec_errors][pbem] = {}
                            n_bootstrap_rounds[n_prec_errors][pbem] = {}
                            n_words_tested[n_prec_errors][pbem] = {}
                        for dp in word['profiles'][pbem]:
                            if dp not in direct_cov[n_prec_errors][pbem]:
                                direct_cov[n_prec_errors][pbem][dp] = {}
                                direct_ideal[n_prec_errors][pbem][dp] = {}
                                n_bootstrap_rounds[n_prec_errors][pbem][dp] = {}
                                n_words_tested[n_prec_errors][pbem][dp] = {}
                            for profiler in word['profiles'][pbem][dp]:
                                if profiler not in direct_cov[n_prec_errors][pbem][dp]:
                                    direct_cov[n_prec_errors][pbem][dp][profiler] = [0 for _ in range(128)]
                                    direct_ideal[n_prec_errors][pbem][dp][profiler] = 0
                                    n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler] = []
                                    n_words_tested[n_prec_errors][pbem][dp][profiler] = 0
                                
                                # couting total number of bits so as to be able to normalize
                                n_words_tested[n_prec_errors][pbem][dp][profiler] += 1
                                
                                # measure the direct coverage
                                bit_set = set()
                                for bit_idx, round_idx in word['profiles'][pbem][dp][profiler]:
                                    if bit_idx in direct_set:
                                        bit_set.add(bit_idx) 
                                        for r_idx in range(round_idx, 128):
                                            direct_cov[n_prec_errors][pbem][dp][profiler][r_idx] += 1
                                direct_ideal[n_prec_errors][pbem][dp][profiler] += len(direct_set)

                                # we evaluate bootstrapping as the time spent in blind search
                                #   for each profiler. we calculate it as the number of rounds required for:
                                #       HARP-PREC: to identify the first prec error within the data bits
                                #       HARP-POSC: to identify the first direct error
                                #       NAIVE: to identify the first posc error
                                #       BEEP:  to identify the first miscorrection (i.e., the first prec error)
                                cur_set = word['profiles'][pbem][dp][profiler]
                                if profiler == 'beep_prec' or profiler == 'beep_posc': # both are more or less the same
                                    if len(posc_set) == 0:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(0)
                                    elif len(cur_set) == 0:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(128)
                                    else:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(min([round_idx for bit_idx, round_idx in cur_set]))
                                elif profiler == 'beep2_prec' or profiler == 'beep2_posc': # both are more or less the same
                                    pass # do nothing for beep2 since it is equivalent to HARP for purposes of direct errors
                                elif profiler == 'naive':
                                    if len(posc_set) == 0:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(0)
                                    elif len(cur_set) == 0:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(128)
                                    else:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(min([round_idx for bit_idx, round_idx in cur_set]))
                                elif profiler == 'harp_prec' or profiler == 'harp_posc': # both look for direct errors only
                                    if len(direct_set) == 0:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(0)
                                    elif len(cur_set) == 0:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(128)
                                    else:
                                        n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler].append(min([round_idx for bit_idx, round_idx in cur_set]))
                                else:  
                                    assert False, "unhandled profiler: " + str(profiler)
                                
                                # HARP-posc can miss some arbitrary postcorrection errors because they are due to UNKNOWN PARITY-CHECK errors
                                #       - these will be caught as indirect errors
                                # HARP-prec can have false positives because a miscorrection may occur in every case that the bit fails
                                if profiler == 'harp_prec' and len(word['profiles'][pbem][dp][profiler]) == 0 and any([i < K for i in prec_set]):
                                    print("[ERROR] HARP_PREC found no bits for:", n_prec_errors, "pbem:", pbem, "dp:", dp, "fname:", fname)
                                    print("    prec_set: ", prec_set)
                                    print("    posc_set: ", posc_set)
                                    print("    direct_set: ", direct_set)
                                    print("    raw profile: ", list(word['profiles'][pbem][dp][profiler]))
                                    assert False
                                if profiler == 'harp_prec' and direct_cov[n_prec_errors][pbem][dp][profiler][-1] != direct_ideal[n_prec_errors][pbem][dp][profiler]:
                                    print("[ERROR] HARP_PREC did not achieve full direct coverage with config n_prec_errors:", n_prec_errors, "pbem:", pbem, "dp:", dp, "fname:", fname)
                                    print("    bit_set: ", bit_set)
                                    print("    direct_set: ", direct_set)
                                    print("    prec_set: ", prec_set)
                                    print("    posc_set: ", posc_set)
                                    print("    raw profile: ", list(word['profiles'][pbem][dp][profiler]))
                                    print("    direct_cov: ", direct_cov[n_prec_errors][pbem][dp][profiler][-1])
                                    print("    direct_cov_ideal: ", direct_ideal[n_prec_errors][pbem][dp][profiler])
                                    assert False
                # print(fname, json_cfg_fname, n_prec_errors, direct_cov[n_prec_errors])
    print("[INFO] as if we are profiling n_words:", n_words, "n_bytes:", n_words * K / 8)

    ########################################################################################################################
    # x-axis is number of prec errors
    ########################################################################################################################
    direct_profilers = ['naive', 'beep_posc', 'harp_prec']
    direct_profiler_labels = {'naive' : "Naive", 'beep_posc' : "BEEP", 'harp_prec' : "HARP-U"}
    profiler_colors = {'naive' : orange, 'beep_posc' : blue, 'harp_prec' : green}
    pbem_to_ax_index = {"PBEM_100" : 3, "PBEM_75" : 2, "PBEM_50" : 1, "PBEM_25" : 0}
    ax_index_to_label = {0 : '25%', 1 : '50%', 2 : '75%', 3 : '100%'}
    n_prec_to_marker = {2 : 'None', 3 : "+", 4 : 'x', 5 : '.'}
    if True:
        # for g_pbem in all_pbem:
        # for g_dp in all_dp:
        # for g_n_prec_errors in all_n_prec_errors:
        fig, ax = plt.subplots(1, 4, figsize=(5, 2.25), sharex=True, sharey=True)
        x = range(1, 129)

        dp_to_marker = {"RANDOM" : 'None', "ALL_ONES" : 'x', "COLSTRIPE" : '1', "COMBINED" : '|', "BEEP" : 'None'}
        markersize = 6
        ncharged_to_shade = {2 : 0.33, 3 : 0.55, 4 : 0.77, 5 : 1.0}
        labels = []
        for profiler in direct_profilers:
            for n_prec_errors in direct_cov:
                # if n_prec_errors != g_n_prec_errors:
                #     continue
                for pbem in direct_cov[n_prec_errors]:
                    # if pbem != g_pbem:
                    #     continue
                    for dp in direct_cov[n_prec_errors][pbem]:
                        # if dp not in ["BEEP", "COMBINED", "RANDOM"]:
                        if dp not in ["BEEP", "RANDOM"]:
                            continue
                        # if dp != g_dp:
                        #     continue
                        if profiler not in direct_cov[n_prec_errors][pbem][dp]:
                            continue

                        # direct cov is the same for prec/posc profilers
                        ideal = direct_ideal[n_prec_errors][pbem][dp][profiler]
                        y = [r / ideal for r in direct_cov[n_prec_errors][pbem][dp][profiler]]
                    
                        if profiler in profiler_colors:
                            ax[pbem_to_ax_index[pbem]].plot(x, y, color=profiler_colors[profiler], lw=1) #, marker=dp_to_marker[dp], ms=markersize)
                            po2_lut = [1, 2, 4, 8, 16, 32, 64, 127]
                            ax[pbem_to_ax_index[pbem]].plot([x[i] for i in po2_lut], [y[i] for i in po2_lut], ls='None', color=profiler_colors[profiler], marker=n_prec_to_marker[n_prec_errors], ms=markersize) #, marker=dp_to_marker[dp], ms=markersize)
                        else:    
                            assert False, "missing profiler"

        labels = []
        for profiler in ['harp_prec', 'naive', 'beep_posc']:
            labels.append(mpatches.Patch(facecolor=profiler_colors[profiler], label=direct_profiler_labels[profiler]))
        labels_shades = []
        # for ncharged in sorted(ncharged_to_shade):
        #     labels_shades.append(mpatches.Patch(facecolor=np.array((0.75, 0.75, 0.75)) * ncharged_to_shade[ncharged], label=str(ncharged)))
        for n_prec in sorted(n_prec_to_marker):
            labels_shades.append(mlines.Line2D([], [], color='k', marker=n_prec_to_marker[n_prec], label=str(n_prec)))
        labels_dps = []
        for dp in ["RANDOM", "ALL_ONES", "COLSTRIPE", "COMBINED"]:
            if dp == 'BEEP':
                continue
            labels_dps.append(mlines.Line2D([], [], color='k', marker=dp_to_marker[dp], ms=markersize, label=dp))


        for i in range(len(ax)):
            ax[i].set_xlabel(ax_index_to_label[i])
            ax[i].xaxis.set_label_position('top') 
            
            # ax[i].set_xscale('log', basex=2) # deprecated since Matploblib 3.3
            try:
                ax[i].set_xscale('log', base=2)
            except:
                ax[i].set_xscale('log', basex=2)

            # ax[i].set_xticks([1] + list(range(32, 128, 32)))
            # xlabels = [1] + list(range(32, 128, 32)) # LINEAR
            # xlabels = [1, 4, 16, 128] # LOG
            # ax[i].set_xticks(xlabels)
            # ax[i].set_xticklabels(xlabels)

            # x_major = mtick.LogLocator(base = 10.0, numticks = 5)
            # ax.xaxis.set_major_locator(x_major)
            ax[i].tick_params(axis='y', which='minor', left=False, right=False)
            ax[i].tick_params(axis='x', which='minor', top=False, bottom=True)
            ax[i].xaxis.set_minor_locator(mtick.AutoMinorLocator())
            ax[i].minorticks_on()
            ax[i].set_xticks(np.array([np.array([1.25, 1.5, 1.75]) * n for n in 2**np.arange(0, 7)]).flatten(), minor=True)
            ax[i].xaxis.set_minor_formatter(mtick.NullFormatter())

            ax[i].set_xticks([1, 2, 4, 8, 16, 32, 64, 128])
            ax[i].set_xticklabels([1, None, 4, None, 16, None, 64, None])
            ax[i].set_xlim([0.75, 192])
            
            # ax[i].tick_params(which='major', width=1.00)
            # ax[i].tick_params(which='major', length=5)
            # ax[i].tick_params(which='minor', width=0.75)
            # ax[i].tick_params(which='minor', length=2.5)

        fig.text(0.58, 0.92, 'Per-Bit Probability of Pre-Correction Error', ha='center')    
        fig.text(0.56, 0.04, 'Number of Profiling Rounds', ha='center')

        ax[0].set_ylabel('Direct Error Coverage')
        ax[0].set_ylim([-0.05, 1.05])
        ax[3].legend(handles=list(labels), labelspacing=0.1, fontsize=6, ncol=1, loc='lower right')
        ax[2].legend(handles=list(reversed(labels_shades)), labelspacing=0.1, fontsize=6, ncol=2, loc='lower right', title="Pre-Correction Errors", title_fontsize=6)
        # ax[1].legend(handles=list(labels_dps), labelspacing=0.1, fontsize=6, ncol=1, loc='lower right')
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.1, bottom=0.22, top=0.82)
        fig.canvas.manager.set_window_title("Figure 6: direct error coverage vs. round - K: " + str(K))
        if not interactive:
            fname = 'fig_6_direct_error_coverage_vs_round.pdf' 
            if savefig_basename != None:
                fname = savefig_basename + fname
            print("[INFO] Saving figure:", fname)
            fig.savefig(fname)

    ########################################################################################################################
    # bootstrapping analysis
    ########################################################################################################################
    prof_to_ax_index = {'naive' : 0, 'beep_posc' : 1, 'harp_prec' : 2}
    ax_index_to_label = {0 : "Naive", 1 : "BEEP", 2 : "HARP-U"}
    pbem_to_color = {'PBEM_25' : default_colors[0], 'PBEM_50' : default_colors[1], 'PBEM_75' : default_colors[2], 'PBEM_100' : default_colors[3]}
    labels_pbem = []
    xval = [0 for _ in pbem_to_ax_index]
    xlabels = [[] for _ in pbem_to_ax_index]
    xticks = [[] for _ in pbem_to_ax_index]
    if True:
        fig, axs = plt.subplots(1, 3, figsize=(5, 2.25), sharex=True, sharey=True)
        for profiler in prof_to_ax_index:
            ax_idx = prof_to_ax_index[profiler]
            for n_prec_errors in [2, 3, 4, 5]:
                for pbem, pbem_label in reversed([("PBEM_100", "100%"), ("PBEM_75", "75%"), ("PBEM_50", "50%"), ("PBEM_25", "25%")]):
                    for dp in n_bootstrap_rounds[n_prec_errors][pbem]:
                        if dp not in ["BEEP", "RANDOM"]: continue
                        if n_prec_errors in n_bootstrap_rounds \
                            and pbem in n_bootstrap_rounds[n_prec_errors] \
                            and dp in n_bootstrap_rounds[n_prec_errors][pbem] \
                            and profiler in n_bootstrap_rounds[n_prec_errors][pbem][dp]:
                            n_bootstrap_rounds_dist = n_bootstrap_rounds[n_prec_errors][pbem][dp][profiler]
                            
                            parts = axs[ax_idx].violinplot(n_bootstrap_rounds_dist, [xval[ax_idx]], widths=0.5, showmeans=False, showmedians=True, showextrema=True)
                            for partname in ('cbars',):
                                vp = parts[partname]
                                vp.set_edgecolor(pbem_to_color[pbem])
                                vp.set_linewidth(0.5)
                            linepairs = []
                            for partname in ('cmins','cmaxes', 'cmedians'):
                                vp = parts[partname]
                                vp.set_edgecolor(pbem_to_color[pbem])
                                for path in vp.get_paths():
                                    start, end = path.vertices
                                    xmed = (start[0] + end[0]) / 2.0
                                    start[0] -= (xmed - start[0])
                                    end[0] += (end[0] - xmed)
                                    linepairs.append((start, end))
                                vp.set_visible(False)
                            lc = mc.LineCollection(linepairs, colors=pbem_to_color[pbem], linewidths=1.5)
                            axs[ax_idx].add_collection(lc)
                            for pc in parts['bodies']:
                                pc.set_facecolor(pbem_to_color[pbem])
                                pc.set_edgecolor('black')
                    xticks[ax_idx].append(xval[ax_idx])                    
                    xlabels[ax_idx].append(str(n_prec_errors))
                    xval[ax_idx] += 0.5
                xval[ax_idx] += 0.5
        
        for pbem, pbem_label in reversed([("PBEM_100", "100%"), ("PBEM_75", "75%"), ("PBEM_50", "50%"), ("PBEM_25", "25%")]):
            labels_pbem.append(mpatches.Patch(facecolor=pbem_to_color[pbem], label=pbem_label))

        for i in range(len(axs)):
            xticks[i] = [np.mean(xticks[i][j:j+4]) for j in range(0, len(xticks[i]), 4)]
            xlabels[i] = xlabels[i][::4]

            axs[i].set_xlabel(ax_index_to_label[i])
            axs[i].xaxis.set_label_position('top') 
            # axs[i].set_xticks([1] + list(range(32, 128, 32)))
            # xlabels = [1] + list(range(32, 128, 32)) # LINEAR
            # xlabels = [1, 2, 4, 8, 20, 64] # LOG
            axs[i].set_xticks(xticks[i])
            axs[i].set_xticklabels(xlabels[i])
            # axs[i].set_xticklabels(xlabels[i], rotation=90, fontsize=6)
            axs[i].set_yticks([0, 32, 64, 96, 128])
        fig.text(0.57, 0.04, 'Number of Pre-Correction Errors Injected', ha='center')

        axs[0].set_xlim([-0.5, xval[0] - 0.5])
        axs[0].set_ylabel('Number of Rounds\nSpent Bootstrapping')
        # axs[3].legend(handles=list(labels), labelspacing=0.1, fontsize=6, ncol=1, loc='lower right')
        axs[2].legend(handles=list(labels_pbem), labelspacing=0.1, fontsize=6, ncol=1, loc='upper right', title="Pre-Correction Error\n  Per-Bit Probability", title_fontsize=6)
        # axs[1].legend(handles=list(labels_dps), labelspacing=0.1, fontsize=6, ncol=1, loc='lower right')
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.1, bottom=0.22)
        fig.canvas.manager.set_window_title("Figure 7: Bootstrapping vs n_prec_errors - K: " + str(K))
        if not interactive:
            fname = 'fig_7_bootstrapping_vs_nprecerrs.pdf' 
            if savefig_basename != None:
                fname = savefig_basename + fname
            print("[INFO] Saving figure:", fname)
            fig.savefig(fname)
        
def plot_figure8_indirect_error_coverage(raw_data, K, interactive, savefig_basename):
    prec_profilers = {'beep_prec', 'harp_prec', 'beep2_prec'}
    posc_profilers = {'naive', 'beep_posc', 'harp_posc', 'beep2_posc'}
    n_words = 0
    missed_indirect_errors = {}
    n_words_tested = {}
    for fname in raw_data:
        for json_cfg_fname in raw_data[fname]:
            if json_database[json_cfg_fname]['k'] != K:
                continue
            N = K + compute_hamming_code_n_parity_bits(K)
            
            for n_prec_errors in raw_data[fname][json_cfg_fname]:
                if n_prec_errors not in missed_indirect_errors:
                    missed_indirect_errors[n_prec_errors] = {}
                    n_words_tested[n_prec_errors] = {}

                for word in raw_data[fname][json_cfg_fname][n_prec_errors]:
                    prec_set = set(word['prec'])
                    posc_set = set(word['posc'])
                    direct_set = {i for i in posc_set if i in prec_set}
                    indirect_set = {i for i in posc_set if i not in prec_set} # posc_set \ type1_set (i.e., set difference; basically everything else)

                    n_words += 1
                    for pbem in word['profiles']:
                        if pbem not in missed_indirect_errors[n_prec_errors]:
                            missed_indirect_errors[n_prec_errors][pbem] = {}
                            n_words_tested[n_prec_errors][pbem] = {}
                        for dp in word['profiles'][pbem]:
                            if dp not in missed_indirect_errors[n_prec_errors][pbem]:
                                missed_indirect_errors[n_prec_errors][pbem][dp] = {}
                                n_words_tested[n_prec_errors][pbem][dp] = {}
                            for profiler in word['profiles'][pbem][dp]:
                                if profiler not in missed_indirect_errors[n_prec_errors][pbem][dp]:
                                    missed_indirect_errors[n_prec_errors][pbem][dp][profiler] = [0 for _ in range(128)]
                                    n_words_tested[n_prec_errors][pbem][dp][profiler] = 0
                                
                                # couting total number of bits so as to be able to normalize
                                n_words_tested[n_prec_errors][pbem][dp][profiler] += 1
                                
                                # counting the number of indirect bits remaining
                                for r in range(128):
                                    if profiler in prec_profilers:
                                        prec_found_set = {bit_idx for bit_idx, round_idx in word['profiles'][pbem][dp][profiler] if round_idx <= r}
                                        posc_found_set = prec_found_set # we simply assume that pre == post (only really counts for HARP+, introduces false positives for HARP)
                                        missed_set = {bit_idx for bit_idx in indirect_set if bit_idx not in posc_found_set}
                                        missed_indirect_errors[n_prec_errors][pbem][dp][profiler][r] += len(missed_set)
                                    elif profiler in posc_profilers:
                                        posc_found_set = {bit_idx for bit_idx, round_idx in word['profiles'][pbem][dp][profiler] if round_idx <= r}
                                        missed_set = {bit_idx for bit_idx in indirect_set if bit_idx not in posc_found_set}
                                        missed_indirect_errors[n_prec_errors][pbem][dp][profiler][r] += len(missed_set)
                                    else:
                                        assert False, "unsupported/unknown profiler: " + str(profiler)
                                # print(indirect_set)
                                # print(missed_indirect_errors[n_prec_errors][pbem][dp][profiler][0])
    
    ########################################################################################################################
    # indirect remaining analysis
    ########################################################################################################################
    indirect_profilers = ['naive', 'beep_posc', 'beep2_posc', 'harp_prec', 'harp_posc']
    indirect_profiler_labels = {'naive' : "Naive", 'beep_posc' : "BEEP", 'beep2_posc' : "HARP-A+BEEP", 'harp_prec' : "HARP-U", 'harp_posc' : "HARP-A"}
    profiler_colors = {'naive' : orange, 'beep_posc' : blue, 'beep2_posc' : purple, 'harp_prec' : green, 'harp_posc' : red}
    nprec_to_ax_index = {2 : 0, 3 : 1, 4 : 2, 5 : 3}
    ax_index_to_label = {0 : "2", 1 : "3", 2 : "4", 3 : "5"}
    # pbem_to_shade = {"PBEM_100" : 0.50, "PBEM_75" : 0.67, "PBEM_50" : 0.83, "PBEM_25" : 1.0}
    pbem_to_marker = {"PBEM_100" : 'None', "PBEM_75" : "+", "PBEM_50" : 'x', "PBEM_25" : '.'}
    markersize={"PBEM_100" : 4, "PBEM_75" : 5, "PBEM_50" : 4, "PBEM_25" : 6}
    marker_indices = [1, 2, 4, 8, 16, 32, 64, 127]
    if True:
        # for g_pbem in all_pbem:
        # for g_dp in all_dp:
        # for g_n_prec_errors in all_n_prec_errors:
        fig, axs = plt.subplots(1, 4, figsize=(5, 2.25), sharex=True, sharey=True)
        x = range(1, 129)

        # dp_to_marker = {"RANDOM" : 'None', "ALL_ONES" : 'x', "COLSTRIPE" : '1', "COMBINED" : '|', "BEEP" : 'None'}
        labels = []
        for profiler in indirect_profilers:
            for n_prec_errors in missed_indirect_errors:
                ax_idx = nprec_to_ax_index[n_prec_errors]
                for pbem in missed_indirect_errors[n_prec_errors]:
                    # if pbem != 'PBEM_50': continue
                    for dp in missed_indirect_errors[n_prec_errors][pbem]:
                        # if dp not in ["BEEP", "COMBINED", "RANDOM"]:
                        if dp not in ["BEEP", "BEEP2", "RANDOM"]:
                            continue
                        if profiler not in missed_indirect_errors[n_prec_errors][pbem][dp]:
                            continue

                        # direct cov is the same for prec/posc profilers
                        y = np.array(missed_indirect_errors[n_prec_errors][pbem][dp][profiler], dtype=float) / n_words_tested[n_prec_errors][pbem][dp][profiler]
                    
                        if profiler in profiler_colors:
                            axs[ax_idx].plot(x, y, color=profiler_colors[profiler], lw=1) #, marker=dp_to_marker[dp], ms=markersize)
                            axs[ax_idx].plot([x[i] for i in marker_indices], [y[i] for i in marker_indices], ls='None', color=profiler_colors[profiler]
                                , marker=pbem_to_marker[pbem], ms=markersize[pbem], mew=0.5)
                            # axs[ax_idx].plot(x, y, color=profiler_colors[profiler] * pbem_to_shade[pbem]) #, marker=dp_to_marker[dp], ms=markersize)
                            # axs[ax_idx].plot(x[::32], y[::32], ls='None', color=profiler_colors[profiler] * pbem_to_shade[pbem], marker=dp_to_marker[dp], ms=markersize) #, marker=dp_to_marker[dp], ms=markersize)
                        else:    
                            assert False, "missing profiler"

        labels = []
        for profiler in ['harp_posc', 'harp_prec', 'naive', 'beep_posc', 'beep2_posc']:
            labels.append(mpatches.Patch(facecolor=profiler_colors[profiler], label=indirect_profiler_labels[profiler]))
        labels_pbem = []
        for pbem, pbem_label in reversed([("PBEM_100", "100%"), ("PBEM_75", "75%"), ("PBEM_50", "50%"), ("PBEM_25", "25%")]):
            labels_pbem.append(mlines.Line2D([], [], color='black', linestyle='-', marker=pbem_to_marker[pbem], ms=markersize[pbem], label=pbem_label))
        # labels_dps = []
        # for dp in ["RANDOM", "ALL_ONES", "COLSTRIPE", "COMBINED"]:
        #     if dp == 'BEEP':
        #         continue
        #     labels_dps.append(mlines.Line2D([], [], color='k', marker=dp_to_marker[dp], ms=markersize[pbem], label=dp))

        axs[0].set_xscale('log')
        for i in range(len(axs)):
            axs[i].set_xlabel(ax_index_to_label[i])
            axs[i].xaxis.set_label_position('top') 

            axs[i].tick_params(axis='y', which='minor', left=False, right=False)
            axs[i].tick_params(axis='x', which='minor', top=False, bottom=True)
            axs[i].xaxis.set_minor_locator(mtick.AutoMinorLocator())
            axs[i].minorticks_on()
            axs[i].set_xticks(np.array([np.array([1.25, 1.5, 1.75]) * n for n in 2**np.arange(0, 7)]).flatten(), minor=True)
            axs[i].xaxis.set_minor_formatter(mtick.NullFormatter())

            axs[i].set_xticks([1, 2, 4, 8, 16, 32, 64, 128])
            axs[i].set_xticklabels([1, None, 4, None, 16, None, 64, None])
            axs[i].set_xlim([0.75, 192])
        
        fig.text(0.58, 0.92, 'Number of Pre-Correction Errors per ECC Word', ha='center')
        fig.text(0.58, 0.04, 'Number of Profiling Rounds', ha='center')

        axs[0].set_ylabel('Missed Indirect Errors\nper ECC Word')
        axs[0].legend(handles=list(labels_pbem), labelspacing=0.1, fontsize=6, ncol=1, loc='upper left', title="Per-Bit Probability\n of Pre-Correction\n         Error\n", title_fontsize=6)
        axs[1].legend(handles=list(labels), handletextpad=0.2, labelspacing=0.1, fontsize=6, ncol=1, loc='upper left')
        # axs[1].legend(handles=list(labels_dps), labelspacing=0.1, fontsize=6, ncol=1, loc='lower right')
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.1, bottom=0.22, top=0.82)
        fig.canvas.manager.set_window_title("Figure 8: Missed indirect errors per ECC word - K: " + str(K))
        if not interactive:
            fname = 'fig_8_missed_indirect_errors_per_ecc_word.pdf' 
            if savefig_basename != None:
                fname = savefig_basename + fname
            print("[INFO] Saving figure:", fname)
            fig.savefig(fname)

def plot_figure9ab_hisograms_of_simultaneous_postcorrection_errors(raw_data, K, interactive, savefig_basename):
    histograms = {}
    n_words_tested = {}
    for fname in raw_data:
        for json_cfg_fname in raw_data[fname]:
            if json_database[json_cfg_fname]['k'] != K:
                continue
            N = K + compute_hamming_code_n_parity_bits(K)
            
            for n_prec_errors in raw_data[fname][json_cfg_fname]:
                if n_prec_errors not in histograms:
                    histograms[n_prec_errors] = {}
                    n_words_tested[n_prec_errors] = {}

                for word in raw_data[fname][json_cfg_fname][n_prec_errors]:
                    prec_set = set(word['prec'])
                    posc_set = set(word['posc'])
                    direct_set = {i for i in posc_set if i in prec_set}
                    indirect_set = {i for i in posc_set if i not in direct_set} # posc_set \ type1_set (i.e., set difference; basically everything else)

                    for pbem in word['histograms']:
                        if pbem not in histograms[n_prec_errors]:
                            histograms[n_prec_errors][pbem] = {}
                            n_words_tested[n_prec_errors][pbem] = {}
                        for dp in word['histograms'][pbem]:
                            if dp not in histograms[n_prec_errors][pbem]:
                                histograms[n_prec_errors][pbem][dp] = {}
                                n_words_tested[n_prec_errors][pbem][dp] = {}
                            for profiler in word['histograms'][pbem][dp]:
                                if profiler not in histograms[n_prec_errors][pbem][dp]:
                                    histograms[n_prec_errors][pbem][dp][profiler] = [[0 for j in range(K)] for i in range(129 + 127)]
                                    n_words_tested[n_prec_errors][pbem][dp][profiler] = 0
                                
                                # couting total number of bits so as to be able to normalize
                                n_words_tested[n_prec_errors][pbem][dp][profiler] += 1
                                
                                # counting the number of indirect bits remaining
                                hist_as_dict = {round_idx : n_simultaneous_errors for round_idx, n_simultaneous_errors in word['histograms'][pbem][dp][profiler]}
                                if len(hist_as_dict) == 0:
                                    hist_as_dict[0] = 0
                                hist_as_dict[129] = max(hist_as_dict.items(), key=lambda x: x[0])[1]
                                round_list = sorted([round_idx for round_idx, n_simultaneous_errors in hist_as_dict.items()])
                                for r_start, r_end in [(round_list[i], round_list[i + 1]) for i in range(len(round_list) - 1)]:
                                    for r in range(r_start, r_end):
                                        n_simultaneous_errors = hist_as_dict[r_start]
                                        histograms[n_prec_errors][pbem][dp][profiler][r][n_simultaneous_errors] += 1
    
    ########################################################################################################################
    # indirect remaining analysis
    ########################################################################################################################
    indirect_profilers = ['naive', 'beep_posc', 'harp_prec', 'harp_posc']
    indirect_profiler_labels = {'naive' : "Naive", 'beep_posc' : "BEEP", 'harp_prec' : "HARP-U", 'harp_posc' : "HARP-A"}
    profiler_colors = {'naive' : orange, 'beep_posc' : blue, 'harp_prec' : green, 'harp_posc' : red}
    nprec_to_axc_index = {2 : 0, 3 : 1, 4 : 2, 5 : 3}
    pbem_to_axr_index = {"PBEM_100" : 0, "PBEM_75" : 1, "PBEM_50" : 2, "PBEM_25" : 3}
    ax_index_to_label = {0 : "2", 1 : "3", 2 : "4", 3 : "5"}
    # pbem_to_shade = {"PBEM_100" : 0.50, "PBEM_75" : 0.67, "PBEM_50" : 0.83, "PBEM_25" : 1.0}
    pbem_to_marker = {"PBEM_100" : 'None', "PBEM_75" : "+", "PBEM_50" : 'x', "PBEM_25" : '.'}
    markersize={"PBEM_100" : 4, "PBEM_75" : 5, "PBEM_50" : 4, "PBEM_25" : 6}
    marker_indices = [1, 2, 4, 8, 16, 32]
    if True:
        # for g_pbem in all_pbem:
        # for g_dp in all_dp:
        # for g_n_prec_errors in all_n_prec_errors:
        # fig, axs = plt.subplots(4, 4, figsize=(15, 3), sharex=True, sharey=True)
        fig, axs = plt.subplots(4, 4, figsize=(6, 3), sharex=True, sharey=True)
        x = np.array(range(0, K), dtype=float)
        bar_width = 1.0 / 6
        bar_offset = {'naive' : -bar_width * 2, 'beep_posc' : -bar_width, 'harp_prec' : 0, 'harp_posc' : bar_width}

        # dp_to_marker = {"RANDOM" : 'None', "ALL_ONES" : 'x', "COLSTRIPE" : '1', "COMBINED" : '|', "BEEP" : 'None'}
        labels = []
        for profiler in indirect_profilers:
            for n_prec_errors in histograms:
                axc_idx = nprec_to_axc_index[n_prec_errors]
                for pbem in histograms[n_prec_errors]:
                    # if pbem != 'PBEM_50': continue
                    axr_idx = pbem_to_axr_index[pbem]
                    for dp in histograms[n_prec_errors][pbem]:
                        # if dp not in ["BEEP", "COMBINED", "RANDOM"]:
                        if dp not in ["BEEP", "RANDOM"]:
                            continue
                        if profiler not in histograms[n_prec_errors][pbem][dp]:
                            continue

                        # direct cov is the same for prec/posc profilers
                        round = 127
                        # y = np.array(histograms[n_prec_errors][pbem][dp][profiler][round], dtype=float) / n_words_tested[n_prec_errors][pbem][dp][profiler]
                        this_hist = histograms[n_prec_errors][pbem][dp][profiler][round]
                        y = np.array(this_hist, dtype=float) / sum(this_hist)
                        
                        if profiler in profiler_colors:
                            axs[axr_idx][axc_idx].bar(x + bar_offset[profiler], y, bar_width, color=profiler_colors[profiler], align='edge', edgecolor='black', linewidth=0.5)
                        else:    
                            assert False, "missing profiler"

        labels = []
        for profiler in ['naive', 'beep_posc', 'harp_prec', 'harp_posc']:
            labels.append(mpatches.Patch(facecolor=profiler_colors[profiler], label=indirect_profiler_labels[profiler]))
        labels_pbem = []
        for pbem, pbem_label in reversed([("PBEM_100", "100%"), ("PBEM_75", "75%"), ("PBEM_50", "50%"), ("PBEM_25", "25%")]):
            labels_pbem.append(mlines.Line2D([], [], color='black', linestyle='-', marker=pbem_to_marker[pbem], ms=markersize[pbem], label=pbem_label))
        # labels_dps = []
        # for dp in ["RANDOM", "ALL_ONES", "COLSTRIPE", "COMBINED"]:
        #     if dp == 'BEEP':
        #         continue
        #     labels_dps.append(mlines.Line2D([], [], color='k', marker=dp_to_marker[dp], ms=markersize[pbem], label=dp))

        # axs[0].set_xscale('log')

        # for pbem in histograms[n_prec_errors]:
        #     axr_idx = pbem_to_axr_index[pbem]
        #     axs[axr_idx][3].set_title(str(pbem))


        for i in range(len(axs)):
            axs[0][i].set_xlabel(ax_index_to_label[i])
            axs[0][i].xaxis.set_label_position('top') 
            # axs[i][0].set_xticks([1] + list(range(32, 128, 32)))
            # xlabels = [1] + list(range(32, 128, 32)) # LINEAR
            # xlabels = [1, 2, 4, 8, 20, 64] # LOG
            xlabels = [0, 1, 2, 3, 4, 5, 6]
            xtickpos = np.array([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
            axs[3][i].set_xticks(xlabels)
            axs[3][i].set_xticklabels(xlabels)
            axs[3][i].xaxis.set_minor_locator(mtick.FixedLocator(xtickpos))
            axs[3][i].tick_params(axis='x', which='major', length=0)
            axs[3][i].tick_params(axis='x', which='minor', length=12)
            axs[3][i].set_xlim([-0.5, 6.5])

            axs[3][i].set_ylim([0, 1.3])
            axs[3][i].set_yticks([0, 0.5, 1])

        fig.text(0.50, 0.93, 'Number of Pre-Correction Errors per ECC Word', ha='center')
        fig.text(0.50, 0.04, 'Maximum Number of Simultaneous Post-Correction Errors Possible', ha='center')
        fig.text(0.03, 0.22, 'Fraction of ECC Words', ha='center', rotation=90)
        fig.text(0.96, 0.28, 'Pre-Correction Error\nPer-Bit Probability', ha='center', rotation=270)
        fig.text(0.92, 0.72, '100%', ha='center', rotation=270)
        fig.text(0.92, 0.55, '75%', ha='center', rotation=270)
        fig.text(0.92, 0.38, '50%', ha='center', rotation=270)
        fig.text(0.92, 0.20, '25%', ha='center', rotation=270)

        # axs[0][0].set_yscale('log')
        # axs[0][0].set_ylabel('Number of ECC Words')
        # axs[3][0].legend(handles=list(labels_pbem), labelspacing=0.1, fontsize=6, ncol=1, loc='upper right', title="Pre-Correction Error\n  Per-Bit Probability", title_fontsize=6)
        axs[0][0].legend(handles=list(labels), labelspacing=0.01, fontsize=6, ncol=1, loc='upper right', frameon=False, borderpad=0.02)

        # axs[1][0].legend(handles=list(labels_dps), labelspacing=0.1, fontsize=6, ncol=1, loc='lower right')
        fig.tight_layout()
        # fig.subplots_adjust(wspace=0.1, hspace=0, bottom=0.16, top=0.86, left=0.06, right=0.94)
        fig.subplots_adjust(wspace=0.1, hspace=0, bottom=0.16, top=0.86, left=0.10, right=0.90)
        fig.canvas.manager.set_window_title("Figure 9a: incomplete coverage after profiling - K: " + str(K) + " subfigure(a)")
        if not interactive:
            fname = 'fig_9a_incomplete_coverage_after_profiling.pdf' 
            if savefig_basename != None:
                fname = savefig_basename + fname
            print("[INFO] Saving figure:", fname)
            fig.savefig(fname)
    
    # Q: how many rounds until the profiler achieves coverage to N-max remaining?  
    if True:
        all_pbem = {"PBEM_100", "PBEM_75", "PBEM_50", "PBEM_25"}
        all_dp = {"RANDOM", "ALL_ONES", "COLSTRIPE", "BEEP"}
        all_n_prec_errors = {2, 3, 4, 5}

        for percentile in [0.9, 0.99, 1.00]:
            # fig, axs = plt.subplots(4, 4, figsize=(15, 3), sharex=True, sharey=True)
            fig, axs = plt.subplots(4, 4, figsize=(6, 3), sharex=True, sharey=True)

            labels = []
            for profiler in indirect_profilers:
                # axc_idx = indirect_profiler_to_axc_idx[profiler]
                for pbem in all_pbem:
                    # if pbem != 'PBEM_50': continue
                    axr_idx = pbem_to_axr_index[pbem]
                    for dp in all_dp:
                        # axc_idx = dp_to_axc_idx[dp]
                        # if dp not in ["BEEP", "COMBINED", "RANDOM"]:
                        if dp not in ["BEEP", "RANDOM"]:
                            continue
                        if profiler not in histograms[n_prec_errors][pbem][dp]:
                            continue

                        # direct cov is the same for prec/posc profilers
                        # print("[DATA] n_prec_errors:", n_prec_errors, "pbem:", pbem, "dp:", dp, "profiler:", profiler, "pmf:", dict([(i, j) for i, j in zip(x, pmf) if j > 0]))

                        # compute the full conditional probability table for this case
                        bar_plot = True
                        for prec_error_count in sorted(nprec_to_axc_index):
                            axc_idx = nprec_to_axc_index[prec_error_count]
                            
                            if prec_error_count in histograms \
                                and pbem in histograms[prec_error_count] \
                                and dp in histograms[prec_error_count][pbem] \
                                and profiler in histograms[prec_error_count][pbem][dp]:

                                x = []
                                y = []
                                # for rber in sorted(rber_to_axc_idx):
                                # for round in range(128):
                                for max_n_max_simultaneous_errors in [1, 2, 3, 4, 5, 6]:
                                    # PMF of max_simultaneous_posc_error_count over the K given N
                                    #   - probability that an ECC word is susceptible to at least one K-bit error pattern given N
                                    
                                    # make sure "round_idx" is NOT zero-indexed; makes no sense to have '0 rounds UNLESS no errors at all
                                    # round_idx == 0 -> happens when there are no errors possible to begin with; not an interesting case; means no uncorrectable errors
                                    if max_n_max_simultaneous_errors > prec_error_count + 1:
                                        round_idx = 0
                                    else:
                                        round_idx = 1 # 1-indexing
                                        for r in range(128):
                                            hist = histograms[prec_error_count][pbem][dp][profiler][r] # r = 0-indexing
                                            cdf = np.cumsum(hist)
                                            k_errors_max = np.argmax(cdf >= percentile * cdf[-1])
                                            if k_errors_max <= max_n_max_simultaneous_errors:
                                                break
                                            round_idx += 1
                                    if round_idx < 129:
                                        x.append(max_n_max_simultaneous_errors)
                                        y.append(round_idx)
                                    elif bar_plot: # off the edge! takes longer than 128 rounds, but we don't know exactly how many
                                        x.append(max_n_max_simultaneous_errors)
                                        y.append(round_idx * 2) # off the top to show > 128

                                if not bar_plot: # LINEPLOT 
                                    axs[axr_idx][axc_idx].plot(x, y, color=profiler_colors[profiler], marker='.', ms=5)
                                    # axs[axr_idx][axc_idx].set_xscale('log')
                                    # axs[axr_idx][axc_idx].set_yscale('log')
                                else: # BARPLOT
                                    bar_width = 1.0 / 6
                                    bar_offset = {'naive' : -bar_width * 2, 'beep_posc' : -bar_width, 'harp_prec' : 0, 'harp_posc' : bar_width}
                                    axs[axr_idx][axc_idx].bar(np.array(x) + bar_offset[profiler], y, width=bar_width, color=profiler_colors[profiler], align='edge', edgecolor='black', linewidth=0.4)
                                print('[DEBUG] max_n_simultaneous_errors: prof:', profiler, 'pbem:', pbem, 'dp:', dp, 'n_prec_errs:', prec_error_count, "percentile", percentile, x, y)

            labels_profilers = []
            for profiler in ['naive', 'beep_posc', 'harp_prec', 'harp_posc']:
                labels_profilers.append(mpatches.Patch(facecolor=profiler_colors[profiler], label=indirect_profiler_labels[profiler]))
            labels_pbem = []
            for pbem, pbem_label in reversed([("PBEM_100", "100%"), ("PBEM_75", "75%"), ("PBEM_50", "50%"), ("PBEM_25", "25%")]):
                labels_pbem.append(mlines.Line2D([], [], color='black', linestyle='-', marker=pbem_to_marker[pbem], ms=markersize[pbem], label=pbem_label))

            for i in range(len(axs)):
                axs[0][i].set_xlabel(ax_index_to_label[i])
                axs[0][i].xaxis.set_label_position('top') 

                # axs[3][i].set_xlim([0.5, 6.5])
                # axs[3][i].set_xticks([1, 2, 3, 4, 5, 6])
                
                xlabels = [1, 2, 3, 4, 5, 6]
                xtickpos = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
                axs[3][i].set_xticks(xlabels)
                axs[3][i].set_xticklabels(xlabels)
                axs[3][i].xaxis.set_minor_locator(mtick.FixedLocator(xtickpos))
                axs[3][i].tick_params(axis='x', which='major', length=0)
                axs[3][i].tick_params(axis='x', which='minor', length=12)
                axs[3][i].set_xlim([0.5, 6.5])

                axs[3][i].set_ylim([-8, 128])
                axs[3][i].set_yticks([0, 32, 64, 96])

                # axs[3][i].set_ylim([1, 138])
                # axs[3][i].set_yscale('log')
                # axs[3][i].set_yticks([1, 4, 16, 64])

            fig.text(0.50, 0.93, 'Number of Pre-Correction Errors per ECC Word', ha='center')
            fig.text(0.50, 0.04, 'Maximum Number of Simultaneous Post-Correction Errors Possible', ha='center')
            fig.text(0.03, 0.19, 'Number of Profiling Rounds', ha='center', rotation=90)
            fig.text(0.96, 0.28, 'Pre-Correction Error\nPer-Bit Probability', ha='center', rotation=270)
            fig.text(0.92, 0.71, '100%', ha='center', rotation=270)
            fig.text(0.92, 0.55, '75%', ha='center', rotation=270)
            fig.text(0.92, 0.38, '50%', ha='center', rotation=270)
            fig.text(0.92, 0.20, '25%', ha='center', rotation=270)

            # axs[0][0].set_yscale('log')
            # axs[0][0].set_ylabel('Number of ECC Words')
            # axs[3][0].legend(handles=list(labels_pbem), labelspacing=0.1, fontsize=6, ncol=1, loc='upper right', title="Pre-Correction Error\n  Per-Bit Probability", title_fontsize=6)
            axs[0][0].legend(handles=list(labels_profilers), labelspacing=0.01, fontsize=6, ncol=1, loc='upper right', frameon=False, borderpad=0.02)
            # axs[1][0].legend(handles=list(labels_dps), labelspacing=0.1, fontsize=6, ncol=1, loc='lower right')
            fig.tight_layout()
            # fig.subplots_adjust(wspace=0.1, hspace=0, bottom=0.16, top=0.86, left=0.06, right=0.94)
            fig.subplots_adjust(wspace=0.1, hspace=0, bottom=0.16, top=0.85, left=0.10, right=0.90)

            fig.canvas.manager.set_window_title("Figure 9b: P" + str(percentile) + " incomplete coverage after profiling - K: " + str(K) + " subfigure(b)")
            if not interactive:
                fname = 'fig_9b_incomplete_covereage_P' + str(percentile) + '_after_profiling.pdf' 
                if savefig_basename != None:
                    fname = savefig_basename + fname
                print("[INFO] Saving figure:", fname)
                fig.savefig(fname)

# ptable2d is assumed to be ptable2d[# precorrection errors][# postcorrection errors]
def calculate_eber_given_rber(K, N, rber, empirical_P_x_posc_ARBs_given_N_prec_ARBs, skipzero=True):
    ev = 0
    ev_err = 0
    for posc_error_count in range(0, K):
        logP_inner_mass_list = []
        logP_inner_error_list = []
        for prec_error_count in range(0, max(empirical_P_x_posc_ARBs_given_N_prec_ARBs.keys()) if skipzero else N):
            logP_of_n_raw_errors = spsta.binom.logpmf(prec_error_count, N, rber)
            if prec_error_count in empirical_P_x_posc_ARBs_given_N_prec_ARBs and \
                posc_error_count in empirical_P_x_posc_ARBs_given_N_prec_ARBs[prec_error_count]:
                logP_inner_mass_list.append(np.log(empirical_P_x_posc_ARBs_given_N_prec_ARBs[prec_error_count][posc_error_count]) + logP_of_n_raw_errors)
            else:
                logP_inner_error_list += logP_of_n_raw_errors # presume that P[K-missed bit pattern | N raw errors] = 1.0 (absolute worst case -> should really take remaining probability :))

        # TODO: perform the multiply-add in the log domain, if possible
        logP_inner_mass = np.logaddexp.reduce(logP_inner_mass_list)
        logP_inner_error = np.logaddexp.reduce(logP_inner_error_list)
        ev += posc_error_count * np.exp(logP_inner_mass) 
        ev_err += posc_error_count * np.exp(logP_inner_error)
    return ev, ev_err

# by law of total probability:
#   P[fail in 64ms window] = 1 - (P[0 raw] + P[1 raw]*P[corrected] + P[2 raw]* ...)
def plot_figure10_effective_bit_error_rate(raw_data, K, interactive, savefig_basename):
    N = K + compute_hamming_code_n_parity_bits(K)

    print("[INFO] counting samples to generate conditional probabilities P[fail_profiling | n_prec_errs]")
    aggregate_ptables = {}
    for fname in raw_data:
        for json_cfg_fname in raw_data[fname]:
            if json_database[json_cfg_fname]['k'] != K:
                continue
            
            for n_prec_errors in raw_data[fname][json_cfg_fname]:
                for word in raw_data[fname][json_cfg_fname][n_prec_errors]:
                    for pbem in word['ptables']:
                        for dp in word['ptables'][pbem]:
                            if dp not in ["BEEP", "BEEP2", "RANDOM"]: # one DP per profiler - nbd
                                continue
                            for profiler in word['ptables'][pbem][dp]:
                                if profiler not in aggregate_ptables:
                                    aggregate_ptables[profiler] = {}
                                if pbem not in aggregate_ptables[profiler]:
                                    aggregate_ptables[profiler][pbem] = {round_idx : {} for round_idx in range(128)}

                                assert len(word['ptables'][pbem][dp][profiler]) > 0
                                ptable_entry_keys = set(word['ptables'][pbem][dp][profiler].keys())
                                cur_ptable_entry = word['ptables'][pbem][dp][profiler][-1] if -1 in ptable_entry_keys else {} # primarily for BEEP2
                                for round_idx in range(128):
                                    if round_idx in ptable_entry_keys:
                                        cur_ptable_entry = word['ptables'][pbem][dp][profiler][round_idx]
                                    for n_missed_posc_errors in cur_ptable_entry:
                                        if n_missed_posc_errors not in aggregate_ptables[profiler][pbem][round_idx]:
                                            aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors] = {}
                                        for n_precorrection_errors in cur_ptable_entry[n_missed_posc_errors]:
                                            if round_idx == 127 and 'harp' in profiler and n_missed_posc_errors >= 2:
                                                print('[WARN] harp incomplete direct profile found!', "fname:", fname, "json_cfg_fname:", json_cfg_fname, "n_prec_errors:", n_prec_errors, "pbem:", pbem, "dp:", dp, "profiler:", profiler, "round_idx:", round_idx, "n_missed_posc_errors:", n_missed_posc_errors, "n_precorrection_errors:", n_precorrection_errors)
                                            if n_precorrection_errors not in aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors]:
                                                aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors][n_precorrection_errors] = {}
                                            for term, wordcount in cur_ptable_entry[n_missed_posc_errors][n_precorrection_errors].items():
                                                if term not in aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors][n_precorrection_errors]:
                                                    aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors][n_precorrection_errors][term] = 0
                                                aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors][n_precorrection_errors][term] += wordcount
    
    # print the probability tables
    # for profiler in aggregate_ptables:
    #     for pbem in aggregate_ptables[profiler]:
    #         for round_idx in aggregate_ptables[profiler][pbem]:
    #             print(profiler, pbem, round_idx, aggregate_ptables[profiler][pbem][round_idx])
                                    
    # evaluate the probability tables for each PBEM
    pbem_to_R = {'PBEM_100' : 1, 'PBEM_75' : 0.75, 'PBEM_50' : 0.5, 'PBEM_25' : 0.25}
    evaluated_probabilities = {}
    for profiler in aggregate_ptables:
        evaluated_probabilities[profiler] = {}
        for pbem in aggregate_ptables[profiler]:
            R = pbem_to_R[pbem]
            evaluated_probabilities[profiler][pbem] = {}
            for round_idx in aggregate_ptables[profiler][pbem]:
                evaluated_probabilities[profiler][pbem][round_idx] = {}
                for n_missed_posc_errors in aggregate_ptables[profiler][pbem][round_idx]:
                    for n_precorrection_errors in aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors]:
                        term_sum = 0
                        for term, count in aggregate_ptables[profiler][pbem][round_idx][n_missed_posc_errors][n_precorrection_errors].items():
                            term_exponents_raw = term[1:].split('o')
                            term_sum += count * math.pow(1-R, int(term_exponents_raw[0])) * math.pow(R, int(term_exponents_raw[1]))

                        if n_precorrection_errors not in evaluated_probabilities[profiler][pbem][round_idx]:
                            evaluated_probabilities[profiler][pbem][round_idx][n_precorrection_errors] = {}
                        assert n_missed_posc_errors not in evaluated_probabilities[profiler][pbem][round_idx][n_precorrection_errors]
                        evaluated_probabilities[profiler][pbem][round_idx][n_precorrection_errors][n_missed_posc_errors] = term_sum
    
    # binarize the probability counts according to the secondary ECC: <=1 error are "no error", >=2 are "error"
    print("Fraction of ECC words with incomplete coverage P[uncorrectable error pattern | n precorrection errors]:")
    binarized_probabilities = {}
    for profiler in evaluated_probabilities:
        binarized_probabilities[profiler] = {}
        for pbem in evaluated_probabilities[profiler]:
            binarized_probabilities[profiler][pbem] = {}
            for round_idx in evaluated_probabilities[profiler][pbem]:
                binarized_probabilities[profiler][pbem][round_idx] = {}
                for n_precorrection_errors in evaluated_probabilities[profiler][pbem][round_idx]:
                    bdist = {0 : 0, 1 : 0} # binary distribution
                    for n_missed_posc_errors in evaluated_probabilities[profiler][pbem][round_idx][n_precorrection_errors]:
                        if n_missed_posc_errors <= 1:
                            bdist[0] += \
                                evaluated_probabilities[profiler][pbem][round_idx][n_precorrection_errors][n_missed_posc_errors]
                        else: 
                            bdist[1] += \
                                evaluated_probabilities[profiler][pbem][round_idx][n_precorrection_errors][n_missed_posc_errors]
                    binarized_probabilities[profiler][pbem][round_idx][n_precorrection_errors] = bdist[1] / float(sum(bdist.values()))
                if round_idx == 127:
                    print("    ", profiler, pbem, round_idx, binarized_probabilities[profiler][pbem][round_idx])

    # let K := number of times that an "uncorrectable" event occurs within "some interval" (i.e., K can be 0, 1, 2, ...) 
    #       Note that we define an "uncorrectable" event as any number of errors that exceed secondary ECC's correction capability (i.e., nerrs >= 2)
    # Each refresh window is independent; the uncorrectable event doesn't depend on previous refresh windows
    # Events cannot occur simultaneously
    # ... then K ~ Poisson distributed 
    #
    # We can calculate:
    #   P[K=k] -> k failures occur in a given interval (e.g., 1B hours worth of refresh windows)
    # Given that:
    #   The average event rate (lambda) is "the average number of errors per refresh window" = E[distribution from above]

    # calculate actual numbers using RBERs
    profilers_to_plot = ['naive', 'beep_posc', 'harp_prec']
    round_indices_to_marker = [0, 1, 2, 4, 8, 16, 32, 64, 127]
    round_indices_to_plot = range(128) # [0, 1, 2, 4, 8, 16, 32, 64, 127]
    rbers = [1e-8, 1e-6, 1e-4]
    FIT_per_Mbit = {} # nfailures in 1 billion device hours - 25K-70K / Mbit
    P_failure = {} # probability that there is 2+ errors in an ECC word after profiling
    E_missed_bits_in_one_word_in_one_REFw = {} # expected number of missed words per Mbit
    for profiler in binarized_probabilities:
        if profiler not in profilers_to_plot: 
            continue
        P_failure[profiler] = {}
        for pbem in binarized_probabilities[profiler]:
            P_failure[profiler][pbem] = {}
            print("[INFO] calculating probability of failure after profiling with profiler:", profiler, "pbem:", pbem)
            # for round_idx in binarized_probabilities[profiler][pbem]:
            for round_idx in round_indices_to_plot:
                for rber in rbers:
                    effective_rber = rber # / pbem_to_R[pbem]
                    logP_mass_list = []
                    logP_error_list = []
                    for prec_error_count in range(0, 6):
                        logP_of_n_raw_errors = spsta.binom.logpmf(prec_error_count, N, effective_rber)
                        if prec_error_count in binarized_probabilities[profiler][pbem][round_idx]:
                            logP_mass_list.append(np.log(binarized_probabilities[profiler][pbem][round_idx][prec_error_count]) + logP_of_n_raw_errors)
                        else:
                            logP_error_list += logP_of_n_raw_errors # presume that P[failure | N raw errors] = 1.0 (absolute worst case -> should really take remaining probability :))

                    # TODO: convert to log-based computation
                    logP_mass = np.logaddexp.reduce(logP_mass_list)
                    logP_error = np.logaddexp.reduce(logP_error_list)
                    
                    # print(profiler, pbem, round_idx, rber, "E[missed bits per word in one refresh window] =", (outer_mass, outer_error))
                    if rber not in P_failure[profiler][pbem]:
                        P_failure[profiler][pbem][rber] = {}
                    P_failure[profiler][pbem][rber][round_idx] = np.exp(logP_mass)

    print("[INFO] generating plot...")
    rber_legend_entries = {}
    profiler_legend_entries = {}
    profiler_labels = {'naive' : "Naive", 'beep_posc' : "BEEP", 'harp_prec' : "HARP-U", 'harp_posc' : "HARP-A", "beep2_posc" : "HARP-A + BEEP"}
    profiler_colors = {'naive' : orange, 'beep_posc' : blue, 'beep2_posc' : purple, 'harp_prec' : green, 'harp_posc' : red}
    pbem_to_ax_index = {"PBEM_100" : 3, "PBEM_75" : 2, "PBEM_50" : 1, "PBEM_25" : 0}
    profiler_to_ax_index = {"naive" : 0, "beep_posc" : 1, "harp_posc" : 2, "harp_prec" : 2}
    ax_index_to_label = {0 : '25%', 1 : '50%', 2 : '75%', 3 : '100%'}
    markers = ['None', 'x', '.', '|']
    rber_to_marker = dict(zip(sorted(rbers), markers[:len(rbers)]))
    markersize = 6

    if 1: # merely scoping this section of code to avoid accidental reuse of these variables later on
        fig, axs = plt.subplots(1, 4, figsize=(5, 2.25), sharex=True, sharey=True)  
        # fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
        for profiler in P_failure:
            for pbem in P_failure[profiler]:
                ax_idx = pbem_to_ax_index[pbem]
                for rber in P_failure[profiler][pbem]:
                    round_idx, probability_mass = zip(*sorted(P_failure[profiler][pbem][rber].items()))
                    print(profiler, pbem, rber, "P[2+-bit error pattern] =", round_idx, probability_mass)
                    axs[ax_idx].plot(round_idx, probability_mass, color=profiler_colors[profiler]) #, markevery=0.2)
                    p = axs[ax_idx].plot([round_idx[i] for i in round_indices_to_marker], [probability_mass[i] for i in round_indices_to_marker]
                        , ls='None', color=profiler_colors[profiler], marker=rber_to_marker[rber], ms=markersize) #, markevery=0.2)
                    if rber not in rber_legend_entries:
                        rber_legend_entries[rber] = mlines.Line2D([], [], color='k'
                            , marker=p[0].get_marker(), ms=markersize, label='$10^{' + str(int(np.log10(rber))) + '}$')
                    if profiler not in profiler_legend_entries:
                        profiler_legend_entries[profiler] = mlines.Line2D([], [], color=p[0].get_color().flatten(), label=profiler_labels[profiler])


        legend1 = axs[-1].legend(handles=sorted(rber_legend_entries.values(), key=lambda x: x.get_label()), ncol=1 #, bbox_to_anchor=(1.075, 0.9)
            , labelspacing=0.1, fontsize=7, loc='lower left', title="RBER", title_fontsize=8)
        legend2 = axs[-2].legend(handles=sorted(profiler_legend_entries.values(), key=lambda x: x.get_label()), ncol=1 #, bbox_to_anchor=(1, 0.1)
            , labelspacing=0.1, fontsize=7, loc='lower left')
        # axs[-1].add_artist(legend1)
        # axs[-2].add_artist(legend2)

        for i in range(len(axs)):
            axs[i].set_xlabel(ax_index_to_label[i])
            axs[i].xaxis.set_label_position('top') 

            axs[i].set_yscale('log')
            try:
                axs[i].set_xscale('log', base=2)
            except:
                axs[i].set_xscale('log', basex=2)

            axs[i].tick_params(axis='y', which='minor', left=False, right=False)
            axs[i].tick_params(axis='x', which='minor', top=False, bottom=True)
            axs[i].xaxis.set_minor_locator(mtick.AutoMinorLocator())
            axs[i].minorticks_on()
            axs[i].set_xticks(np.array([np.array([1.25, 1.5, 1.75]) * n for n in 2**np.arange(0, 7)]).flatten(), minor=True)
            axs[i].xaxis.set_minor_formatter(mtick.NullFormatter())

            axs[i].set_xticks([1, 2, 4, 8, 16, 32, 64, 128])
            axs[i].set_xticklabels([1, None, 4, None, 16, None, 64, None])
            axs[i].set_xlim([0.75, 192])

        fig.text(0.58, 0.92, 'Per-Bit Probability of Pre-Correction Error', ha='center')    
        fig.text(0.56, 0.04, 'Number of Profiling Rounds', ha='center')    
        axs[0].set_ylabel('Probability of observing\n$\geq$2 errors post-repair')
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.1, bottom=0.22, top=0.82)
        fig.canvas.manager.set_window_title("Figure 10.1: P[2+ error patterns] per word per trefw")
        if not interactive:
            fname = 'fig_10_1_probability_of_uncorrectable_error.pdf' 
            if savefig_basename != None:
                fname = savefig_basename + fname
            print("[INFO] Saving figure:", fname)
            fig.savefig(fname)

    ########################################################################
    # BER post-repair, pre-secondary-ECC
    # BER post-secondary-ECC (apply effects of secondary ECC)
    ########################################################################
    # normalize the distributions for each 'n_precorrection_errors' value
    pdist_pre_secondary_ecc = evaluated_probabilities
    print("[INFO] normalizing the distribution of P[X-bit missed post-correction pattern | N raw bit errors] over X for each N")
    for profiler in pdist_pre_secondary_ecc:
        for pbem in pdist_pre_secondary_ecc[profiler]:
            for round_idx in pdist_pre_secondary_ecc[profiler][pbem]:
                for n_precorrection_errors in pdist_pre_secondary_ecc[profiler][pbem][round_idx]:
                    distribution = pdist_pre_secondary_ecc[profiler][pbem][round_idx][n_precorrection_errors]
                    total_sum = sum(distribution.values()) 
                    if total_sum < 1e-15: # tolerance
                        normalized_distribution = {0 : 1.0}
                    else:
                        normalized_distribution = {i : j / total_sum for i, j in distribution.items()}
                    pdist_pre_secondary_ecc[profiler][pbem][round_idx][n_precorrection_errors] = normalized_distribution

            print("    ", profiler, pbem)
            for round_idx in [0, 31, 63, 127]:
                print("        ", round_idx, "E[distribution]:", {i : sum([k*l for k, l in j.items()]) for i, j in pdist_pre_secondary_ecc[profiler][pbem][round_idx].items()})
                for n_precorrection_errors in pdist_pre_secondary_ecc[profiler][pbem][round_idx]:
                    print("            n_prec_errors:", n_precorrection_errors, "distribution:", pdist_pre_secondary_ecc[profiler][pbem][round_idx][n_precorrection_errors])

    # adjust the probability counts according to the secondary ECC: any samples with 1 (correctable) error move to the '0' bucket
    # assume an SEC Hamming code
    # pre_ecc_dist is pre_ecc_dist[n_raw_errors][n_errrs_pre_secondary_ecc] = P[Y|X]
    # post_ecc_dist is post_ecc_dist[n_raw_errors][n_errrs_post_secondary_ecc] = P[Z|X]
    # transform is transform[n_errrs_pre_secondary_ecc][n_errrs_post_secondary_ecc] = P[Z|Y]
    def transform_dist(pre_ecc_dist):
        # {i : {k : l / sum(j.values()) for k, l in j.items()} for i, j in a.items()}
        transform = {0: {0: 1.0}, 1: {0: 1.0}, 2: {0: 0.00397, 1: 0.09392, 2: 0.500023, 3: 0.402087}, 3: {0: 0.000293, 1: 0.012106, 2: 0.135846, 3: 0.496422, 4: 0.355333}, 4: {0: 1.7e-05, 1: 0.001889, 2: 0.029969, 3: 0.191712, 4: 0.470278, 5: 0.306135}, 5: {0: 1e-06, 1: 0.000184, 2: 0.004401, 3: 0.047213, 4: 0.225531, 5: 0.452494, 6: 0.270176}, 6: {1: 1.6e-05, 2: 0.000459, 3: 0.008363, 4: 0.066086, 5: 0.255339, 6: 0.433011, 7: 0.236726}, 7: {2: 4.1e-05, 3: 0.001046, 4: 0.013584, 5: 0.086816, 6: 0.278132, 7: 0.412242, 8: 0.208139}, 8: {2: 4e-06, 3: 0.000103, 4: 0.002065, 5: 0.02027, 6: 0.107588, 7: 0.297611, 8: 0.390341, 9: 0.182018}, 9: {3: 8e-06, 4: 0.00025, 5: 0.003611, 6: 0.029089, 7: 0.129718, 8: 0.311364, 9: 0.367382, 10: 0.158578}, 10: {4: 2e-05, 5: 0.000454, 6: 0.0055, 7: 0.039073, 8: 0.150736, 9: 0.322639, 10: 0.34341, 11: 0.138168}}
        post_ecc_dist = {}
        for n_raw_errors in pre_ecc_dist: # within the 71/136-bit raw memory bits
            if n_raw_errors not in post_ecc_dist:
                post_ecc_dist[n_raw_errors] = {}
            for n_posc_errors in range(20): # within the 57/120-bit secondary ECC dataword 
                # P[Z|X] = SUM_Y P[Z,Y|X] = SUM_Y P[Z|Y,X]*P[Y|X]  (by chain rule)
                # HOWEVER -> We only have P[Z|Y], and assuming P[Z|Y] == P[Z|Y,X] is assuming conditional independence (Z \perp Y | X)
                P_n_posc_errors_given_n_raw_errors = 0
                for n_prec_errors in pre_ecc_dist[n_raw_errors]: # within the 64/128-bit on-die ECC dataword
                    if n_raw_errors in transform and n_posc_errors in transform[n_prec_errors]: # else the contribution is 0 :)
                        P_n_posc_errors_given_n_raw_errors += pre_ecc_dist[n_raw_errors][n_prec_errors] * transform[n_prec_errors][n_posc_errors]
                post_ecc_dist[n_raw_errors][n_posc_errors] = P_n_posc_errors_given_n_raw_errors
        return post_ecc_dist

    def transform_dist_naive(pre_ecc_dist):
        post_ecc_dist = {}
        for n_precorrection_errors in pre_ecc_dist:
            post_ecc_dist[n_precorrection_errors] = {0 : 0}
            for n_posc_errors in pre_ecc_dist[n_precorrection_errors]:
                if n_posc_errors not in post_ecc_dist[n_precorrection_errors]:
                    post_ecc_dist[n_precorrection_errors][n_posc_errors] = 0
                if n_posc_errors == 1:
                    post_ecc_dist[n_precorrection_errors][0] += pre_ecc_dist[n_precorrection_errors][n_posc_errors]
                else:
                    post_ecc_dist[n_precorrection_errors][n_posc_errors] += pre_ecc_dist[n_precorrection_errors][n_posc_errors]
        return post_ecc_dist


    # calculate the post-secondary-ECC distribution
    pdist_post_secondary_ecc = {}
    for profiler in pdist_pre_secondary_ecc:
        pdist_post_secondary_ecc[profiler] = {}
        for pbem in pdist_pre_secondary_ecc[profiler]:
            pdist_post_secondary_ecc[profiler][pbem] = {}
            for round_idx in pdist_pre_secondary_ecc[profiler][pbem]:
                pre_ecc_dist = pdist_pre_secondary_ecc[profiler][pbem][round_idx]
                post_ecc_dist = transform_dist(pre_ecc_dist)
                pdist_post_secondary_ecc[profiler][pbem][round_idx] = post_ecc_dist
                # print("[DEBUG] transform for profiler:", profiler, "pbem:", pbem, "round_idx:", round_idx, ":")
                # print("        pre: ", pre_ecc_dist)
                # print("        post:", post_ecc_dist)

    # compute the BER (i.e., expectation values) for each configuration and rber
    profilers_to_plot = ['naive', 'beep_posc', 'harp_prec', 'harp_posc']
    round_indices_to_marker = [0, 1, 2, 4, 8, 16, 32, 64, 127]
    round_indices_to_plot = range(128)
    rbers = [1e-8, 1e-6, 1e-4]

    BER_pre_secondary_ECC = {} # expected number of missed bits in a single ECC word (basically the BER)
    for profiler in pdist_pre_secondary_ecc:
        if profiler not in profilers_to_plot: 
            continue
        BER_pre_secondary_ECC[profiler] = {}
        for pbem in pdist_pre_secondary_ecc[profiler]:
            BER_pre_secondary_ECC[profiler][pbem] = {}
            print("[INFO] calculating E[missed bits / word] after profiling with profiler:", profiler, "pbem:", pbem)
            for round_idx in round_indices_to_plot:
                for rber in rbers:
                    effective_rber = rber # / pbem_to_R[pbem]
                    ptable = pdist_pre_secondary_ecc[profiler][pbem][round_idx]
                    E_n_posc_errors_per_word, e_mass = calculate_eber_given_rber(K, N, effective_rber, ptable)

                    if rber not in BER_pre_secondary_ECC[profiler][pbem]:
                        BER_pre_secondary_ECC[profiler][pbem][rber] = {}
                    BER_pre_secondary_ECC[profiler][pbem][rber][round_idx] = E_n_posc_errors_per_word / float(K)

                    print("[DEBUG] profiler:", profiler, "pbem:", pbem, "nrounds:", round_idx
                        , "E[# errors pre-secondary-ECC | rber=" + str(rber) + "] =", E_n_posc_errors_per_word)

    profilers_to_plot = ['naive', 'beep_posc', 'harp_prec'] # HARP-PREC is identical to HARP-POSC here
    round_indices_to_marker = [0, 1, 2, 4, 8, 16, 32, 64, 127]
    round_indices_to_plot = range(128)
    rbers = [1e-8, 1e-6, 1e-4]

    BER_post_secondary_ECC = {} # expected number of missed bits in a single ECC word (basically the BER)
    for profiler in pdist_post_secondary_ecc:
        if profiler not in profilers_to_plot: 
            continue
        BER_post_secondary_ECC[profiler] = {}
        for pbem in pdist_post_secondary_ecc[profiler]:
            BER_post_secondary_ECC[profiler][pbem] = {}
            print("[INFO] calculating E[missed bits / word] after profiling with profiler:", profiler, "pbem:", pbem)
            for round_idx in round_indices_to_plot:
                for rber in rbers:
                    effective_rber = rber # / pbem_to_R[pbem]
                    ptable = pdist_post_secondary_ecc[profiler][pbem][round_idx]
                    E_n_posc_errors_per_word, e_mass = calculate_eber_given_rber(K, N, effective_rber, ptable)

                    if rber not in BER_post_secondary_ECC[profiler][pbem]:
                        BER_post_secondary_ECC[profiler][pbem][rber] = {}
                    BER_post_secondary_ECC[profiler][pbem][rber][round_idx] = E_n_posc_errors_per_word / float(K)

                    print("[DEBUG] profiler:", profiler, "pbem:", pbem, "nrounds:", round_idx
                        , "E[# errors post-secondary-ECC | rber=" + str(rber) + "] =", E_n_posc_errors_per_word)

    print("[INFO] generating plot...")
    profiler_labels = {'naive' : "Naive", 'beep_posc' : "BEEP", 'harp_prec' : "HARP-U", 'harp_posc' : "HARP-A", "beep2_posc" : "HARP-A + BEEP"}
    profiler_colors = {'naive' : orange, 'beep_posc' : blue, 'beep2_posc' : purple, 'harp_prec' : green, 'harp_posc' : red}
    pbem_to_ax_index = {"PBEM_100" : 3, "PBEM_75" : 2, "PBEM_50" : 1, "PBEM_25" : 0}
    profiler_to_ax_index = {"naive" : 0, "beep_posc" : 1, "harp_posc" : 2, "harp_prec" : 2}
    ax_index_to_label = {0 : '25%', 1 : '50%', 2 : '75%', 3 : '100%'}
    markers = ['None', 'x', '.', '|']
    rber_to_marker = dict(zip(sorted(rbers), markers[:len(rbers)]))
    markersize = 6

    for data_stream, y_label, plot_fname, window_title in [
          (BER_pre_secondary_ECC, "BER Before\nReactive Profiling", 'fig_10_2_ber_before_secondary_ecc.pdf', "Figure 10.2: BER before secondary ECC")
        , (BER_post_secondary_ECC, "BER After\nReactive Profiling", 'fig_10_3_ber_after_secondary_ecc.pdf', "Figure 10.3: BER after secondary ECC")]:
        rber_legend_entries = {}
        profiler_legend_entries = {}
    
        fig, axs = plt.subplots(1, 4, figsize=(5, 2.25), sharex=True, sharey=True)  
        # fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
        for profiler in data_stream:
            for pbem in data_stream[profiler]:
                ax_idx = pbem_to_ax_index[pbem]
                for rber in data_stream[profiler][pbem]:
                    round_idx, BERs = zip(*sorted(data_stream[profiler][pbem][rber].items()))
                    print(plot_fname, profiler, pbem, rber, "E[missed bits per word in one refresh window] =", round_idx, BERs)
                    axs[ax_idx].plot(round_idx, BERs, color=profiler_colors[profiler]) #, markevery=0.2)
                    p = axs[ax_idx].plot([round_idx[i] for i in round_indices_to_marker], [BERs[i] for i in round_indices_to_marker]
                        , ls='None', color=profiler_colors[profiler], marker=rber_to_marker[rber], ms=markersize) #, markevery=0.2)
                    if rber not in rber_legend_entries:
                        rber_legend_entries[rber] = mlines.Line2D([], [], color='k'
                            , marker=p[0].get_marker(), ms=markersize, label='$10^{' + str(int(np.log10(rber))) + '}$')
                    if profiler not in profiler_legend_entries:
                        profiler_legend_entries[profiler] = mlines.Line2D([], [], color=p[0].get_color().flatten(), label=profiler_labels[profiler])

        legend1 = axs[-1].legend(handles=sorted(rber_legend_entries.values(), key=lambda x: x.get_label()), ncol=1 #, bbox_to_anchor=(1.075, 0.9)
            , labelspacing=0.1, fontsize=7, loc='lower left', title="RBER", title_fontsize=8)
        legend2 = axs[-2].legend(handles=sorted(profiler_legend_entries.values(), key=lambda x: x.get_label()), ncol=1 #, bbox_to_anchor=(1, 0.1)
            , labelspacing=0.1, fontsize=7, loc='lower left')
        # axs[-1].add_artist(legend1)
        # axs[-2].add_artist(legend2)

        for i in range(len(axs)):
            axs[i].set_xlabel(ax_index_to_label[i])
            axs[i].xaxis.set_label_position('top') 

            axs[i].set_yscale('log')
            try:
                axs[i].set_xscale('log', base=2)
            except:
                axs[i].set_xscale('log', basex=2)

            axs[i].tick_params(axis='y', which='minor', left=False, right=False)
            axs[i].tick_params(axis='x', which='minor', top=False, bottom=True)
            axs[i].xaxis.set_minor_locator(mtick.AutoMinorLocator())
            axs[i].minorticks_on()
            axs[i].set_xticks(np.array([np.array([1.25, 1.5, 1.75]) * n for n in 2**np.arange(0, 7)]).flatten(), minor=True)
            axs[i].xaxis.set_minor_formatter(mtick.NullFormatter())

            axs[i].set_xticks([1, 2, 4, 8, 16, 32, 64, 128])
            axs[i].set_xticklabels([1, None, 4, None, 16, None, 64, None])
            axs[i].set_xlim([0.75, 192])

        fig.text(0.58, 0.92, 'Per-Bit Probability of Pre-Correction Error', ha='center')    
        fig.text(0.56, 0.04, 'Number of Profiling Rounds', ha='center')    
        axs[0].set_ylabel(y_label)
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.1, bottom=0.22, top=0.82)
        fig.canvas.manager.set_window_title(window_title)
        if not interactive:
            fname = plot_fname 
            if savefig_basename != None:
                fname = savefig_basename + fname
            print("[INFO] Saving figure:", fname)
            fig.savefig(fname)

def analyze_data(raw_data, interactive, savefig_basename):
    print("[INFO] analyzing raw_data")

    all_k = {json_database[json_cfg_fname]['k'] for fname in raw_data for json_cfg_fname in raw_data[fname]}
    if savefig_basename and os.path.isdir(savefig_basename):
        savefig_basename = os.path.join(os.path.normpath(savefig_basename), '')

    if True:
        print("") # newline
        print("[INFO] #######################################################")
        print("[INFO] Plot 1/4 - Figure 6,7: Coverage of direct errors, Bootstrapping")
        print("[INFO] #######################################################")
        try:
            for K in all_k:
                plot_figure6_direct_error_coverage_and_figure7_bootstrapping(raw_data, K, interactive, savefig_basename)
        except Exception as e:
            print(e)
            print(traceback.format_exc())

    if True:
        print("") # newline
        print("[INFO] #######################################################")
        print("[INFO] Plot 2/4 - Figure 8: Indirect errror coverage")
        print("[INFO] #######################################################")
        try:
            for K in all_k:
                plot_figure8_indirect_error_coverage(raw_data, K, interactive, savefig_basename)
        except Exception as e:
            print(e)
            print(traceback.format_exc())

    if True:
        print("") # newline
        print("[INFO] #######################################################")
        print("[INFO] Plot 3/4 - Figure 9a,9b: Impact of incomplete direct error coverage")
        print("[INFO] #######################################################")
        try:
            for K in all_k:
                plot_figure9ab_hisograms_of_simultaneous_postcorrection_errors(raw_data, K, interactive, savefig_basename)
        except Exception as e:
            print(e)
            print(traceback.format_exc())

    if True:
        print("") # newline
        print("[INFO] #######################################################")
        print("[INFO] Plot 4/4 - Figure 10: BER analysis")
        print("[INFO] #######################################################")
        try:
            for K in all_k:
                plot_figure10_effective_bit_error_rate(raw_data, K, interactive, savefig_basename)
        except Exception as e:
            print(e)
            print(traceback.format_exc())

    if interactive:
        plt.show()
    # print(raw_data)

def extract_bracelist(spt, start_pos, open_brace='{', close_brace='}'):
    assert spt[start_pos] == open_brace
    epos = start_pos + 1
    while epos < len(spt) and spt[epos] != close_brace:
        epos += 1
    assert epos < len(spt), "mismatched brace list!"
    return spt[start_pos + 1 : epos], epos

def load_json_cfg_file(fname):
    assert json_file_list != None, "[ERROR] JSON ECC configuration file database not initialized"    
    assert fname in json_file_list, "[ERROR] missing JSON config file: " + fname
    fpath = json_file_list[fname]
    with open(fpath) as f:
        json_data = json.load(f)
    return json_data

def parse_file(fname):
    data = {}
    with open(fname, 'r') as f:
        line_num = 1
        timeout_sample = False
        harp_p_augmentation_warning_given = False
        try:
            for line in f:
                if line.startswith("[INFO]"):
                    if line[7:].startswith("evaluating JSON ECC code configuration file"):
                        spt = line.strip().split(' ')
                        json_cfg_fname = os.path.basename(spt[7])
                        if json_cfg_fname not in json_database:
                            json_database[json_cfg_fname] = load_json_cfg_file(json_cfg_fname)
                            assert rseed == json_database[json_cfg_fname]['p'], "mismatch!"
                        assert rseed == p, "mismatch! missing output?"
                        assert json_cfg_fname not in data, "saw same ECC code twice"
                        assert len(data) == code_idx, "missing output?"
                        data[json_cfg_fname] = {}
                    elif line[7:].startswith("evaluating ECC code:"):
                        spt = line.strip().split(' ')
                        code_idx = int(spt[4])
                        rseed = int(spt[8])
                    elif line[7:].startswith("generating Hamming code"):
                        spt = line.strip().split(' ')
                        k = int(spt[6])
                        p = int(spt[8])
                    elif line[7:].startswith("n_errors:"):
                        spt = line.strip().split(' ')
                        n_prec_errors = int(spt[2])
                        data[json_cfg_fname][n_prec_errors] = []
                elif line.startswith("[DATA]"):
                    if line[7:].startswith("word:"):
                        spt = line.strip().split(' ')
                        # word_idx = int(spt[2]) # IGNORED - idx can be off when we skip words
                        all_prec_error_positions, end_idx = extract_bracelist(spt, 4)
                        all_posc_error_positions_charged_only, end_idx = extract_bracelist(spt, end_idx + 2)
                        all_posc_error_positions, end_idx = extract_bracelist(spt, end_idx + 2)
                        assert len(all_prec_error_positions) == n_prec_errors, "mismatch!"
                        # assert len(data[json_cfg_fname][n_prec_errors]) == word_idx, "missing output? exp: " + str(word_idx) + " got: " + str(len(data[json_cfg_fname][n_prec_errors]))
                        data[json_cfg_fname][n_prec_errors].append({
                              'prec' : [int(i) for i in all_prec_error_positions]
                            , 'posc_charged' : [int(i) for i in all_posc_error_positions_charged_only]
                            , 'posc' : [int(i) for i in all_posc_error_positions]
                            , 'profiles' : {}
                            , 'histograms' : {}
                            , 'ptables' : {}
                            , 'line_num' : line_num
                        })
                    elif line[7:].startswith("pbem:") or line[7:].startswith("profiles"):
                        spt = line.strip().split(' ')

                        if spt[1] == 'pbem:':
                            assert len(spt) > 6, "truncated output"
                            
                            # new sample - reset the timeout state
                            timeout_sample = False 

                            pbem = spt[2]
                            dp = spt[4]
                            n_rounds = spt[6] # unused

                            # extract profiles from each profiler
                            if spt[7] != 'profiles': # THIS IS A SAT TIMEOUT - skip the entire sample
                                timeout_sample = True
                            else:
                                if pbem not in data[json_cfg_fname][n_prec_errors][-1]['profiles']:
                                    data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem] = {}
                                if dp not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem]:
                                    data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp] = {}
                                sample_spt_idx = 8
                        elif spt[1] == 'profiles':
                            assert timeout_sample == True, "should never get here otherwise!"
                            sample_spt_idx = 2
                        else:
                            assert False, "unexpected output"
                        
                        if not timeout_sample:
                            if(spt[sample_spt_idx] == 'b:'):
                                beep_prec, end_idx = extract_bracelist(spt, sample_spt_idx + 1)
                                assert spt[end_idx + 1] == '+:'
                                beep_posc, end_idx = extract_bracelist(spt, end_idx + 2)
                                assert 'beep_prec' not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['beep_prec'] = tuple(tuple(map(int, i.split(','))) for i in beep_prec)
                                assert 'beep_posc' not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['beep_posc'] = tuple(tuple(map(int, i.split(','))) for i in beep_posc)
                            elif(spt[sample_spt_idx] == '2:'):
                                beep2_prec, end_idx = extract_bracelist(spt, sample_spt_idx + 1)
                                assert spt[end_idx + 1] == '+:'
                                beep2_posc, end_idx = extract_bracelist(spt, end_idx + 2)
                                assert 'beep2_prec' not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['beep2_prec'] = tuple(tuple(map(int, i.split(','))) for i in beep2_prec)
                                assert 'beep2_posc' not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['beep2_posc'] = tuple(tuple(map(int, i.split(','))) for i in beep2_posc)
                            elif(spt[sample_spt_idx] == 'n:'):
                                naive, end_idx = extract_bracelist(spt, sample_spt_idx + 1)
                                assert spt[end_idx + 1] == 'h:'
                                harp_prec, end_idx = extract_bracelist(spt, end_idx + 2)
                                assert spt[end_idx + 1] == '+:'
                                harp_posc, end_idx = extract_bracelist(spt, end_idx + 2)
                                assert 'naive' not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['naive'] = tuple(tuple(map(int, i.split(','))) for i in naive)
                                assert 'harp_prec' not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]
                                harp_prec_set = tuple(tuple(map(int, i.split(','))) for i in harp_prec)
                                data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['harp_prec'] = harp_prec_set
                                assert 'harp_posc' not in data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]
                                harp_posc_set_potentially_missing_harp_prec_set = tuple(tuple(map(int, i.split(','))) for i in harp_posc)
                                all_harp_bits = {bit_idx : round_idx for bit_idx, round_idx in harp_prec_set}
                                for bit_idx, round_idx in harp_posc_set_potentially_missing_harp_prec_set:
                                    if bit_idx in all_harp_bits:
                                        all_harp_bits[bit_idx] = min(round_idx, all_harp_bits[bit_idx])
                                    else:
                                        all_harp_bits[bit_idx] = round_idx
                                data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['harp_posc'] = tuple((bit_idx, round_idx) for bit_idx, round_idx in all_harp_bits.items())
                                if len(data[json_cfg_fname][n_prec_errors][-1]['profiles'][pbem][dp]['harp_posc']) != len(harp_posc_set_potentially_missing_harp_prec_set) and harp_p_augmentation_warning_given == False:
                                    print("[WARN] HARP+ set was missing HARP bits - suspect old-style dump - added HARP set to HARP+ set")
                                    harp_p_augmentation_warning_given = True
                                # print("[WARN] BUG: HARP+ MUST INCLUDE ALL BITS IN HARP FOR CORRECTNESS - YES, FALSE POSITIVES")
                            else:
                                assert False, "problem with line: " + str(line)
                    elif line[7:].startswith("histograms"):
                        spt = line.strip().split() # no arguments splits on any amount of whitespace 

                        # if SAT timeout for just the histogram, there will be no data here
                        # if timeout based on the profile computations, we don't add an entry at all
                        if len(spt) > 2 and not timeout_sample:
                            if pbem not in data[json_cfg_fname][n_prec_errors][-1]['histograms']:
                                data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem] = {}
                            if dp not in data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem]:
                                data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp] = {}

                            if spt[2] == '+:':
                                beep_posc_hist, end_idx = extract_bracelist(spt, 3)
                                assert 'beep_posc' not in data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]['beep_posc'] = tuple(tuple(map(int, i.split(','))) for i in beep_posc_hist)
                            elif spt[2] == '2:':
                                beep2_posc_hist, end_idx = extract_bracelist(spt, 3)
                                assert 'beep2_posc' not in data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]['beep2_posc'] = tuple(tuple(map(int, i.split(','))) for i in beep2_posc_hist)
                            elif spt[2] == 'n:':
                                naive_hist, end_idx = extract_bracelist(spt, 3)
                                assert spt[end_idx + 1] == 'h:'
                                harp_prec_hist, end_idx = extract_bracelist(spt, end_idx + 2)
                                assert spt[end_idx + 1] == 'h+:'
                                harp_posc_hist, end_idx = extract_bracelist(spt, end_idx + 2)
                                assert 'naive' not in data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]['naive'] = tuple(tuple(map(int, i.split(','))) for i in naive_hist)
                                assert 'harp_prec' not in data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]['harp_prec'] = tuple(tuple(map(int, i.split(','))) for i in harp_prec_hist)
                                assert 'harp_posc' not in data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]
                                data[json_cfg_fname][n_prec_errors][-1]['histograms'][pbem][dp]['harp_posc'] = tuple(tuple(map(int, i.split(','))) for i in harp_posc_hist)
                            else:
                                assert False, "unexpected data in line: " + str(spt)
                    elif line[7:].startswith("P[x|n]"):
                        if pbem not in data[json_cfg_fname][n_prec_errors][-1]['ptables']:
                            data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem] = {}
                        if dp not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem]:
                            data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp] = {}

                        spt = line.strip().split() # no arguments splits on any amount of whitespace 
                        assert spt[2] == "profiler:", "Corrupted output!"
                        profiler = spt[3]
                        assert spt[4] == "n_prec_errors:", "Corrupted output!"
                        n_prec_errors = int(spt[5])
                        assert spt[6] == "within_round_idx:", "Corrupted output!"
                        round_idx = int(spt[7])
                        if round_idx == 18446744073709551615: # -1u64 represents BEEP2's initial knowledege (i.e., before round_idx=0) 
                            round_idx = -1
                        assert spt[8] == "known_bits:", "Corrupted output!"
                        known_bits, end_idx = extract_bracelist(spt, 9)
                        assert spt[end_idx + 1] == 'log2(normalizer):', "Corrupted output!"
                        normalizer = spt[end_idx + 2]
                        assert spt[end_idx + 3] == 'data:', "Corrupted output!"
                        datapoints = spt[end_idx + 4:]
                        raw_probability_table = {}
                        for entry in datapoints:
                            assert entry[:2] == 'x=', "Corrupted output!"
                            n_missed_posc_errors_string, data_string = entry[:-1].split('{')
                            n_missed_posc_errors = int(n_missed_posc_errors_string[2:])

                            n_prec_errors_entries = data_string.replace(']','[').split('[')
                            for entry_idx in range(len(n_prec_errors_entries) // 2):
                                n = int(n_prec_errors_entries[2 * entry_idx + 0][:-1])
                                entry = [i.split(':') for i in n_prec_errors_entries[2 * entry_idx + 1].split(',')[:-1]]
                                if n_missed_posc_errors not in raw_probability_table:
                                    raw_probability_table[n_missed_posc_errors] = {}
                                raw_probability_table[n_missed_posc_errors][n] = {i[0] : int(i[1]) for i in entry}
                        
                        if profiler == 'b':
                            if 'beep_posc' not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]:
                                data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['beep_posc'] = {}
                            assert round_idx not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['beep_posc']
                            data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['beep_posc'][round_idx] = raw_probability_table
                        elif profiler == '2':
                            if 'beep2_posc' not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]:
                                data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['beep2_posc'] = {}
                            assert round_idx not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['beep2_posc']
                            data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['beep2_posc'][round_idx] = raw_probability_table
                        elif profiler == 'n':
                            if 'naive' not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]:
                                data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['naive'] = {}
                            assert round_idx not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['naive']
                            data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['naive'][round_idx] = raw_probability_table
                        elif profiler == 'h':
                            if 'harp_prec' not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]:
                                data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['harp_prec'] = {}
                            assert round_idx not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['harp_prec']
                            data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['harp_prec'][round_idx] = raw_probability_table
                        elif profiler == '+':
                            if 'harp_posc' not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]:
                                data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['harp_posc'] = {}
                            assert round_idx not in data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['harp_posc']
                            data[json_cfg_fname][n_prec_errors][-1]['ptables'][pbem][dp]['harp_posc'][round_idx] = raw_probability_table
                        else:
                            assert False, "unexpected data in line: " + str(spt)

                line_num += 1
        except:
            # print("[DEBUG] data:", data)
            print("[ERROR] failed parsing at fname:", fname, "line:", line_num)
            raise

    # post-process the BEEP2 profile values so that they follow the HARP+ timeline
    for json_cfg_fname in data:
        beep2_warning_random_suppressed = False
        beep2_warning_harp_suppressed = False
        K = json_database[json_cfg_fname]['k']
        for n_prec_errors in data[json_cfg_fname]:
            for word in data[json_cfg_fname][n_prec_errors]:
                prec_set = set(word['prec'])
                posc_set = set(word['posc'])
                direct_set = {i for i in posc_set if i in prec_set}
                for pbem in word['profiles']:
                    if 'BEEP2' in word['profiles'][pbem]:
                        if 'RANDOM' not in word['profiles'][pbem]:
                            if not beep2_warning_random_suppressed: 
                                print("[WARN] cannot compute profiling rounds for BEEP2 without RANDOM pattern tested! pbem:", pbem, "prec:", prec_set, "posc:", posc_set)
                                print("[WARN] this warning suppressed in future")
                                beep2_warning_random_suppressed = True
                            del word['profiles'][pbem]['BEEP2']
                            continue
                        if 'harp_prec' not in word['profiles'][pbem]['RANDOM']:
                            if not beep2_warning_harp_suppressed: 
                                print("[WARN] cannot compute profiling rounds for BEEP2 without HARP present also! pbem:", pbem, "prec:", prec_set, "posc:", posc_set)
                                print("[WARN] this warning suppressed in future")
                                beep2_warning_harp_suppressed = True
                            del word['profiles'][pbem]['BEEP2']
                            continue
                        
                        # figure out how many rounds the same config took HARP to achive full coverage of direct errors
                        all_bit_idxs_found = {bit_idx for bit_idx, round_idx in word['profiles'][pbem]['RANDOM']['harp_prec']}
                        if len(direct_set) == 0:
                            n_rounds_for_full_harp_direct_coverage = 0
                        elif all_bit_idxs_found == direct_set:
                            n_rounds_for_full_harp_direct_coverage = max([round_idx for bit_idx, round_idx in word['profiles'][pbem]['RANDOM']['harp_prec']])
                        else:
                            print("WARN: incomplete HARP direct coverage:", direct_set, word['profiles'][pbem]['RANDOM']['harp_prec'])
                            n_rounds_for_full_harp_direct_coverage = 128

                        # update BEEP2 round counts accordingly
                        new_bracelist = []
                        for bit_idx, round_idx in word['profiles'][pbem]['BEEP2']['beep2_prec']:
                            if round_idx == 18446744073709551615: # (uint64_t)-1
                                round_idx = 0
                            round_idx += n_rounds_for_full_harp_direct_coverage
                            new_bracelist.append((bit_idx, round_idx))
                        word['profiles'][pbem]['BEEP2']['beep2_prec'] = tuple(new_bracelist)
                        new_bracelist = []
                        for bit_idx, round_idx in word['profiles'][pbem]['BEEP2']['beep2_posc']:
                            if round_idx == 18446744073709551615: # (uint64_t)-1
                                round_idx = 0
                            round_idx += n_rounds_for_full_harp_direct_coverage
                            new_bracelist.append((bit_idx, round_idx))
                        word['profiles'][pbem]['BEEP2']['beep2_posc'] = tuple(new_bracelist)
                        # print(word['profiles'][pbem]['BEEP2'])

                        # update BEEP2 histograms accordingly
                        if 'BEEP2' in word['histograms'][pbem]:
                            new_bracelist = []
                            for round_idx, n_simultaneous_errors in word['histograms'][pbem]['BEEP2']['beep2_posc']:
                                new_bracelist.append((round_idx + n_rounds_for_full_harp_direct_coverage, n_simultaneous_errors))
                            word['histograms'][pbem]['BEEP2']['beep2_posc'] = tuple(new_bracelist)
                            # print(word['histograms'][pbem]['BEEP2']['beep2_posc'])
                        else:
                            print("[WARN] found profile but NOT histogram for BEEP2")


    # add the 'ALL' data pattern that swtiches between the individual patterns:
    # {RANDOM, ALL_ONES, COLSTRIPE}
    n_words = 0
    for json_cfg_fname in data:
        K = json_database[json_cfg_fname]['k']
        for n_prec_errors in data[json_cfg_fname]:
            for word in data[json_cfg_fname][n_prec_errors]:
                prec_set = set(word['prec'])
                posc_set = set(word['posc'])
                n_words += 1
                for pbem in word['profiles']:
                    assert 'COMBINED' not in word['profiles'][pbem], "wat. when did we start precomputing/running this elsewhere?"
                    combined_profiles = {}
                    for profiler in {'naive', 'harp_prec', 'harp_posc'}:
                        missing_sample = False
                        cur_profile = []
                        for combined_round_idx in range(128):
                            single_pattern_round_idx = combined_round_idx // 3
                            dp = ['RANDOM', 'ALL_ONES', 'COLSTRIPE'][combined_round_idx % 3]
                            if pbem not in word['profiles'] or dp not in word['profiles'][pbem]: # this sample is missing
                                missing_sample = True
                                break
                            for bit_idx, round_idx in word['profiles'][pbem][dp][profiler]:
                                if round_idx == single_pattern_round_idx:
                                    cur_profile.append((bit_idx, combined_round_idx))
                        if not missing_sample:
                            all_bits = {i for i, j in cur_profile}
                            combined_profiles[profiler] = tuple(min([(i, j) for i, j in cur_profile if i == this_bit], key=lambda x: x[1]) for this_bit in all_bits)
                    word['profiles'][pbem]['COMBINED'] = combined_profiles
                                                    

    # print("[DEBUG]", data)
    print("[DEBUG] fname:", fname, "data size (MiB):", total_size(data) / 1048576)
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

def main(interactive, savefig_basename, input_files_and_dirs):
    # sanitize input filenames
    validated_fpaths = []
    for f_or_d in input_files_and_dirs:
        if os.path.isfile(f_or_d):
            validated_fpaths.append(f_or_d)
        elif os.path.isdir(f_or_d):
            for root, dirs, files in os.walk(f_or_d):
                validated_fpaths += [os.path.join(root, f) for f in files]
        else:
            print("[ERROR] invalid file/dir:", f_or_d)
            sys.exit(-1)

    # segregate data and JSON files
    cleaned_fpaths = []
    for fpath in validated_fpaths:
        if fpath.endswith('.json'):
            fname = os.path.basename(fpath)
            json_file_list[fname] = fpath
        else:
            cleaned_fpaths.append(fpath)

    print("[INFO] parsing", len(cleaned_fpaths), "input data files given", len(json_file_list), "input JSON configuration files...")
    all_data = parse_files(cleaned_fpaths)
    analyze_data(all_data, interactive, savefig_basename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--interactive", help="draw interactive plot with matplotlib (else save plot to file)", action='store_true')
    parser.add_argument("-o", "--output-filename-base", help="custom filename base of plots to save - will be suffixed with -figN.pdf", type=str, default=None)
    parser.add_argument("input_files_and_dirs", metavar="input-files-and-dirs", help="JSON + data files (or directory containing files) to parse", nargs='+')
    args = parser.parse_args()

    if not display_valid and args.interactive:
        print("[ERROR] cannot create interactive plots without valid X server display")
        sys.exit(-1)
    main(args.interactive, args.output_filename_base, args.input_files_and_dirs)    
