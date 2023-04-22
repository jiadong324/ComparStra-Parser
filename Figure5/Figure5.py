#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/18

'''

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from matplotlib_venn import venn2

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Constant import *
from Helpers.Functions import get_survivor_supp, get_survivor_suppvec

sns.set_theme(style="ticks", font="Arial", font_scale=1.0)

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2

def plot_5a5b(figdir):

    df_wgs = pd.read_csv(f'{WORKDIR}/truvari/assembly_read_truvari.tsv', sep='\t', header=[0])

    stra_order = ['Assembly', 'Read']
    stra_color = [STRACOLORS[ele] for ele in stra_order]
    stra_legends = [Patch(label=ele, color=STRACOLORS[ele]) for ele in stra_order]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 3))
    sns.boxplot(data=df_wgs, x='recall', y='dataset', hue='strategy', hue_order=stra_order, palette=stra_color, ax=axes[0])
    axes[0].set_xlabel('Recall', fontsize=13)

    sns.boxplot(data=df_wgs, x='precision', y='dataset', hue='strategy', hue_order=stra_order, palette=stra_color, ax=axes[1])
    axes[1].set_xlabel('Precision', fontsize=13)

    for i, ax in enumerate(axes):
        ax.legend(handles=stra_legends)
        ax.set_ylabel('')
        ax.set_yticks(np.arange(6))
        ax.set_yticklabels(['HiFi-10kb', 'HiFi-15kb', 'HiFi-18kb', 'ONT-9kb', 'ONT-19kb', 'ONT-30kb'], fontsize=13)

        if i == 0:
            ax.set_xlim(0.7, 1)
            ax.set_xticks(np.linspace(0.7, 1, 4))
            ax.set_xticklabels([int(ele * 100) for ele in np.linspace(0.7, 1, 4)], fontsize=12)
        else:
            ax.set_xlim(0.4, 1)
            ax.set_xticks(np.linspace(0.4, 1, 4))
            ax.set_xticklabels([int(ele * 100) for ele in np.linspace(0.4, 1, 4)], fontsize=12)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()

    fig.savefig(f'{figdir}/fig5a.pdf')

    df_cmrg = pd.read_csv(f'{WORKDIR}/CMRGs/assembly_read_cmrgs_truvari.tsv', sep='\t', header=[0])

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 3))
    sns.boxplot(data=df_cmrg, x='recall', y='dataset', hue='strategy', hue_order=stra_order, palette=stra_color,
                ax=axes[0])
    axes[0].set_xlabel('Recall', fontsize=13)
    sns.boxplot(data=df_cmrg, x='precision', y='dataset', hue='strategy', hue_order=stra_order, palette=stra_color,
                ax=axes[1])
    axes[1].set_xlabel('Precision', fontsize=13)

    for i, ax in enumerate(axes):
        ax.set_ylabel('')
        ax.set_yticks(np.arange(6))
        ax.set_yticklabels(['HiFi-10kb', 'HiFi-15kb', 'HiFi-18kb', 'ONT-9kb', 'ONT-19kb', 'ONT-30kb'], fontsize=13)

        if i == 0:
            ax.set_xlim(0.7, 1)
            ax.set_xticks(np.linspace(0.7, 1, 4))
            ax.set_xticklabels([int(ele * 100) for ele in np.linspace(0.7, 1, 4)], fontsize=12)
        else:
            ax.set_xlim(0.4, 1)
            ax.set_xticks(np.linspace(0.4, 1, 4))
            ax.set_xticklabels([int(ele * 100) for ele in np.linspace(0.4, 1, 4)], fontsize=12)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.legend(handles=stra_legends)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig5b.pdf')
    # plt.show()


def plot_5c5d(figdir, merge_type):

    plats = {'HiFi': ('HiFi-10kb', 'HiFi-15kb', 'HiFi-18kb'), 'ONT':('ONT-9kb', 'ONT-19kb', 'ONT-30kb')}
    stras = ['read', 'assm']

    supp_color = {'HiFi':{'1': '#fbd6b6', '2': '#f8b67f', '3': '#F2984E'},
                  'ONT':{'1': '#d1ddd2', '2': '#afc4b1', '3': '#8EAB8F'}}

    supp_legend = [Patch(label='HiFi', color='#F2984E'), Patch(label='ONT', color='#8EAB8F')]

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    plat_idx = 0
    barwidth = 0.3
    r1 = np.arange(len(stras))
    r2 = [r + barwidth + 0.02 for r in r1]
    rs = [r1, r2]

    for plat, set_label in plats.items():
        supp_counts = [[0, 0], [0, 0], [0, 0]]
        for col_idx, stra in enumerate(stras):


            merged = f'{WORKDIR}/CMRGs/plat_compare/{stra}.{plat}.{merge_type}.merged.vcf'
            supp, merged_total = get_survivor_supp(merged)

            supp_counts[0][col_idx] += supp[3]
            supp_counts[1][col_idx] += supp[2]
            supp_counts[2][col_idx] += supp[1]

        total = [i + j + k for i, j , k in zip(supp_counts[0], supp_counts[1], supp_counts[2])]
        supp3_pcrt = [i / j * 100 for i, j in zip(supp_counts[0], total)]
        supp2_pcrt = [i / j * 100 for i, j in zip(supp_counts[1], total)]
        supp1_pcrt = [i / j * 100 for i, j in zip(supp_counts[2], total)]

        ax.bar(rs[plat_idx], supp3_pcrt, width=barwidth, color=supp_color[plat]['3'])
        ax.bar(rs[plat_idx], supp2_pcrt, bottom=supp3_pcrt, width=barwidth, color=supp_color[plat]['2'])
        bottom = [i + j for i, j in zip(supp3_pcrt, supp2_pcrt)]

        ax.bar(rs[plat_idx], supp1_pcrt, bottom=bottom, width=barwidth, color=supp_color[plat]['1'])

        plat_idx += 1

    ax.set_ylabel(f'% of {merge_type}', fontsize=13)

    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks([r + (barwidth + 0.02) / 2 for r in r1])
    ax.set_xticklabels(['Read', 'Assembly'], fontsize=13)

    ax.legend(handles=supp_legend)

    fig.tight_layout()

    fig.savefig(f'{figdir}/fig5c-{merge_type}.pdf')

    for plat in plats:
        fig, ax = plt.subplots(1, 1, figsize=(4, 3))
        ax.set_title(plat, fontsize=13)
        merged = f'{WORKDIR}/CMRGs/stra_compare/stra.{plat}.{merge_type}s.merged.vcf'
        supp_vec, merged_total = get_survivor_suppvec(merged)
        n_10, n_01, n_11 = 0, 0, 0
        if '10' in supp_vec:
            n_10 = supp_vec['10']

        if '01' in supp_vec:
            n_01 = supp_vec['01']

        if '11' in supp_vec:
            n_11 = supp_vec['11']

        venn2(subsets={'10': n_10, '01': n_01, '11': n_11}, set_labels=('Read', 'Assembly'), ax=ax)
        fig.tight_layout()
        fig.savefig(f'{figdir}/fig5d-{plat}-{merge_type}s.venn.pdf')



def plot_5e(figdir):

    df_diffcov = pd.read_csv(f'{WORKDIR}/CMRGs/stra_diffcov_truvari.tsv', sep='\t', header=[0])

    plat_stra_rp = {'HiFi': {'Assembly': {'35X':[[], []], '5X':[[], []]}, 'Read': {'35X':[[], []], '5X':[[], []]}},
               'ONT': {'Assembly': {'35X':[[], []], '5X':[[], []]}, 'Read': {'35X':[[], []], '5X':[[], []]}}}

    plat_caller_rp = {'HiFi': {TOOLMAP[caller]: [0, 0] for caller in READCALLERS},
                      'ONT': {TOOLMAP[caller]: [0, 0] for caller in READCALLERS}}

    read_based_callers = [TOOLMAP[ele] for ele in READCALLERS]
    for idx, row in df_diffcov.iterrows():
        cov, stra, recall, precision, plat, caller = row['coverage'], row['stra'], float(row['recall']), float(row['precision']), row['plat'], row['caller']
        if cov in ['35X', '5X']:
            plat_stra_rp[plat][stra][cov][0].append(recall)
            plat_stra_rp[plat][stra][cov][1].append(precision)

        if cov == '5X' and caller in read_based_callers:
            plat_caller_rp[plat][caller][0] += recall
            plat_caller_rp[plat][caller][1] += precision


    stra_legends_line = [Line2D([0], [0], label='Assembly', lw=2, color=STRACOLORS['Assembly']),
                        Line2D([0], [0], label='Read', lw=2, color=STRACOLORS['Read']),
                         Line2D([0], [0], label='Precision', lw=2, ls='--', color='black'),
                         Line2D([0], [0], label='Recall', lw=2, color='black')]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))
    df_hifi_diffcov = df_diffcov[df_diffcov['plat'] == 'HiFi']
    df_ont_diffcov = df_diffcov[df_diffcov['plat'] == 'ONT']

    sns.lineplot(data=df_hifi_diffcov, x='coverage', y='recall',
                hue='stra', hue_order=['Assembly', 'Read'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']], lw=2.5, marker='o',
                 markersize=8, ci=None, ax=axes[0])

    sns.lineplot(data=df_hifi_diffcov, x='coverage', y='precision',
                 hue='stra', hue_order=['Assembly', 'Read'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']],
                 lw=2.5, marker='o', ls='--', markersize=8, ci=None, ax=axes[0])

    axes[0].set_title('HiFi-18kb', fontsize=13)

    sns.lineplot(data=df_ont_diffcov, x='coverage', y='recall',
                hue='stra', hue_order=['Assembly', 'Read'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']],
                 ci=None, lw=2.5, marker='o', markersize=8, ax=axes[1])

    sns.lineplot(data=df_ont_diffcov, x='coverage', y='precision',
                 hue='stra', hue_order=['Assembly', 'Read'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']],
                 ci=None, lw=2.5, ls='--', marker='o', markersize=8, ax=axes[1])

    axes[1].set_title('ONT-30kb', fontsize=13)

    for i, ax in enumerate(axes):

        ax.legend(handles=stra_legends_line)

        ax.set_ylabel('')

        ax.set_ylim(0, 1)
        ax.set_yticks(np.linspace(0, 1, 5))
        ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)], fontsize=12)

        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['35X', '20X', '10X', '5X'], fontsize=13)
        ax.set_xlabel('')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig5e.pdf')
    # plt.show()


def plot_5f(figdir):

    df_diffcov = pd.read_csv(f'{WORKDIR}/CMRGs/stra_diffcov_truvari.tsv', sep='\t', header=[0])

    plat_stra_rp = {'HiFi': {'Assembly': {'35X':[[], []], '5X':[[], []]}, 'Read': {'35X':[[], []], '5X':[[], []]}},
               'ONT': {'Assembly': {'35X':[[], []], '5X':[[], []]}, 'Read': {'35X':[[], []], '5X':[[], []]}}}

    plat_caller_rp = {'HiFi': {TOOLMAP[caller]: [0, 0] for caller in READCALLERS},
                      'ONT': {TOOLMAP[caller]: [0, 0] for caller in READCALLERS}}

    read_based_callers = [TOOLMAP[ele] for ele in READCALLERS]
    for idx, row in df_diffcov.iterrows():
        cov, stra, recall, precision, plat, caller = row['coverage'], row['stra'], float(row['recall']), float(row['precision']), row['plat'], row['caller']
        if cov in ['35X', '5X']:
            plat_stra_rp[plat][stra][cov][0].append(recall)
            plat_stra_rp[plat][stra][cov][1].append(precision)

        if cov == '5X' and caller in read_based_callers:
            plat_caller_rp[plat][caller][0] += recall
            plat_caller_rp[plat][caller][1] += precision

    caller_legends = [Line2D([0], [0], color='white', mec='black', marker=TOOLMARKERS[ele], label=TOOLMAP[ele], ms=10) for ele in READCALLERS]

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]
    caller_legends.extend(plat_legends)

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    for plat, call_rp in plat_caller_rp.items():
        for caller, (recall, precision) in call_rp.items():
            ax.scatter(recall, precision, color=PLATCOLORS[plat], marker=TOOLMARKERS2[caller], s=150)

    ax.legend(handles=caller_legends)

    ax.set_xlabel('Recall (%)', fontsize=13)
    ax.set_xlim(0.6, 1)
    ax.set_xticks(np.linspace(0.6, 1, 5))
    ax.set_xticklabels([int(val * 100) for val in np.linspace(0.6, 1, 5)], fontsize=12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylabel('Precision (%)', fontsize=13)
    ax.set_ylim(0.2, 1)
    ax.set_yticks(np.linspace(0.2, 1, 5))
    ax.set_yticklabels([int(val * 100) for val in np.linspace(0.2, 1, 5)], fontsize=12)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig5f.pdf')
    # plt.show()

def main():
    print('\n===== Creating Figure 5 =====')

    if not os.path.exists(f'{FIGDIR}/Fig5'):
        os.mkdir(f'{FIGDIR}/Fig5')

    plot_5a5b(f'{FIGDIR}/Fig5')
    print(f'Figures saved to {FIGDIR}/Fig5/fig5a.pdf; {FIGDIR}/Fig5/fig5b.pdf')

    plot_5c5d(f'{FIGDIR}/Fig5', 'fp')
    print(f'Figures saved to {FIGDIR}/Fig5/fig5c-fp.pdf; {FIGDIR}/Fig5/fig5b-HiFi-fps.venn.pdf; {FIGDIR}/Fig5/fig5b-ONT-fps.venn.pdf')

    plot_5c5d(f'{FIGDIR}/Fig5', 'fn')
    print(f'Figures saved to {FIGDIR}/Fig5/fig5c-fn.pdf; {FIGDIR}/Fig5/fig5b-HiFi-fns.venn.pdf; {FIGDIR}/Fig5/fig5b-ONT-fns.venn.pdf')

    plot_5e(f'{FIGDIR}/Fig5')
    print(f'Figures saved to {FIGDIR}/Fig5/fig5e.pdf')

    plot_5f(f'{FIGDIR}/Fig5')
    print(f'Figures saved to {FIGDIR}/Fig5/fig5f.pdf')

if __name__ == '__main__':
    main()
