#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/2/21

'''
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib_venn import venn2, venn3

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Functions import *
from Helpers.Constant import *

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

def plot_s7a7b(figdir):
    plats = {'HiFi': ('HiFi-10kb', 'HiFi-15kb', 'HiFi-18kb'), 'ONT':('ONT-9kb', 'ONT-19kb', 'ONT-30kb')}
    stras = ['read', 'assm']

    stra_dict = {'read': 'Read', 'assm': 'Assembly'}
    for merge_type in ['fp', 'fn']:
        for plat, set_label in plats.items():
            for col_idx, stra in enumerate(stras):

                merged = f'{WORKDIR}/CMRGs/plat_compare/{stra}.{plat}.{merge_type}.merged.vcf'
                supp_vec, merged_total = get_survivor_suppvec(merged)

                supp1, supp2, supp3 = 0, 0, 0
                n_100, n_010, n_001, n_110, n_101, n_011, n_111 = 0, 0, 0, 0, 0, 0, 0
                if '100' in supp_vec:
                    n_100 = supp_vec['100']
                    supp1 += supp_vec['100']
                if '010' in supp_vec:
                    n_010 = supp_vec['010']
                    supp1 += supp_vec['010']
                if '001' in supp_vec:
                    n_001 = supp_vec['001']
                    supp1 += supp_vec['001']

                if '110' in supp_vec:
                    n_110 = supp_vec['110']
                    supp2 += supp_vec['110']
                if '101' in supp_vec:
                    n_101 = supp_vec['101']
                    supp2 += supp_vec['101']
                if '011' in supp_vec:
                    n_011 = supp_vec['011']
                    supp2 += supp_vec['011']

                if '111' in supp_vec:
                    n_111 = supp_vec['111']
                    supp3 += supp_vec['111']

                fig, ax = plt.subplots(1, 1, figsize=(4, 3))
                ax.set_title(f'{stra_dict[stra]}-{merge_type}', fontsize=13)
                venn3(subsets={'100': n_100, '010': n_010, '001': n_001, '111': n_111, '110': n_110, '101': n_101, '011': n_011},
                      set_labels=set_label, ax=ax)
                fig.tight_layout()
                fig.savefig(f'{figdir}/figs7ab-{stra}-{merge_type}-{plat}.pdf')

    # plt.show()


def plot_s7c(figdir):

    df_diffcov = pd.read_csv(f'{WORKDIR}/CMRGs/stra_diffcov_truvari.tsv', sep='\t', header=[0])

    df_assm_diffcov = df_diffcov[df_diffcov['stra'] == 'Assembly']

    legends_lines = [Line2D([0], [0], label='shasta', lw=2, color=ASSMBLERCOLOR['shasta']),
                         Line2D([0], [0], label='flye', lw=2, color=ASSMBLERCOLOR['flye']),
                         Line2D([0], [0], label='hifiasm', lw=2, color=ASSMBLERCOLOR['hifiasm']),
                         Line2D([0], [0], label='SVIM-asm', lw=2, ls='--', color='black'),
                     Line2D([0], [0], label='PAV', lw=2, color='black')]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))

    sns.lineplot(data=df_assm_diffcov[df_assm_diffcov['plat'] == 'HiFi'], x='coverage', y='recall', hue='aa',
                 hue_order=['hifiasm', 'flye'], palette=[ASSMBLERCOLOR['hifiasm'], ASSMBLERCOLOR['flye']], style='caller',
                 lw=2, marker='o', markersize=8, ax=axes[0])
    axes[0].set_title('HiFi', fontsize=13)

    sns.lineplot(data=df_assm_diffcov[df_assm_diffcov['plat'] == 'ONT'], x='coverage', y='recall', hue='aa',
                 hue_order=['shasta', 'flye'], palette=[ASSMBLERCOLOR['shasta'], ASSMBLERCOLOR['flye']], style='caller',
                 lw=2, marker='o', markersize=8, ax=axes[1])
    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Recall', fontsize=13)
            ax.legend(handles=legends_lines)
        else:
            ax.legend('', frameon=False)
            ax.set_ylabel('')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        ax.set_xlabel('')
        ax.set_ylim(0, 1)
        ax.set_yticks(np.linspace(0, 1, 5))
        ax.set_yticklabels([int(val * 100) for val in np.linspace(0, 1, 5)], fontsize=12)

        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['35X', '20X', '10X', '5X'], fontsize=13)


    fig.tight_layout()
    fig.savefig(f'{figdir}/figs7c-1.pdf')

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))

    sns.lineplot(data=df_assm_diffcov[df_assm_diffcov['plat'] == 'HiFi'], x='coverage', y='precision', hue='aa',
                 hue_order=['hifiasm', 'flye'], palette=[ASSMBLERCOLOR['hifiasm'], ASSMBLERCOLOR['flye']],
                 style='caller', lw=2, marker='o', markersize=8, ax=axes[0])
    axes[0].set_title('HiFi', fontsize=13)

    sns.lineplot(data=df_assm_diffcov[df_assm_diffcov['plat'] == 'ONT'], x='coverage', y='precision', hue='aa',
                 hue_order=['shasta', 'flye'], palette=[ASSMBLERCOLOR['shasta'], ASSMBLERCOLOR['flye']], style='caller',
                 lw=2, marker='o', markersize=8, ax=axes[1])

    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Precision (%)', fontsize=13)
            ax.legend(handles=legends_lines)
        else:
            ax.legend('', frameon=False)
            ax.set_ylabel('')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        ax.set_xlabel('')
        ax.set_ylim(0, 1)
        ax.set_yticks(np.linspace(0, 1, 5))
        ax.set_yticklabels([int(val * 100) for val in np.linspace(0, 1, 5)], fontsize=12)

        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['35X', '20X', '10X', '5X'], fontsize=13)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs7c-2.pdf')

    plat_stra_rp = {
        'HiFi': {'Assembly': {'35X': [[], []], '5X': [[], []]}, 'Read': {'35X': [[], []], '5X': [[], []]}},
        'ONT': {'Assembly': {'35X': [[], []], '5X': [[], []]}, 'Read': {'35X': [[], []], '5X': [[], []]}}}

    plat_caller_rp = {'HiFi': {TOOLMAP[caller]: [0, 0] for caller in READCALLERS},
                      'ONT': {TOOLMAP[caller]: [0, 0] for caller in READCALLERS}}

    read_based_callers = [TOOLMAP[ele] for ele in READCALLERS]
    for idx, row in df_diffcov.iterrows():
        cov, stra, recall, precision, plat, caller = row['coverage'], row['stra'], float(row['recall']), float(
            row['precision']), row['plat'], row['caller']
        if cov in ['35X', '5X']:
            plat_stra_rp[plat][stra][cov][0].append(recall)
            plat_stra_rp[plat][stra][cov][1].append(precision)

        if cov == '5X' and caller in read_based_callers:
            plat_caller_rp[plat][caller][0] += recall
            plat_caller_rp[plat][caller][1] += precision

    fig, axes = plt.subplots(1, 2, figsize=(6, 4))
    df_read_diffcov = df_diffcov[df_diffcov['stra'] == 'Read']
    sns.lineplot(data=df_read_diffcov[df_read_diffcov['plat'] == 'ONT'], x='coverage', y='recall', hue='caller',
                 hue_order=[TOOLMAP[ele] for ele in READCALLERS],
                 palette=[TOOLCOLORS2[TOOLMAP[ele]] for ele in READCALLERS],
                 lw=2, marker='o', markersize=8, ax=axes[0])

    axes[0].set_ylabel('Recall of ONT (%)', fontsize=13)

    sns.lineplot(data=df_read_diffcov[df_read_diffcov['plat'] == 'ONT'], x='coverage', y='precision', hue='caller',
                 hue_order=[TOOLMAP[ele] for ele in READCALLERS],
                 palette=[TOOLCOLORS2[TOOLMAP[ele]] for ele in READCALLERS],
                 lw=2, marker='o', markersize=8, ax=axes[1])

    axes[1].set_ylabel('Precision of ONT (%)', fontsize=13)

    for i, ax in enumerate(axes):
        ax.set_ylim(0.4, 1)
        ax.set_yticks(np.linspace(0.4, 1, 4))
        ax.set_yticklabels([int(val * 100) for val in np.linspace(0.4, 1, 4)], fontsize=12)

        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['35X', '20X', '10X', '5X'], rotation=60, fontsize=13)
        ax.set_xlabel('')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs7c-3.pdf')

    fig, axes = plt.subplots(1, 2, figsize=(6, 4))
    df_read_diffcov = df_diffcov[df_diffcov['stra'] == 'Read']
    sns.lineplot(data=df_read_diffcov[df_read_diffcov['plat'] == 'HiFi'], x='coverage', y='recall', hue='caller',
                 hue_order=[TOOLMAP[ele] for ele in READCALLERS],
                 palette=[TOOLCOLORS2[TOOLMAP[ele]] for ele in READCALLERS],
                 lw=2, marker='o', markersize=8, ax=axes[0])

    axes[0].set_ylabel('Recall of HiFi (%)', fontsize=13)

    sns.lineplot(data=df_read_diffcov[df_read_diffcov['plat'] == 'HiFi'], x='coverage', y='precision', hue='caller',
                 hue_order=[TOOLMAP[ele] for ele in READCALLERS],
                 palette=[TOOLCOLORS2[TOOLMAP[ele]] for ele in READCALLERS],
                 lw=2, marker='o', markersize=8, ax=axes[1])

    axes[1].set_ylabel('Precision of HiFi (%)', fontsize=13)

    for i, ax in enumerate(axes):
        ax.set_ylim(0.4, 1)
        ax.set_yticks(np.linspace(0.4, 1, 4))
        ax.set_yticklabels([int(val * 100) for val in np.linspace(0.4, 1, 4)], fontsize=12)

        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['35X', '20X', '10X', '5X'], rotation=60, fontsize=13)
        ax.set_xlabel('')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs7c-4.pdf')
    # plt.show()


def main():
    print('\n===== Creating Extended Data Fig 7 =====')

    if not os.path.exists(f'{FIGDIR}/FigS7'):
        os.mkdir(f'{FIGDIR}/FigS7')

    plot_s7a7b(f'{FIGDIR}/FigS7')
    print(f'Figures saved to {FIGDIR}/Fig5/figs7ab-*-*-*.pdf')

    plot_s7c(f'{FIGDIR}/FigS7')
    print(f'Figures saved to {FIGDIR}/Fig5/figs7c-1.pdf; {FIGDIR}/Fig5/figs7c-2.pdf; {FIGDIR}/Fig5/figs7c-3.pdf; {FIGDIR}/Fig5/figs7c-4.pdf')

if __name__ == '__main__':
    main()