#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/16

'''
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Constant import *
from Reader.Read import process_read_calls
from Reader.Assembly import process_assembly_calls

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

def plot_2a(figdir):

    ins_pcrt = []

    sv_count_list = []
    align_sv_count = pd.read_csv(f'{WORKDIR}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in align_sv_count.iterrows():
        plat = 'HiFi'
        all_sv_num, ins_num, del_num, aligner, dataset, caller = int(row['all_num']), int(row['ins_num']), int(row['del_num']), row['aligner'], row['dataset'], row['caller']
        if 'ont' in dataset:
            plat = 'ONT'

        ins_pcrt.append((TOOLMAP[caller], aligner, PLATMAP[dataset], ins_num * 100 / all_sv_num, plat, 'Read-based'))
        # del_pcrt.append((TOOLMAP[caller], aligner, PLATMAP[dataset], del_num * 100 / all_sv_num, plat, 'Read-based'))
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], all_sv_num, plat, 'Read-based'))

    pav_sv_count = pd.read_csv(f'{WORKDIR}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        svcount = ins_num + del_num + int(row['inv_num'])
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))

        ins_pcrt.append((TOOLMAP[caller], assembler, PLATMAP[dataset], ins_num * 100 / svcount, plat, 'Assembly-based'))
        # del_pcrt.append((TOOLMAP[caller], assembler, PLATMAP[dataset], del_num * 100 / svcount, plat, 'Assembly-based'))

    svimasm_sv_count = pd.read_csv(f'{WORKDIR}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))

        ins_pcrt.append((TOOLMAP[caller], assembler, PLATMAP[dataset], ins_num * 100 / svcount, plat, 'Assembly-based'))

    df_ins_pcrt = pd.DataFrame(ins_pcrt, columns=['caller', 'aa', 'dataset', 'pcrt', 'plat', 'stra'])
    df_sv_counts = pd.DataFrame(sv_count_list, columns=['caller', 'dataset', 'count', 'plat', 'stra'])

    hue_order = ['Assembly-based', 'Read-based']
    hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]
    stra_legend = [Patch(label='Assembly', color=STRACOLORS['Assembly']), Patch(label='Read', color=STRACOLORS['Read'])]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))
    sns.boxplot(data=df_sv_counts, y='dataset', x='count', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[0])
    sns.boxplot(data=df_ins_pcrt, y='dataset', x='pcrt', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[1])
    for i, ax in enumerate(axes):

        ax.set_ylabel('')
        ax.set_yticks([0, 1, 2, 3, 4, 5])
        ax.set_yticklabels(['HiFi-10kb', 'HiFi-15kb', 'HiFi-18kb', 'ONT-9kb', 'ONT-19kb', 'ONT-30kb'], fontsize=13)
        if i == 0:
            ax.set_xlim(10000, 40000)
            ax.set_xticks(np.linspace(10000, 40000, 4))
            ax.set_xticklabels([int(val) for val in np.linspace(10, 40, 4)], fontsize=12)

            ax.set_xlabel('# of SVs (x$10^3$)', fontsize=13)
        else:
            ax.set_xlim(0, 100)
            ax.set_xticks(np.linspace(0, 100, 5))
            ax.set_xticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

            ax.set_xlabel('% of INS', fontsize=13)

        ax.legend(handles=stra_legend)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='x', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig2a.pdf')
    # plt.show()




def main():

    print('\n===== Prepare data for Figure 2d and 2e =====')
    print(f'Intermediate files include: '
          f'\n{WORKDIR}/caller_sv_counts_region.tsv'
          f'\n{WORKDIR}/caller_sv_count.tsv'
          f'\n{WORKDIR}/pav_sv_counts_region.tsv'
          f'\n{WORKDIR}/pav_sv_counts.tsv'
          f'\n{WORKDIR}/svimasm_sv_counts_region.tsv'
          f'\n{WORKDIR}/svimasm_sv_counts.tsv\n')

    process_read_calls()
    process_assembly_calls()

    if not os.path.exists(FIGDIR):
        os.mkdir(FIGDIR)

    if not os.path.exists(f'{FIGDIR}/Fig2'):
        os.mkdir(f'{FIGDIR}/Fig2')

    plot_2a(f'{FIGDIR}/Fig2')
    print(f'Figures saved to {FIGDIR}/Fig2/fig2a.pdf')

if __name__ == '__main__':
    main()