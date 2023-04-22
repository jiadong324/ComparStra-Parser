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
from matplotlib.patches import Patch

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


def plot_s6a(figdir):
    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'], 'ONT': ['ont_9kb', 'ont_19kb' , 'ont_30kb']}

    datasets_labels = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb', 'ont_9kb', 'ont_19kb', 'ont_30kb']
    hq_svnums = []

    for region_type in ['WGS', 'ExTD']:
        for plat, datasets in datasets_dict.items():
            assembly_unique = {'del': [], 'ins': []}
            overlaps = {'del': [], 'ins': []}
            read_unique = {'del': [], 'ins': []}

            read_pcrt = {'del': [], 'ins': []}
            assembly_pcrt = {'del': [], 'ins': []}

            for dataset in datasets:
                for svtype in ['ins', 'del']:
                    suppvec_dict = {}
                    for line in open(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/stra_hq_compare/{region_type}/{dataset}.{svtype}.suppvec_info.{region_type}.tsv'):
                        entries = line.strip().split('\t')
                        suppvec, count = entries[0], int(entries[1])
                        suppvec_dict[suppvec] = count
                        if suppvec == '10':
                            assembly_unique[svtype].append(count)
                        elif suppvec == '11':
                            overlaps[svtype].append(count)
                        else:
                            read_unique[svtype].append(count)

                    hq_svnums.append((region_type, SVTYPEMAP[svtype], 'Read', suppvec_dict['01'] + suppvec_dict['11'], dataset, f'{plat}-{SVTYPEMAP[svtype]}'))
                    hq_svnums.append((region_type, SVTYPEMAP[svtype], 'Assembly', suppvec_dict['10'] + suppvec_dict['11'], dataset, f'{plat}-{SVTYPEMAP[svtype]}'))

                    read_pcrt[svtype].append(suppvec_dict['11'] / (suppvec_dict['01'] + suppvec_dict['11']))
                    assembly_pcrt[svtype].append(suppvec_dict['11'] / (suppvec_dict['10'] + suppvec_dict['11']))

    df_svnum = pd.DataFrame(hq_svnums, columns=['region', 'svtype', 'stra', 'count', 'dataset', 'plat'])

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(7, 4))

    sns.lineplot(data=df_svnum[df_svnum['svtype'] == 'INS'], x='dataset', y='count', hue='stra', style='region', lw=2,
                palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], marker='o', ax=axes[0])

    sns.lineplot(data=df_svnum[df_svnum['svtype'] == 'DEL'], x='dataset', y='count', hue='stra', style='region', lw=2,
                palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], marker='o', ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('# of SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.legend(title='')

        ax.set_yticks(np.linspace(0, 12000, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0,12000, 5)], fontsize=12)

        ax.set_xlabel('')
        ax.set_xticks(np.arange(len(datasets_labels)))
        ax.set_xticklabels([PLATMAP[ele] for ele in datasets_labels], rotation=90, fontsize=13)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs6a.pdf')
    # plt.show()

def plot_s6b(figdir):
    # datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'], 'ONT': ['ont_9kb', 'ont_19kb' , 'ont_30kb']}

    hq_svnums = []

    for region_type in ['WGS', 'ExTD']:
        for plat, datasets in DATASET_DICT.items():
            assembly_unique = {'del': [], 'ins': []}
            overlaps = {'del': [], 'ins': []}
            read_unique = {'del': [], 'ins': []}

            read_pcrt = {'del': [], 'ins': []}
            assembly_pcrt = {'del': [], 'ins': []}

            for dataset in datasets:
                for svtype in ['ins', 'del']:
                    suppvec_dict = {}
                    for line in open(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/stra_hq_compare/{region_type}/{dataset}.{svtype}.suppvec_info.{region_type}.tsv'):
                        entries = line.strip().split('\t')
                        suppvec, count = entries[0], int(entries[1])
                        suppvec_dict[suppvec] = count
                        if suppvec == '10':
                            assembly_unique[svtype].append(count)
                        elif suppvec == '11':
                            overlaps[svtype].append(count)
                        else:
                            read_unique[svtype].append(count)

                    hq_svnums.append((region_type, SVTYPEMAP[svtype], 'Read', suppvec_dict['01'] + suppvec_dict['11'], dataset, f'{plat}-{SVTYPEMAP[svtype]}'))
                    hq_svnums.append((region_type, SVTYPEMAP[svtype], 'Assembly', suppvec_dict['10'] + suppvec_dict['11'], dataset, f'{plat}-{SVTYPEMAP[svtype]}'))

                    read_pcrt[svtype].append(suppvec_dict['11'] / (suppvec_dict['01'] + suppvec_dict['11']))
                    assembly_pcrt[svtype].append(suppvec_dict['11'] / (suppvec_dict['10'] + suppvec_dict['11']))

            barwidth = 0.3

            r1 = np.arange(len(datasets))
            r2 = [r + barwidth + 0.04 for r in r1]
            rs = [r1, r2]
            xticks = [r + (barwidth + 0.04) / 2 for r in r1]

            legends = [Patch(label='Assembly specific', edgecolor='k', facecolor=STRACOLORS['Assembly']),
                       Patch(label='Read specific', edgecolor='k', facecolor=STRACOLORS['Read']),
                       Patch(label='Shared', facecolor='white', edgecolor='k', hatch='///')]

            fig, ax = plt.subplots(1, 1, figsize=(5, 4))
            for idx, svtype in enumerate(['del', 'ins']):
                ax.bar(rs[idx], assembly_unique[svtype], width=barwidth, color=STRACOLORS['Assembly'], edgecolor='black')
                ax.bar(rs[idx], overlaps[svtype], width=barwidth, color='white', hatch='//', edgecolor='black',
                       bottom=assembly_unique[svtype])
                ax.bar(rs[idx], read_unique[svtype], width=barwidth, color=STRACOLORS['Read'], edgecolor='black',
                       bottom=[i + j for i, j in zip(assembly_unique[svtype], overlaps[svtype])])

            # ax.set_yticks(np.linspace(0, 6000, 5))
            # ax.set_yticklabels([int(val) for val in np.linspace(0, 6000, 5)], fontsize=12)
            ax.set_ylabel(f'# of {region_type}-SVs', fontsize=13)

            ax.set_xticks(xticks)
            ax.set_xticklabels([PLATMAP[ele] for ele in datasets], fontsize=13, rotation=90)

            ax.legend(handles=legends)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            fig.tight_layout()
            fig.savefig(f'{figdir}/figs6b-{region_type}-{plat}.pdf')

    # plt.show()

def plot_s6c(figdir):

    datasets_dict = {'HiFi': 'hifi_18kb', 'ONT': 'ont_30kb'}

    invalid_counter = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}
    valid_total = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}
    uniques_total = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}

    samples = ['HG003', 'HG004']

    for plat, dataset in datasets_dict.items():
        fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))
        for svtype_idx, svtype in enumerate(['ins', 'del']):
            trio_scores = []
            trio_valid_qs = {}
            total_valid_svs = 0
            for sample_idx, sample in enumerate(samples):
                with open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/vapor_{svtype}_{sample}.tsv', 'r') as f:
                    next(f)
                    for line in f:
                        entries = line.strip().split('\t')
                        uniques_total[plat][svtype] += 1
                        if entries[6] == 'NA':
                            continue
                        sv_id, vapor_gs = entries[4], float(entries[6])
                        trio_scores.append((vapor_gs, sample, plat))

                        if sv_id in trio_valid_qs:
                            trio_valid_qs[sv_id].append(vapor_gs)
                        else:
                            trio_valid_qs[sv_id] = []
                            trio_valid_qs[sv_id].append(vapor_gs)

            invalid_count = 0
            for sv_id, scores in trio_valid_qs.items():
                if len(scores) == 1:
                    continue
                total_valid_svs += 1
                if scores[0] == 0 and scores[1] == 0:
                    invalid_count += 1

            error_rate = round(invalid_count / total_valid_svs * 100, 2)

            invalid_counter[plat][svtype] = invalid_count
            valid_total[plat][svtype] = total_valid_svs

            # print(f'{plat}: {svtype} error rate: {error_rate}. Valid total SVs: {total_valid_svs}')

            df_scores = pd.DataFrame(trio_scores, columns=['score', 'sample', 'plat'])

            sns.histplot(data=df_scores, x='score', hue='sample', ax=axes[svtype_idx])

            axes[svtype_idx].set_xlabel('VaPoR_GS', fontsize=13)
            axes[svtype_idx].set_ylim(0, 85)
            axes[svtype_idx].set_yticks(np.linspace(0, 80, 5))
            axes[svtype_idx].set_yticklabels([int(val) for val in np.linspace(0, 80, 5)], fontsize=12)

            axes[svtype_idx].set_ylabel(f'# of read specific SVs', fontsize=13)

            axes[svtype_idx].spines['top'].set_visible(False)
            axes[svtype_idx].spines['right'].set_visible(False)
            axes[svtype_idx].text(.25, .6, f'{SVTYPEMAP[svtype]} FDR={error_rate}%', fontsize=12, transform=axes[svtype_idx].transAxes)

        fig.tight_layout()
        fig.savefig(f'{figdir}/figs6c-{plat}.pdf')

    # plt.show()

def main():
    print('\n==== Creating Extended Data Fig 6 =====')

    if not os.path.exists(f'{FIGDIR}/FigS6'):
        os.mkdir(f'{FIGDIR}/FigS6')

    print(f'All panels are under: {FIGDIR}/FigS6')

    plot_s6a(f'{FIGDIR}/FigS6')
    plot_s6b(f'{FIGDIR}/FigS6')
    plot_s6c(f'{FIGDIR}/FigS6')

if __name__ == '__main__':
    main()
