#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/17

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

def plot_s1a(figdir):
    ins_pcrt = []
    sv_count_list = []
    align_sv_count = pd.read_csv(f'{WORKDIR}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in align_sv_count.iterrows():
        plat = 'HiFi'
        all_sv_num, ins_num, del_num, aligner, dataset, caller = int(row['all_num']), int(row['ins_num']), int(
            row['del_num']), row['aligner'], row['dataset'], row['caller']
        if 'ont' in dataset:
            plat = 'ONT'

        ins_pcrt.append((TOOLMAP[caller], aligner, PLATMAP[dataset], ins_num * 100 / all_sv_num, plat, 'Read-based'))
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

    svimasm_sv_count = pd.read_csv(f'{WORKDIR}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))

        ins_pcrt.append((TOOLMAP[caller], assembler, PLATMAP[dataset], ins_num * 100 / svcount, plat, 'Assembly-based'))

    df_ins_pcrt = pd.DataFrame(ins_pcrt, columns=['caller', 'aa', 'dataset', 'pcrt', 'plat', 'stra'])
    assemblers = ['hifiasm', 'flye', 'shasta']

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    read_ins_pcrt = df_ins_pcrt[df_ins_pcrt['stra'] == 'Read-based']
    xticks = np.arange(len(DATASETS))

    for caller in READCALLERS:
        sns.stripplot(data=read_ins_pcrt[read_ins_pcrt['caller'] == TOOLMAP[caller]], x='dataset', y='pcrt',
                      marker=TOOLMARKERS[caller], size=8, edgecolor='black', linewidth=1, color='white', ax=axes[0])

    sns.stripplot(data=read_ins_pcrt, x='dataset', y='pcrt', hue='aa', hue_order=READALIGNERS,
                  palette=[ALIGNERCOLOR[ele] for ele in READALIGNERS], size=7, edgecolor='black', linewidth=1, ax=axes[1])

    for i, ax in enumerate(axes):
        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in DATASETS], rotation=60, fontsize=13)
        ax.set_xlabel('')

        ax.set_ylim(20, 80)
        ax.set_yticks(np.linspace(20, 80, 4))
        ax.set_yticklabels([int(val) for val in np.linspace(20, 80, 4)], fontsize=12)
        if i == 0:
            ax.set_ylabel('% of INS', fontsize=13)
            ax.legend('', frameon=False)
        else:
            ax.set_ylabel('')
            ax.legend(title='', loc='best')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs1a-1.pdf')

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    assm_ins_pcrt = df_ins_pcrt[df_ins_pcrt['stra'] == 'Assembly-based']
    for caller in ASMCALLERS:
        sns.stripplot(data=assm_ins_pcrt[assm_ins_pcrt['caller'] == TOOLMAP[caller]], x='dataset', y='pcrt', marker=TOOLMARKERS[caller],
                      edgecolor='black', linewidth=1, color='white', size=10, ax=axes[0])

    sns.stripplot(data=assm_ins_pcrt, x='dataset', y='pcrt', hue='aa',
                  hue_order=assemblers, palette=[ASSMBLERCOLOR[ele] for ele in assemblers], size=7, edgecolor='black', linewidth=1, ax=axes[1])

    for i, ax in enumerate(axes):
        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in DATASETS], fontsize=13, rotation=60)
        ax.set_xlabel('')

        ax.set_ylim(20, 80)
        ax.set_yticks(np.linspace(20, 80, 4))
        ax.set_yticklabels([int(val) for val in np.linspace(20, 80, 4)], fontsize=12)
        if i == 0:
            ax.set_ylabel('% of INS', fontsize=13)
            ax.legend('', frameon=False)
        else:
            ax.set_ylabel('')
            ax.legend(title='', loc='best')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs1a-2.pdf')
    # plt.show()


def plot_s1b(figdir):
    assm_svs_small = []
    assm_svs_large = []

    read_svs_small = []
    read_svs_large = []

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for aligner in READALIGNERS:
            for read_caller in READCALLERS:
                read_calls = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/HG002.{read_caller}.exbnd.bed'
                with open(read_calls, 'r') as f:
                    for line in f:
                        entries = line.strip().split('\t')
                        svsize = abs(int(entries[4]))
                        if svsize <= 1000:
                            read_svs_small.append((svsize, entries[3], read_caller, plat, 'Read-based'))
                        elif svsize <= 10000:
                            read_svs_large.append((svsize, entries[3], read_caller, plat, 'Read-based'))

            if aligner == 'minimap2':
                for asm_caller in ASMCALLERS:
                    asm_calls = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{asm_caller}.hifiasm.bed'
                    if 'ont' in dataset:
                        asm_calls = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{asm_caller}.flye.bed'

                    with open(asm_calls, 'r') as f:
                        for line in f:
                            entries = line.strip().split('\t')
                            svsize = abs(int(entries[4]))
                            if svsize <= 1000:
                                assm_svs_small.append((svsize, entries[3], asm_caller, plat, 'Assembly-based'))
                            elif svsize <= 10000:
                                assm_svs_large.append((svsize, entries[3], asm_caller, plat, 'Assembly-based'))

    df_readsvs_small = pd.DataFrame(read_svs_small, columns=['size', 'svtype', 'caller', 'plat', 'stra'])
    df_readsvs_large = pd.DataFrame(read_svs_large, columns=['size', 'svtype', 'caller', 'plat', 'stra'])

    df_assmsvs_small = pd.DataFrame(assm_svs_small, columns=['size', 'svtype', 'caller', 'plat', 'stra'])
    df_assmsvs_large = pd.DataFrame(assm_svs_large, columns=['size', 'svtype', 'caller', 'plat', 'stra'])

    svtype_order = ['DEL', 'INS']
    svtype_colors = [SVTYPECOLORS[ele] for ele in svtype_order]
    svtype_legend = [Patch(facecolor=SVTYPECOLORS[ele], label=ele) for ele in svtype_order]

    fig, axes = plt.subplots(2, 1, figsize=(7, 5))

    sns.histplot(data=df_assmsvs_small, x='size', hue='svtype', bins=70, hue_order=svtype_order, palette=svtype_colors,
                 alpha=1, ax=axes[0])

    axes[0].set_xlim(50, 1000)
    axes[0].set_xticks([50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axes[0].set_xticklabels(['50', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1,000'], fontsize=12,
                            rotation=60)

    axes[0].set_yticks([0, 10000, 20000])
    axes[0].set_yticklabels(['0', '10,000', '20,000'], fontsize=12)

    sns.histplot(data=df_assmsvs_large, x='size', hue='svtype', bins=70, hue_order=svtype_order, palette=svtype_colors,
                 alpha=1, ax=axes[1])

    axes[1].set_xlim(1000, 10000)
    axes[1].set_xticks([1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
    axes[1].set_xticklabels(['1,000', '2,000', '3,000', '4,000', '5,000', '6,000', '7,000', '8,000', '9,000', '10,000'],
                            fontsize=12, rotation=60)

    axes[1].set_yticks([0, 1000, 2000])
    axes[1].set_yticklabels(['0', '1,000', '2,000'])

    for ax in axes:
        ax.tick_params(axis='both', which='both', length=0)
        ax.set_xlabel('SV length(bp)', fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_ylabel('# of SVs', fontsize=13)
        ax.legend(handles=svtype_legend, title='Type', frameon=False)

    fig.tight_layout()

    fig.savefig(f'{figdir}/figs1b-1.pdf')

    fig, axes = plt.subplots(2, 1, figsize=(7, 5))
    sns.histplot(data=df_readsvs_small, x='size', hue='svtype', bins=70, hue_order=svtype_order, palette=svtype_colors,
                 alpha=1, ax=axes[0])

    axes[0].set_xlim(50, 1000)
    axes[0].set_xticks([50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axes[0].set_xticklabels(['50', '100', '200', '300', '400', '500', '600', '700', '800', '900', '1,000'], fontsize=12,
                            rotation=60)

    axes[0].set_yticks([0, 75000, 150000])
    axes[0].set_yticklabels(['0', '75,000', '150,000'], fontsize=12)

    sns.histplot(data=df_readsvs_large, x='size', hue='svtype', bins=70, hue_order=svtype_order, palette=svtype_colors,
                 alpha=1, ax=axes[1])

    axes[1].set_xlim(1000, 10000)
    axes[1].set_xticks([1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
    axes[1].set_xticklabels(['1,000', '2,000', '3,000', '4,000', '5,000', '6,000', '7,000', '8,000', '9,000', '10,000'],
                            fontsize=12, rotation=60)

    axes[1].set_yticks([0, 7500, 15000])
    axes[1].set_yticklabels(['0', '7,500', '15,000'], fontsize=12)

    for ax in axes:
        ax.tick_params(axis='both', which='both', length=0)
        ax.set_xlabel('SV length(bp)', fontsize=12)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        ax.set_ylabel('# of SVs', fontsize=13)
        ax.legend(handles=svtype_legend, title='Type', frameon=False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs1b-2.pdf')
    # plt.show()

def plot_s1c(figdir):
    assembly_counts = {region: [] for region in GENOMICREGIONS}
    read_counts = {region: [] for region in GENOMICREGIONS}

    read_insdel_counts = {'INS': {region: [] for region in GENOMICREGIONS},
                          'DEL': {region: [] for region in GENOMICREGIONS}}
    assm_insdel_counts = {'INS': {region: [] for region in GENOMICREGIONS},
                          'DEL': {region: [] for region in GENOMICREGIONS}}

    align_sv_count = pd.read_csv(f'{WORKDIR}/caller_sv_count.tsv', header=[0], sep='\t')

    caller_sv_totals = {}
    for idx, row in align_sv_count.iterrows():
        all_sv_num, ins_num, del_num, aligner, dataset, caller = int(row['all_num']), int(row['ins_num']), int(
            row['del_num']), row['aligner'], row['dataset'], row['caller']
        key = f'{caller}-{aligner}-{dataset}'
        caller_sv_totals[key] = all_sv_num

    pav_sv_count = pd.read_csv(f'{WORKDIR}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(
            row['total']), int(row['ins_num']), int(row['del_num'])
        key = f'{caller}-{assembler}-{dataset}'
        caller_sv_totals[key] = svcount

    svimasm_sv_count = pd.read_csv(f'{WORKDIR}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(
            row['total']), int(row['ins_num']), int(row['del_num'])
        key = f'{caller}-{assembler}-{dataset}'
        caller_sv_totals[key] = svcount

    sv_region_info = []
    df_sv_counts = pd.read_csv(f'{WORKDIR}/caller_sv_counts_region.tsv', sep='\t', header=0)
    read_caller_sv_region = {}
    for idx, row in df_sv_counts.iterrows():
        caller, dataset, aligner, region, count, svtype = row['caller'], row['dataset'], row['aligner'], row[
            'region'], int(row['count']), row['svtype']
        read_counts[region].append(count)

        if svtype in read_insdel_counts:
            read_insdel_counts[svtype][region].append(count)

        key = f'{caller}-{aligner}-{dataset}'
        if key in read_caller_sv_region:
            read_caller_sv_region[key][region] += count
        else:
            read_caller_sv_region[key] = {region: 0 for region in GENOMICREGIONS}
            read_caller_sv_region[key][region] += count

    for key, count_region in read_caller_sv_region.items():
        caller, aligner, dataset = key.split('-')
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for region, count in count_region.items():
            sv_region_info.append((caller, plat, dataset, aligner, count / caller_sv_totals[key] * 100, region, 'Read'))

    assm_caller_sv_region = {}
    pav_sv_count = pd.read_csv(f'{WORKDIR}/pav_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        dataset, region, svcount, assembler, svtype = row['dataset'], row['region'], int(row['count']), row[
            'assembler'], row['svtype']
        assembly_counts[region].append(svcount)

        if svtype in assm_insdel_counts:
            assm_insdel_counts[svtype][region].append(svcount)

        key = f'pav-{assembler}-{dataset}'
        if key in assm_caller_sv_region:
            assm_caller_sv_region[key][region] += svcount
        else:
            assm_caller_sv_region[key] = {region: 0 for region in GENOMICREGIONS}
            assm_caller_sv_region[key][region] += svcount

    svimasm_sv_count = pd.read_csv(f'{WORKDIR}/svimasm_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        dataset, region, svcount, assembler, svtype = row['dataset'], row['region'], int(row['count']), row[
            'assembler'], row['svtype']
        assembly_counts[region].append(svcount)
        key = f'svimasm-{assembler}-{dataset}'
        if key in assm_caller_sv_region:
            assm_caller_sv_region[key][region] += svcount
        else:
            assm_caller_sv_region[key] = {region: 0 for region in GENOMICREGIONS}
            assm_caller_sv_region[key][region] += svcount

    for key, count_region in assm_caller_sv_region.items():
        caller, assembler, dataset = key.split('-')
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for region, count in count_region.items():
            sv_region_info.append(
                (caller, plat, dataset, assembler, count / caller_sv_totals[key] * 100, region, 'Assembly'))

    region_avg = {}

    for region in GENOMICREGIONS:
        region_avg[region] = [np.mean(read_counts[region]), np.mean(assembly_counts[region])]

    size = 0.4

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    read_ax = axes[0]
    read_ax.set_title('Read-based', fontsize=13)
    read_insdel_avg = {}
    for region in GENOMICREGIONS:
        read_insdel_avg[region] = [np.mean(read_insdel_counts['INS'][region]),
                                   np.mean(read_insdel_counts['DEL'][region])]

    df_read_insdel_avg = pd.DataFrame(read_insdel_avg)

    totals = [i + j + k + m for i, j, k, m in
              zip(df_read_insdel_avg['Tandem Repeats'], df_read_insdel_avg['Repeat Masked'],
                  df_read_insdel_avg['Segment Dup'], df_read_insdel_avg['Simple Region'])]
    sim_pcrt = [i / j * 100 for i, j in zip(df_read_insdel_avg['Tandem Repeats'], totals)]
    rm_pcrt = [i / j * 100 for i, j in zip(df_read_insdel_avg['Repeat Masked'], totals)]
    sd_pcrt = [i / j * 100 for i, j in zip(df_read_insdel_avg['Segment Dup'], totals)]
    uni_pcrt = [i / j * 100 for i, j in zip(df_read_insdel_avg['Simple Region'], totals)]

    r = [0.2, 0.8]
    read_ax.bar(r, sim_pcrt, width=size, color=REGIONCOLORS['Tandem Repeats'])
    read_ax.bar(r, rm_pcrt, width=size, color=REGIONCOLORS['Repeat Masked'], bottom=sim_pcrt)
    read_ax.bar(r, sd_pcrt, width=size, color=REGIONCOLORS['Segment Dup'],
                bottom=[i + j for i, j in zip(sim_pcrt, rm_pcrt)])
    read_ax.bar(r, uni_pcrt, width=size, color=REGIONCOLORS['Simple Region'],
                bottom=[i + j + k for i, j, k in zip(sim_pcrt, rm_pcrt, sd_pcrt)])

    assm_ax = axes[1]
    assm_ax.set_title('Assembly-based', fontsize=13)
    assm_insdel_avg = {}
    for region in GENOMICREGIONS:
        assm_insdel_avg[region] = [np.mean(assm_insdel_counts['INS'][region]),
                                   np.mean(assm_insdel_counts['DEL'][region])]

    df_assm_insdel_avg = pd.DataFrame(assm_insdel_avg)

    totals = [i + j + k + m for i, j, k, m in
              zip(df_assm_insdel_avg['Tandem Repeats'], df_assm_insdel_avg['Repeat Masked'],
                  df_assm_insdel_avg['Segment Dup'], df_assm_insdel_avg['Simple Region'])]
    sim_pcrt = [i / j * 100 for i, j in zip(df_assm_insdel_avg['Tandem Repeats'], totals)]
    rm_pcrt = [i / j * 100 for i, j in zip(df_assm_insdel_avg['Repeat Masked'], totals)]
    sd_pcrt = [i / j * 100 for i, j in zip(df_assm_insdel_avg['Segment Dup'], totals)]
    uni_pcrt = [i / j * 100 for i, j in zip(df_assm_insdel_avg['Simple Region'], totals)]

    assm_ax.bar(r, sim_pcrt, width=size, color=REGIONCOLORS['Tandem Repeats'])
    assm_ax.bar(r, rm_pcrt, width=size, color=REGIONCOLORS['Repeat Masked'], bottom=sim_pcrt)
    assm_ax.bar(r, sd_pcrt, width=size, color=REGIONCOLORS['Segment Dup'],
                bottom=[i + j for i, j in zip(sim_pcrt, rm_pcrt)])
    assm_ax.bar(r, uni_pcrt, width=size, color=REGIONCOLORS['Simple Region'],
                bottom=[i + j + k for i, j, k in zip(sim_pcrt, rm_pcrt, sd_pcrt)])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.set_xticks(r)
        ax.set_xticklabels(['INS', 'DEL'], fontsize=13)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels(['0', '25', '50', '75', '100'], fontsize=12)

        ax.grid(axis='y', ls='--', color='grey', lw=1.5)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()

    fig.savefig(f'{figdir}/figs1c.pdf')
    # plt.show()

def plot_s1d(figdir):
    SVTYPECOLORS = {'DEL': '#dc9049', 'DUP': '#54b06d', 'INS': '#5090dd', 'INV': '#f2ed4f', 'BND': '#c265e7', 'Others': '#CDB699', 'DUP:TANDEM': '#CDB699'}

    df_sv_count = pd.read_csv(f'{WORKDIR}/caller_sv_count.tsv', header=[0], sep='\t')

    caller_sv_counts = {'INS': [0 for i in range(len(DATASETS))],
                        'DEL': [0 for i in range(len(DATASETS))],
                        'INV': [0 for i in range(len(DATASETS))],
                        'DUP': [0 for i in range(len(DATASETS))],
                        'Others': [0 for i in range(len(DATASETS))]}

    for idx, row in df_sv_count.iterrows():
        all_sv_num, ins_num, del_num, inv_num, dup_num, aligner, dataset, caller = int(row['all_num']), \
                                                                                   int(row['ins_num']), int(
            row['del_num']), int(row['inv_num']), int(row['dup_num']), row['aligner'], \
                                                                                   row['dataset'], row['caller']

        dataset_idx = DATASETS.index(dataset)

        caller_sv_counts['INS'][dataset_idx] += ins_num
        caller_sv_counts['DEL'][dataset_idx] += del_num
        caller_sv_counts['INV'][dataset_idx] += inv_num
        caller_sv_counts['DUP'][dataset_idx] += dup_num
        caller_sv_counts['Others'][dataset_idx] += all_sv_num - ins_num - del_num - inv_num - dup_num

    svtype_legends = [Patch(facecolor=SVTYPECOLORS[svtype], edgecolor='black', label=svtype) for svtype in
                      ['INS', 'DEL', 'INV', 'DUP', 'Others']]

    barwidth = 0.5
    xticks = np.arange(len(DATASETS))
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    ax.bar(xticks, caller_sv_counts['INS'], color=SVTYPECOLORS['INS'], edgecolor='black', width=barwidth)
    ax.bar(xticks, caller_sv_counts['DEL'], color=SVTYPECOLORS['DEL'], bottom=caller_sv_counts['INS'],
           edgecolor='black',
           width=barwidth)
    ax.bar(xticks, caller_sv_counts['DUP'], color=SVTYPECOLORS['DUP'],
           bottom=[i + j for i, j in zip(caller_sv_counts['INS'], caller_sv_counts['DEL'])], edgecolor='black',
           width=barwidth)
    ax.bar(xticks, caller_sv_counts['INV'], color=SVTYPECOLORS['INV'],
           bottom=[i + j + k for i, j, k in
                   zip(caller_sv_counts['INS'], caller_sv_counts['DEL'], caller_sv_counts['DUP'])],
           edgecolor='black', width=barwidth)
    ax.bar(xticks, caller_sv_counts['Others'], color=SVTYPECOLORS['Others'], bottom=[i + j + k + m for i, j, k, m in
                                                                                     zip(caller_sv_counts['INS'],
                                                                                         caller_sv_counts['DEL'],
                                                                                         caller_sv_counts['INV'],
                                                                                         caller_sv_counts['DUP'])],
           edgecolor='black', width=barwidth)

    ax.legend(handles=svtype_legends, ncol=2)
    ax.set_xticks(xticks)
    ax.set_xticklabels([PLATMAP[ele] for ele in DATASETS], rotation=60, fontsize=13)

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(axis='y', ls='--', color='grey', lw=2)

    ax.set_ylim(0, 500000)
    ax.set_yticks(np.linspace(0, 500000, 6))
    ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 50, 6)], fontsize=12)
    ax.set_ylabel('Number of SVs (x$10^4$)', fontsize=13)

    plt.tight_layout()
    fig.savefig(f'{figdir}/figs1d.pdf')
    # plt.show()

def main():
    print('\n==== Creating Extended Data Fig 1 =====')

    if not os.path.exists(f'{FIGDIR}/FigS1'):
        os.mkdir(f'{FIGDIR}/FigS1')

    plot_s1a(f'{FIGDIR}/FigS1')
    print(f'Figures saved to {FIGDIR}/FigS1/figs1a.pdf')

    plot_s1b(f'{FIGDIR}/FigS1')
    print(f'Figures saved to {FIGDIR}/FigS1/figs1b.pdf')

    plot_s1c(f'{FIGDIR}/FigS1')
    print(f'Figures saved to {FIGDIR}/FigS1/figs1c.pdf')

    plot_s1d(f'{FIGDIR}/FigS1')
    print(f'Figures saved to {FIGDIR}/FigS1/figs1d.pdf')


if __name__ == '__main__':
    main()
