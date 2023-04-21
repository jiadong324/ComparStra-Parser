#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/16

'''

import math
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pysam
from matplotlib.patches import Patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Constant import *
from Helpers.Functions import get_survivor_supp, get_overlaps


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


def prepare_data():
    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        outdir = f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for caller in ASMCALLERS:

            tmp_file = f'{outdir}/{caller}.{dataset}.txt'
            tmp_file_writer = open(tmp_file, 'w')
            for assembler in plat_assemblers[plat]:
                vcf_file = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                print(vcf_file, file=tmp_file_writer)
            tmp_file_writer.close()

            merged_out_vcf = f'{outdir}/{caller}.{dataset}.jasmine.merged.vcf'
            print(f'Obtain {caller} {dataset} assembler uniques and overlaps ...')
            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {outdir}/jasmine.assm.merge.log'
            os.system(cmd)
            get_assembler_compare_info(plat_assemblers[plat], merged_out_vcf, outdir, caller, simple_reps, rmsk, sds)
            os.remove(tmp_file)

        for caller in READCALLERS:

            tmp_file = f'{outdir}/{caller}.{dataset}.txt'
            tmp_file_writer = open(tmp_file, 'w')

            for aligner in READALIGNERS:
                vcf_file = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/HG002.{caller}.vcf'
                print(f'{vcf_file}', file=tmp_file_writer)
            tmp_file_writer.close()

            merged_out_vcf = f'{outdir}/{caller}.{dataset}.jasmine.merged.vcf'
            print(f'Obtain {caller} {dataset} aligner uniques and overlaps ...')

            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1 > {outdir}/jasmine.read.merge.log'
            os.system(cmd)

            get_aligner_compare_info(READALIGNERS, merged_out_vcf, outdir, caller, simple_reps, rmsk, sds)
            os.remove(tmp_file)

def plot_2d(figdir):

    concordants_pcrt = []

    for caller in READCALLERS:
        wgs_supp4_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}
        extd_supp4_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}

        for dataset_idx, dataset in enumerate(DATASETS):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            this_dataset_index = DATASET_DICT[plat].index(dataset)
            unique_info = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.aligner-unique.exTD.tsv', header=[0], sep='\t')
            matched_info = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.aligner-concordant.exTD.tsv', header=[0], sep='\t')

            merged_total = len(unique_info) + len(matched_info)

            supp_dict = {}
            for idx, row in matched_info.iterrows():
                supp = int(row['SUPP'])
                if supp in supp_dict:
                    supp_dict[supp] += 1
                else:
                    supp_dict[supp] = 1

            concordants_pcrt.append((supp_dict[4] / merged_total * 100, caller, plat, 'ExTD'))
            extd_supp4_pcrt[plat][this_dataset_index] += supp_dict[4] / merged_total * 100

            wgs_merged_vcf = f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(wgs_merged_vcf)

            concordants_pcrt.append((wgs_supp_dict[4] / wgs_merged_total * 100, caller, plat, 'WGS'))
            wgs_supp4_pcrt[plat][this_dataset_index] += wgs_supp_dict[4] / wgs_merged_total * 100


    df_concordants = pd.DataFrame(concordants_pcrt, columns=['pcrt', 'caller', 'plat', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    fig, ax = plt.subplots(1, 1, figsize=(3, 4))

    sns.barplot(data=df_concordants, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                 palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=ax)


    ax.set_ylabel('% of aligner concordant SVs', fontsize=13)

    ax.set_xlabel('')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['ExTD', 'WGS'], fontsize=13)

    ax.legend(handles=plat_legends)
    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

    ax.grid(axis='y', ls='--', color='grey', lw=1.5)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()

    fig.savefig(f'{figdir}/fig2d-1.pdf')

    unique_pcrt = []
    concordants_pcrt = []

    for caller in ASMCALLERS:
        wgs_supp2_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}
        extd_supp2_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}

        for dataset_idx, dataset in enumerate(DATASETS):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            unique_info = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.assembler-unique.exTD.tsv', header=[0], sep='\t')
            matched_info = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.assembler-concordant.exTD.tsv', header=[0], sep='\t')

            merged_total = len(unique_info) + len(matched_info)
            supp_dict = {}

            for idx, row in matched_info.iterrows():
                supp = int(row['SUPP'])
                if supp in supp_dict:
                    supp_dict[supp] += 1
                else:
                    supp_dict[supp] = 1

            this_dataset_index = DATASET_DICT[plat].index(dataset)
            unique_pcrt.append((len(unique_info) / merged_total * 100, plat, caller, 'ExTD'))
            concordants_pcrt.append((supp_dict[2] / merged_total * 100, plat, caller, 'ExTD'))

            extd_supp2_pcrt[plat][this_dataset_index] += supp_dict[2] / merged_total * 100

            merged_vcf = f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(merged_vcf)

            unique_pcrt.append((wgs_supp_dict[1] / wgs_merged_total * 100, plat, caller, 'WGS'))
            concordants_pcrt.append((wgs_supp_dict[2] / wgs_merged_total * 100, plat, caller, 'WGS'))
            wgs_supp2_pcrt[plat][this_dataset_index] += wgs_supp_dict[2] / wgs_merged_total * 100

    df_concordants = pd.DataFrame(concordants_pcrt, columns=['pcrt', 'plat', 'caller', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    fig, ax = plt.subplots(1, 1, figsize=(3, 4))
    sns.barplot(data=df_concordants, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=ax)

    ax.set_ylabel('% of assembler concordant SVs', fontsize=13)
    ax.set_xlabel('')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['ExTD', 'WGS'], fontsize=13)

    ax.legend(handles=plat_legends)
    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

    ax.grid(axis='y', ls='--', color='grey', lw=1.5)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig2d-2.pdf')
    # plt.show()

def plot_2e(figdir):
    leftbp_pcrt_list = []
    shift_labels = ['0', '0,10']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']

    for platform_idx, dataset in enumerate(DATASETS):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for caller_idx, caller in enumerate(READCALLERS):
            wgs_shift_label = {shift: 0 for shift in shift_labels}
            df_matched = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.aligner-concordant.tsv', sep='\t', header=[0])
            this_total = 0

            for idx, row in df_matched.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 4:
                    this_total += 1
                    if start_std == 0:
                        wgs_shift_label['0'] += 1
                    elif start_std <= 10:
                        wgs_shift_label['0,10'] += 1

            df_extd_match_info = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.aligner-concordant.exTD.tsv', sep='\t', header=[0])
            extd_total = 0
            extd_shift_label = {shift: 0 for shift in shift_labels}
            for idx, row in df_extd_match_info.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 4:
                    extd_total += 1
                    if start_std == 0:
                        extd_shift_label['0'] += 1
                    elif start_std <= 10:
                        extd_shift_label['0,10'] += 1

            for shift in shift_labels:
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller], wgs_shift_label[shift] / this_total * 100, shift, 'WGS'))
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller], extd_shift_label[shift] / extd_total * 100, shift, 'ExTD'))

        for caller in ASMCALLERS:
            wgs_shift_label = {shift: 0 for shift in shift_labels}
            df_matched = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.assembler-concordant.tsv', sep='\t', header=[0])
            this_total = 0

            for idx, row in df_matched.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 2:
                    this_total += 1
                    if start_std == 0:
                        wgs_shift_label['0'] += 1
                    elif start_std <= 10:
                        wgs_shift_label['0,10'] += 1

            df_extd_matched = pd.read_csv( f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.assembler-concordant.exTD.tsv', sep='\t', header=[0])
            extd_total = 0
            extd_shift_label = {shift: 0 for shift in shift_labels}
            for idx, row in df_extd_matched.iterrows():
                start_std = abs(float(row['START_STD']))
                if int(row['SUPP']) == 2:
                    extd_total += 1
                    if start_std == 0:
                        extd_shift_label['0'] += 1
                    elif start_std <= 10:
                        extd_shift_label['0,10'] += 1

            for shift in shift_labels:
                leftbp_pcrt_list.append(
                    (plat, PLATMAP[dataset], TOOLMAP[caller], wgs_shift_label[shift] / this_total * 100, shift, 'WGS'))
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller],
                                         extd_shift_label[shift] / extd_total * 100, shift, 'ExTD'))

    df_bpstd = pd.DataFrame(leftbp_pcrt_list, columns=['plat', 'dataset', 'caller', 'pcrt', 'shift', 'region'])
    df_bpstd_subset1 = df_bpstd[df_bpstd['shift'] == '0']

    # caller_legends = [Line2D([0], [0], color='white', mfc=TOOLCOLORS[ele], marker='o', markersize=8, label=TOOLMAP[ele])
    #                   for ele in callers]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))
    sns.stripplot(data=df_bpstd_subset1[df_bpstd_subset1['region'] == 'WGS'], x='dataset', y='pcrt', hue='caller', marker='o', lw=2, size=8,
                  hue_order=[TOOLMAP[ele] for ele in callers], palette=[TOOLCOLORS[ele] for ele in callers], ax=axes[0])
    sns.stripplot(data=df_bpstd_subset1[df_bpstd_subset1['region'] == 'ExTD'], x='dataset', y='pcrt', hue='caller', marker='o', lw=2, size=8,
                  hue_order=[TOOLMAP[ele] for ele in callers], palette=[TOOLCOLORS[ele] for ele in callers], ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of BSD-0 SVs', fontsize=13)
            # ax.legend(handles=caller_legends)
            ax.legend('', frameon=False)
        else:
            ax.set_ylabel('')
            ax.legend('', frameon=False)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

        ax.set_xlabel('')
        ax.set_xticks(np.arange(len(DATASETS)))
        ax.set_xticklabels([PLATMAP[ele] for ele in DATASETS], fontsize=13, rotation=60)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig2e.pdf')
    # plt.show()

def plot_2f2g(figdir):

    aligner_unique_counts = []
    aligner_unique_info = []
    aligner_uniques = {aligner: [] for aligner in READALIGNERS}

    unique_size_in_range = {aligner: 0 for aligner in READALIGNERS}

    for dataset_idx, dataset in enumerate(DATASETS):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for j, caller in enumerate(READCALLERS):

            unique_info_out = f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.aligner-unique.exTD.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])
            aligner_unique_dict = {aligner: 0 for aligner in READALIGNERS}
            aligner_unique_region_dict = {aligner: [0 for i in range(len(GENOMICREGIONS))] for aligner in READALIGNERS}

            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row['REGION_TYPE']
                region_idx = GENOMICREGIONS.index(sv_region)
                aligner_unique_dict[aligner] += 1
                aligner_unique_region_dict[aligner][region_idx] += 1
                if svlen >= 50:
                    aligner_unique_info.append((svlen, aligner, sv_region, svtype, plat))

                if svlen >= 100 and svlen <= 1000:
                    unique_size_in_range[aligner] += 1

            for aligner, count in aligner_unique_dict.items():
                aligner_unique_counts.append((count, plat, PLATMAP[dataset], caller, aligner))
                aligner_uniques[aligner].append(count)

    df_aligner_unique = pd.DataFrame(aligner_unique_info, columns=['svlen', 'aligner', 'region', 'svtype', 'plat'])
    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(6, 4))

    aligner_order = ['winnowmap', 'lra', 'ngmlr', 'minimap2']
    aligner_legends = [Patch(facecolor=ALIGNERCOLOR[aligner], label=aligner) for aligner in aligner_order]

    sns.histplot(data=df_aligner_unique[df_aligner_unique['plat'] == 'HiFi'], x='svlen', hue='aligner', bins=50,
                 alpha=1,hue_order=['minimap2', 'winnowmap', 'lra', 'ngmlr'], palette=[ALIGNERCOLOR[ele] for ele in ['minimap2', 'winnowmap', 'lra', 'ngmlr']], log_scale=True,
                 ax=axes[0])

    axes[0].set_ylabel('HiFi (x$10^3$)', fontsize=13)

    sns.histplot(data=df_aligner_unique[df_aligner_unique['plat'] == 'ONT'], x='svlen', hue='aligner', bins=50, alpha=1,
                 hue_order=['winnowmap', 'minimap2', 'ngmlr', 'lra'], palette=[ALIGNERCOLOR[ele] for ele in ['winnowmap', 'minimap2', 'ngmlr', 'lra']], log_scale=True, ax=axes[1])

    axes[1].set_ylabel('ONT (x$10^3$)', fontsize=13)

    for ax in axes:
        ax.legend(handles=aligner_legends)
        ax.tick_params(axis='both', which='both', length=0)

        ax.set_xlabel('SV length(bp)', fontsize=12)
        ax.set_xlim(50, 100000)
        ax.set_xticks([50, 100, 300, 1000, 6000, 100000])
        ax.set_xticklabels(['50', '100', '300', '1,000', '6,000', '100,000'], fontsize=12, rotation=60)

        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_ylabel('Number of SVs', fontsize=13)

        ax.set_xlabel('')
        # ax.set_ylim(0, 1000)
        ax.set_yticks(np.linspace(0, 1000, 3))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 1000, 3)], fontsize=12)

    fig.tight_layout()
    # plt.show()
    fig.savefig(f'{figdir}/fig2f.pdf')

    svtypes = ['INS', 'DUP', 'DEL']

    svtypes_pcrt = {aligner: {region: {svtype: [] for svtype in svtypes} for region in GENOMICREGIONS} for aligner in
                    READALIGNERS}

    for dataset_idx, dataset in enumerate(DATASETS):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(READCALLERS):
            unique_info_out = f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.aligner-unique.exTD.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])

            this_total = {aligner: 0 for aligner in READALIGNERS}
            svtypes_count = {aligner: {region: [0, 0, 0] for region in GENOMICREGIONS} for aligner in READALIGNERS}
            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row[
                    'REGION_TYPE']

                if svlen >= 100 and svlen <= 1000:
                    this_total[aligner] += 1
                    if svtype in svtypes:
                        svtype_idx = svtypes.index(svtype)
                        svtypes_count[aligner][sv_region][svtype_idx] += 1

            for aligner, region_counts in svtypes_count.items():
                for region, counts in region_counts.items():
                    for i, val in enumerate(counts):
                        svtype = svtypes[i]
                        svtypes_pcrt[aligner][region][svtype].append(val / this_total[aligner] * 100)

    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(7, 4))
    xticks = np.arange(len(svtypes))
    bar_width = 0.6

    for col_idx, aligner in enumerate(READALIGNERS):
        ax = axes[col_idx]
        for j, region in enumerate(GENOMICREGIONS):
            this_pcrt_avg = []
            for svtype in svtypes:
                this_pcrt_avg.append(np.mean(svtypes_pcrt[aligner][region][svtype]))

            if j == 0:
                ax.bar(xticks, this_pcrt_avg, color=REGIONCOLORS[region], width=bar_width, label=region)

            else:
                bottoms = []
                for k in range(0, j):
                    this_avg = []
                    for svtype in svtypes:
                        this_avg.append(np.mean(svtypes_pcrt[aligner][GENOMICREGIONS[k]][svtype]))
                    bottoms.append(this_avg)
                bottom_sum = [sum(x) for x in zip(*bottoms)]
                ax.bar(xticks, this_pcrt_avg, bottom=bottom_sum, width=bar_width, color=REGIONCOLORS[region],
                       label=region)

        ax.set_xticks(xticks)
        ax.set_xticklabels(svtypes, fontsize=13, rotation=60)

        ax.set_title(aligner, fontsize=13)

        if col_idx == 0:
            ax.set_ylabel('% of aligner unique SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(0, 80)
        ax.set_yticks(np.linspace(0, 80, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 80, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig2g.pdf')

def get_assembler_compare_info(assemblers, merged_vcf, compare_outdir, caller, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    extd_matched_list = []
    extd_unique_list = []
    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')

            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            merged_type = info_dict['SVTYPE']
            merged_id = entries[2]
            supp_vec = info_dict['SUPP_VEC']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])


            # region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 50, simple_reps, rmsk, sds)
            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 50, simple_reps, rmsk, sds)

            if supp > 1:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_assembler = assemblers[supp_vec.index('1')]
                unique_list.append((merged_id, merged_type, unique_assembler, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_unique_list.append((merged_id, merged_type, unique_assembler, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

    ## At whole genome scale, comparing SV detect by a caller on different assembler
    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{caller}.assembler-concordant.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'ASSEMBLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{caller}.assembler-unique.tsv', sep='\t', header=True, index=False)

    ## Outside of tandem repeat regions, comparing SV detect by a caller on different assembler
    df_extd_matched = pd.DataFrame(extd_matched_list,columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{caller}.assembler-concordant.exTD.tsv', header=True, sep='\t', index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list,columns=['ID_MATCH', 'TYPE_MATCH', 'ASSEMBLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{caller}.assembler-unique.exTD.tsv', sep='\t', header=True, index=False)


def get_aligner_compare_info(aligners, merged_vcf, compare_outdir, caller, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    extd_matched_list = []
    extd_unique_list = []

    svs_by_regions = {rep: [0 for i in range(len(aligners))] for rep in GENOMICREGIONS}

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')

            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            merged_type = info_dict['SVTYPE']
            merged_id = entries[2]
            supp_vec = info_dict['SUPP_VEC']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region1(entries[0], int(entries[1]), end, 50, simple_reps, rmsk, sds)

            if supp > 1:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_aligner = aligners[supp_vec.index('1')]
                unique_list.append((merged_id, merged_type, unique_aligner, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))
                if region_label != 'Tandem Repeats':
                    extd_unique_list.append((merged_id, merged_type, unique_aligner, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

                if supp_vec == '01':
                    svs_by_regions[region_label][2] += 1
                if supp_vec == '10':
                    svs_by_regions[region_label][0] += 1

    ## At whole genome scale, comparing SV detect by a caller on different alignments
    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{caller}.aligner-concordant.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{caller}.aligner-unique.tsv', sep='\t', header=True, index=False)

    ## Outside of the tandem repeat regions, comparing SV detect by a caller on different alignments
    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{caller}.aligner-concordant.exTD.tsv', header=True, sep='\t', index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{caller}.aligner-unique.exTD.tsv', sep='\t', header=True, index=False)


def annotate_sv_region(chrom, start, end, pcrt_thresh, simreps_tabix, rmsk_tabix, sd_tabix):
    if start > end:
        start, end = end, start
    size = end - start + 1
    annotations = []

    if 'chr' not in chrom:
        chrom = f'chr{chrom}'

    for simrep in simreps_tabix.fetch(chrom, start, end):
        entries = simrep.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)

        if overlap_size > 0:
            motif = rp_info.split(',')[-1]
            overlap_pcrt = min(overlap_size / size * 100, 100)
            subtype = 'VNTR' if len(motif) >= 7 else 'STR'
            if overlap_pcrt >= pcrt_thresh:
                annotations.append(('Tandem Repeats', subtype, overlap_pcrt))


    for rmsk in rmsk_tabix.fetch(chrom, start, end):
        entries = rmsk.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[4]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            rptype = rp_info.split(',')[11]
            if overlap_pcrt >= pcrt_thresh:
                if rptype == 'Simple_repeat':
                    motif = rptype[1: -2]
                    subtype = 'VNTR' if len(motif) >= 7 else 'STR'
                    annotations.append(('Tandem Repeats', subtype, overlap_pcrt))
                    continue
                annotations.append(('Repeat Masked', rptype, overlap_pcrt))

    for sd in sd_tabix.fetch(chrom, start, end):
        entries = sd.strip().split('\t')
        sd_start, sd_end, sd_mate_coord = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, sd_start, sd_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            annotations.append(('Segment Dup', 'SegDup', overlap_pcrt))

    if len(annotations) == 0:
        return ('Simple Region', 'None', 0)

    sorted_annotations = sorted(annotations, key=lambda x:x[1], reverse=True)

    return sorted_annotations[0]

def annotate_sv_region1(chrom, start, end, pcrt_thresh, simreps_tabix, rmsk_tabix, sd_tabix):
    if start > end:
        start, end = end, start
    size = end - start + 1
    annotations = []

    if 'chr' not in chrom:
        chrom = f'chr{chrom}'

    for simrep in simreps_tabix.fetch(chrom, start, end):
        entries = simrep.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)

        if overlap_size > 0:
            motif = rp_info.split(',')[-1]
            overlap_pcrt = min(overlap_size / size * 100, 100)
            subtype = 'VNTR' if len(motif) >= 7 else 'STR'
            if overlap_pcrt >= pcrt_thresh:
                annotations.append(('Tandem Repeats', subtype, overlap_pcrt))
                return ('Tandem Repeats', subtype, overlap_pcrt)

    for rmsk in rmsk_tabix.fetch(chrom, start, end):
        entries = rmsk.strip().split('\t')
        rp_start, rp_end, rp_info = int(entries[1]), int(entries[2]), entries[4]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            rptype = rp_info.split(',')[11]
            if overlap_pcrt >= pcrt_thresh:
                if rptype == 'Simple_repeat':
                    motif = rptype[1: -2]
                    subtype = 'VNTR' if len(motif) >= 7 else 'STR'
                    annotations.append(('Tandem Repeats', subtype, overlap_pcrt))
                    return ('Tandem Repeats', subtype, overlap_pcrt)
                annotations.append(('Repeat Masked', rptype, overlap_pcrt))

    for sd in sd_tabix.fetch(chrom, start, end):
        entries = sd.strip().split('\t')
        sd_start, sd_end, sd_mate_coord = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, sd_start, sd_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            annotations.append(('Segment Dup', 'SegDup', overlap_pcrt))

    if len(annotations) == 0:
        return ('Simple Region', 'None', 0)

    sorted_annotations = sorted(annotations, key=lambda x:x[1], reverse=True)

    return sorted_annotations[0]


def main():
    print('\n===== Prepare data for Figure 2d and 2e =====')
    print(f'Intermediate file directory: {WORKDIR}/dataset_fig2d2e_tmpfile')
    prepare_data()


    if not os.path.exists(f'{FIGDIR}/Fig2'):
        os.mkdir(f'{FIGDIR}/Fig2')

    print('\n===== Creating Figure2d and Figure2e =====')

    plot_2d(f'{FIGDIR}/Fig2')
    print(f'Figures saved to {FIGDIR}/Fig2/fig2d-1.pdf; {FIGDIR}/Fig2/fig2d-2.pdf')

    plot_2e(f'{FIGDIR}/Fig2')
    print(f'Figures saved to {FIGDIR}/Fig2/fig2e.pdf')

    print('\n===== Creating Figure2f and Figure2g =====')

    plot_2f2g(f'{FIGDIR}/Fig2')
    print(f'Figures saved to {FIGDIR}/Fig2/fig2f.pdf; {FIGDIR}/Fig2/fig2g.pdf')

if __name__ == '__main__':

    main()