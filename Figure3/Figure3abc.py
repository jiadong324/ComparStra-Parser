#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/13

'''

import math
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pysam
from matplotlib.lines import Line2D

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from Helpers.Constant import *
from Helpers.Functions import get_overlaps, read_suppvec_info, write_suppvec_info

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


def plot_3a3b(figdir):

    wgs_stra_pcrts = []
    extd_stra_pcrts = []

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for col_idx, aligner in enumerate(READALIGNERS):
        stra_matched = [[0 for i in range(len(DATASETS))] for j in range(3)]
        for dataset_idx, dataset in enumerate(DATASETS):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for caller in READCALLERS:
                for asm_caller in ASMCALLERS:
                    for assembler in plat_assemblers[plat]:
                        path = 'fig3a3b_tmpfile'
                        suppvec_info = f'{WORKDIR}/{plat}/{aligner}_{dataset}/{path}/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.jasmine.suppinfo.tsv'
                        suppvec_dict, merged_total = read_suppvec_info(suppvec_info)

                        for suppvec, count in suppvec_dict.items():
                            if suppvec == '10':
                                stra_matched[0][dataset_idx] += count
                            elif suppvec == '01':
                                stra_matched[2][dataset_idx] += count
                            else:
                                stra_matched[1][dataset_idx] += count

                        wgs_stra_pcrts.append((plat, suppvec_dict['11'] / (suppvec_dict['11'] + suppvec_dict['10']) * 100, 'Assembly', PLATMAP[dataset], aligner, assembler, TOOLMAP[caller], TOOLMAP[asm_caller]))
                        wgs_stra_pcrts.append((plat, suppvec_dict['11'] / (suppvec_dict['11'] + suppvec_dict['01']) * 100, 'Read', PLATMAP[dataset], aligner, assembler, TOOLMAP[caller], TOOLMAP[asm_caller]))

                        df_extd_matched = pd.read_csv(f'{WORKDIR}/{plat}/{aligner}_{dataset}/{path}/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.exTD.tsv', header=[0], sep='\t')
                        df_extd_uniques = pd.read_csv(f'{WORKDIR}/{plat}/{aligner}_{dataset}/{path}/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.uniques.info.exTD.tsv', header=[0], sep='\t')

                        suppvec_dict = {'10': 0, '01': 0}
                        for idx, row in df_extd_uniques.iterrows():
                            caller = row['CALLER']

                            if caller in ASMCALLERS:
                                suppvec_dict['10'] += 1
                            else:
                                suppvec_dict['01'] += 1

                        extd_stra_pcrts.append((plat, len(df_extd_matched) / (len(df_extd_matched) + suppvec_dict['10']) * 100, 'Assembly', PLATMAP[dataset], aligner, assembler, TOOLMAP[caller], TOOLMAP[asm_caller]))
                        extd_stra_pcrts.append((plat, len(df_extd_matched) / (len(df_extd_matched) + suppvec_dict['01']) * 100, 'Read', PLATMAP[dataset], aligner, assembler, TOOLMAP[caller], TOOLMAP[asm_caller]))

    df_wgs_stra_pcrt = pd.DataFrame(wgs_stra_pcrts, columns=['plat', 'pcrt', 'stra', 'dataset', 'aligner', 'assembler', 'caller','asmcaller'])
    df_extd_stra_pcrt = pd.DataFrame(extd_stra_pcrts, columns=['plat', 'pcrt', 'stra', 'dataset', 'aligner', 'assembler', 'caller', 'asmcaller'])

    ## Figure 3a

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    stra_legends = [Line2D([0], [0], label='ExTD', color='black', ls='--', lw=2),
                    Line2D([0], [0], label='WGS', color='black', lw=2),
                    Line2D([0], [0], label='Assembly', lw=2, color=STRACOLORS['Assembly']),
                    Line2D([0], [0], label='Read', lw=2, color=STRACOLORS['Read'])]

    sns.lineplot(data=df_wgs_stra_pcrt, x='dataset', y='pcrt', hue='stra', hue_order=['Assembly', 'Read'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']],
                 ci=None, lw=2, marker='o', ax=ax)
    sns.lineplot(data=df_extd_stra_pcrt, x='dataset', y='pcrt', hue='stra', hue_order=['Assembly', 'Read'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']],
                 ci=None, lw=2, ls='--', marker='o', ax=ax)

    ax.set_ylabel('% of overlapped SVs', fontsize=13)
    ax.set_ylim(50, 90)
    ax.set_yticks(np.linspace(50, 90, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(50, 90, 5)], fontsize=12)

    ax.set_xlabel('')
    ax.set_xticks(np.arange(len(DATASETS)))
    ax.set_xticklabels([PLATMAP[ele] for ele in DATASETS], fontsize=13, rotation=90)

    ax.legend(handles=stra_legends)

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig3a.pdf')

    ## Figure 3b
    fig, axes = plt.subplots(2, 4, sharex='col', sharey='row', figsize=(7, 4))

    fig_idx = 0
    for assembler in plat_assemblers['HiFi']:
        axes[0][fig_idx].set_title(assembler, fontsize=13)
        wgs_hifi_pcrt = df_wgs_stra_pcrt[df_wgs_stra_pcrt['plat'] == 'HiFi']
        sns.lineplot(data=wgs_hifi_pcrt[wgs_hifi_pcrt['assembler'] == assembler], y='pcrt', x='dataset', hue='aligner',
                     hue_order=READALIGNERS, palette=[ALIGNERCOLOR[ele] for ele in READALIGNERS], lw=2, marker='o', ci=None,
                     ax=axes[0][fig_idx])

        extd_hifi_pcrt = df_extd_stra_pcrt[df_extd_stra_pcrt['plat'] == 'HiFi']
        sns.lineplot(data=extd_hifi_pcrt[extd_hifi_pcrt['assembler'] == assembler], y='pcrt', x='dataset',
                     hue='aligner',
                     hue_order=READALIGNERS, palette=[ALIGNERCOLOR[ele] for ele in READALIGNERS], lw=2, ls='--', marker='o',
                     ci=None, ax=axes[1][fig_idx])

        fig_idx += 1

    for assembler in plat_assemblers['ONT']:
        axes[0][fig_idx].set_title(assembler, fontsize=13)
        wgs_ont_pcrt = df_wgs_stra_pcrt[df_wgs_stra_pcrt['plat'] == 'ONT']
        sns.lineplot(data=wgs_ont_pcrt[wgs_ont_pcrt['assembler'] == assembler], y='pcrt', x='dataset', hue='aligner',
                     hue_order=READALIGNERS, palette=[ALIGNERCOLOR[ele] for ele in READALIGNERS], lw=2, marker='o', ci=None,
                     ax=axes[0][fig_idx])

        extd_ont_pcrt = df_extd_stra_pcrt[df_extd_stra_pcrt['plat'] == 'ONT']
        sns.lineplot(data=extd_ont_pcrt[extd_ont_pcrt['assembler'] == assembler], y='pcrt', x='dataset', hue='aligner',
                     hue_order=READALIGNERS, palette=[ALIGNERCOLOR[ele] for ele in READALIGNERS], lw=2, ls='--', marker='o',
                     ci=None, ax=axes[1][fig_idx])

        fig_idx += 1

    for row_idx in range(axes.shape[0]):
        for col_idx in range(axes.shape[1]):
            ax = axes[row_idx][col_idx]

            ax.set_xlabel('')
            ax.set_xticks([0, 1, 2])

            if col_idx == 0:
                if row_idx == 0:
                    ax.set_ylabel('WGS', fontsize=13)
                else:
                    ax.set_ylabel('ExTD', fontsize=13)
            else:
                ax.set_ylabel('')

            if col_idx < 2:
                ax.set_xticklabels(['HiFi-10kb', 'HiFi-15kb', 'HiFi-18kb'], rotation=60, fontsize=13)
            else:
                ax.set_xticklabels(['ONT-9kb', 'ONT-19kb', 'ONT-30kb'], rotation=60, fontsize=13)

            ax.legend('', frameon=False)
            ax.set_ylim(50, 90)
            ax.set_yticks(np.linspace(50, 90, 3))
            ax.set_yticklabels([int(val) for val in np.linspace(50, 90, 3)], fontsize=12)

            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig3b.pdf')

    # plt.show()

def plot_3c(figdir):

    dataset = 'hifi_18kb'
    aligners = ['minimap2', 'winnowmap', 'ngmlr', 'lra']

    bpstd_pcrt_info = []
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    plat = 'HiFi'
    if 'ont' in DATASETS:
        plat = 'ONT'
    for col_idx, aligner in enumerate(aligners):
        for assembler in plat_assemblers[plat]:
            for caller in READCALLERS:
                for asm_caller in ASMCALLERS:
                    ins_bpstd_dict = {'Acc': 0, 'Inacc': 0}
                    del_bpstd_dict = {'Acc': 0, 'Inacc': 0}

                    ins_total = 0
                    del_total = 0

                    matched_file = f'{WORKDIR}/{plat}/{aligner}_{dataset}/fig3a3b_tmpfile/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.exTD.tsv'
                    df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                    for idx, row in df_matched.iterrows():
                        start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                        if start_std <= 10:
                            if svtype == 'INS':
                                ins_total += 1
                                ins_bpstd_dict['Acc'] += 1
                            if svtype == 'DEL':
                                del_total += 1
                                del_bpstd_dict['Acc'] += 1
                        else:
                            if svtype == 'INS':
                                ins_total += 1
                            if svtype == 'DEL':
                                del_total += 1

                    bpstd_pcrt_info.append(('INS', ins_bpstd_dict['Acc'] / ins_total * 100, assembler, aligner, caller, asm_caller))
                    bpstd_pcrt_info.append(('DEL', del_bpstd_dict['Acc'] / del_total * 100, assembler, aligner, caller, asm_caller))


    df_bpstd = pd.DataFrame(bpstd_pcrt_info, columns=['svtype', 'pcrt', 'assembler', 'aligner', 'caller', 'asmcaller'])
    df_ins_bpstd = df_bpstd[df_bpstd['svtype'] == 'INS']
    df_del_bpstd = df_bpstd[df_bpstd['svtype'] == 'DEL']


    fig, axes = plt.subplots(2, 4, sharey='row', sharex='col', figsize=(7, 4))


    for i, aligner in enumerate(aligners):

        sns.lineplot(data=df_ins_bpstd[df_ins_bpstd['aligner'] == aligner], x='caller', y='pcrt', marker='o',
                     hue='asmcaller', hue_order=['pav', 'svimasm'], palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']],
                     style='assembler', ax=axes[0][i])

        axes[0][i].set_title(aligner, fontsize=13)
        axes[0][i].spines['left'].set_visible(False)
        axes[0][i].spines['right'].set_visible(False)
        axes[0][i].legend('', frameon=False)

        axes[0][i].set_ylabel('% of BSD-10 INS', fontsize=13)
        axes[0][i].set_ylim(80, 100)
        axes[0][i].set_yticks(np.linspace(80, 100, 3))
        axes[0][i].set_yticklabels([int(ele) for ele in np.linspace(80, 100, 3)], fontsize=12)

        sns.lineplot(data=df_del_bpstd[df_del_bpstd['aligner'] == aligner], x='caller', y='pcrt', marker='o',
                     hue='asmcaller', hue_order=['pav', 'svimasm'], palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']],
                     style='assembler', ax=axes[1][i])

        axes[1][i].set_xlabel('')
        axes[1][i].set_xticks(np.arange(len(READCALLERS)))
        axes[1][i].set_xticklabels([TOOLMAP[ele] for ele in READCALLERS], rotation=90, fontsize=13)

        axes[1][i].spines['left'].set_visible(False)
        axes[1][i].spines['right'].set_visible(False)
        axes[1][i].legend('', frameon=False)

        axes[1][i].set_ylim(80, 100)
        axes[1][i].set_ylabel('% of BSD-10 DEL', fontsize=13)
        axes[1][i].set_yticks(np.linspace(80, 100, 3))
        axes[1][i].set_yticklabels([int(ele) for ele in np.linspace(80, 100, 3)], fontsize=12)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig3c.pdf')

    # plt.show()


def prepare_data():
    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')
    regioned_svs_counts = []

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for assembler in plat_assemblers[plat]:
            for asm_method in ASMCALLERS:
                asm_calls = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{asm_method}.{assembler}.vcf'
                for caller in READCALLERS:
                    for aligner in READALIGNERS:
                        print(f'Comparing {asm_method}-{assembler} to {caller}-{aligner} on {dataset} ...')
                        read_calls = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/HG002.{caller}.vcf'
                        compare_outdir = f'{WORKDIR}/{plat}/{aligner}_{dataset}/fig3a3b_tmpfile'

                        if not os.path.exists(compare_outdir):
                            os.mkdir(compare_outdir)
                        tmp_file = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.txt'
                        tmp_file_writer = open(tmp_file, 'w')
                        print(asm_calls, file=tmp_file_writer)
                        print(read_calls, file=tmp_file_writer)

                        tmp_file_writer.close()

                        merged_vcf = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{dataset}.jasmine.merged.vcf'
                        cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {compare_outdir}/jasmine.caller.merge.log'
                        os.system(cmd)

                        os.remove(tmp_file)
                        supp_info = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{dataset}.jasmine.suppinfo.tsv'
                        write_suppvec_info(merged_vcf, supp_info)

                        svs_by_regions = get_caller_compare_info([asm_method, caller], merged_vcf, compare_outdir, caller, asm_method, aligner, assembler, dataset, simple_reps, rmsk, sds)

                        for region_label, counts in svs_by_regions.items():
                            regioned_svs_counts.append((caller, asm_method, dataset, aligner, assembler, region_label, counts[0], counts[1], counts[2]))

    df_regioned_svs = pd.DataFrame(regioned_svs_counts, columns=['caller', 'asm_method', 'dataset', 'aligner', 'assembler', 'region', 'assm_unique','intersects', 'align_unique'])
    df_regioned_svs.to_csv(f'{WORKDIR}/strategy_compare_byregions.tsv', header=True, sep='\t', index=False)


def get_caller_compare_info(callers, merged_vcf, compare_outdir, read_caller, asm_caller, aligner, assembler, dataset, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    extd_matched_list = []
    extd_unique_list = []

    svs_by_regions = {repeat: [0, 0, 0] for repeat in GENOMICREGIONS}


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

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
            if dataset == 'ont_9kb':
                region_label, rptype, pcrt = annotate_sv_region_without_priTD(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if supp == len(callers):
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_caller = callers[supp_vec.index('1')]
                unique_list.append((merged_id, merged_type, unique_caller, entries[0], int(entries[1]), int(info_dict['END']), info_dict['SVLEN'], region_label, rptype, pcrt))
                if region_label != 'Tandem Repeats':
                    extd_unique_list.append((merged_id, merged_type, unique_caller, entries[0], int(entries[1]), int(info_dict['END']), info_dict['SVLEN'], region_label, rptype, pcrt))

                if supp_vec == '01':
                    svs_by_regions[region_label][2] += 1
                if supp_vec == '10':
                    svs_by_regions[region_label][0] += 1


    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.uniques.info.tsv', sep='\t', header=True, index=False)

    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.exTD.tsv',header=True, sep='\t', index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{read_caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.uniques.info.exTD.tsv', sep='\t', header=True, index=False)

    return svs_by_regions

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

def annotate_sv_region_without_priTD(chrom, start, end, pcrt_thresh, simreps_tabix, rmsk_tabix, sd_tabix):

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

    print('\nPrepare data for Figure 3a, 3b and 3c =====')
    print(f'Intermediate file directory: {WORKDIR}/aligner_dataset_fig3a3b_tmpfile')

    prepare_data()

    if not os.path.exists(f'{FIGDIR}/Fig3'):
        os.mkdir(f'{FIGDIR}/Fig3')

    print('\n==== Creating Figure3a, 3b and 3c =====')

    plot_3a3b(f'{FIGDIR}/Fig3')
    print(f'Figures saved to {FIGDIR}/Fig3/fig3a.pdf; {FIGDIR}/Fig3/fig3b.pdf')

    plot_3c(f'{FIGDIR}/Fig3')
    print(f'Figures saved to {FIGDIR}/Fig3/fig3c.pdf')

if __name__ == '__main__':
    main()