#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/2/22

'''
import os
import sys
import math
import pysam
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


def plot_s5a(figdir):

    unique_pcrt = []

    for caller in ASMCALLERS:
        wgs_supp2_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}
        extd_supp2_pcrt = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}

        for dataset_idx, dataset in enumerate(DATASETS):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            unique_info = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs5_tmpfile/{caller}.assembler-unique.exTD.tsv', header=[0], sep='\t')
            matched_info = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs5_tmpfile/{caller}.assembler-concordant.exTD.tsv', header=[0], sep='\t')

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

            extd_supp2_pcrt[plat][this_dataset_index] += supp_dict[2] / merged_total * 100

            merged_vcf = f'{WORKDIR}/{plat}/{dataset}_figs5_tmpfile/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(merged_vcf)

            unique_pcrt.append((wgs_supp_dict[1] / wgs_merged_total * 100, plat, caller, 'WGS'))
            wgs_supp2_pcrt[plat][this_dataset_index] += wgs_supp_dict[2] / wgs_merged_total * 100

    df_unique = pd.DataFrame(unique_pcrt, columns=['pcrt', 'plat', 'caller', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    fig, ax = plt.subplots(1, 1, figsize=(3, 4))
    sns.barplot(data=df_unique, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=ax)


    ax.set_ylabel('% of assembler unique SVs', fontsize=13)

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
    fig.savefig(f'{figdir}/figs5a.pdf')
    # plt.show()



def plot_s5b(figdir):

    regions = ['Repeat Masked', 'Segment Dup', 'Simple Region', 'Tandem Repeats']
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}

    svtypes_info = []

    for dataset_idx, dataset in enumerate(DATASETS):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(ASMCALLERS):
            unique_info_out = f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.assembler-unique.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])
            svtypes_dict = {assm: {'INS': 0, 'DEL': 0} for assm in plat_assemblers[plat]}
            this_total = {assm: 0 for assm in plat_assemblers[plat]}
            assembler_unique_region_dict = {assembler: [0 for i in range(len(regions))] for assembler in
                                            plat_assemblers[plat]}
            assembler_unique_dict = {assembler: 0 for assembler in plat_assemblers[plat]}

            for idx, row in df_unique.iterrows():
                assembler, svtype, svlen, sv_region = row['ASSEMBLER'], row['TYPE_MATCH'], abs(int(row['SVLEN'])), row[
                    'REGION_TYPE']
                if svtype in ['INS', 'DEL']:
                    svtypes_dict[assembler][svtype] += 1
                    this_total[assembler] += 1
                    region_idx = regions.index(sv_region)
                    assembler_unique_dict[assembler] += 1
                    assembler_unique_region_dict[assembler][region_idx] += 1

            for assm, svtype_count in svtypes_dict.items():
                for svtype, count in svtype_count.items():
                    svtypes_info.append((svtype, count / this_total[assm] * 100, assm, plat, dataset, caller))


    df_svtypes_info = pd.DataFrame(svtypes_info, columns=['svtype', 'pcrt', 'assembler', 'plat', 'dataset', 'caller'])

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    sns.boxplot(data=df_svtypes_info[df_svtypes_info['plat'] == 'HiFi'], x='svtype', y='pcrt', hue='assembler',
                hue_order=['flye', 'hifiasm'], palette=[ASSMBLERCOLOR['flye'], ASSMBLERCOLOR['hifiasm']], ax=axes[0])

    sns.boxplot(data=df_svtypes_info[df_svtypes_info['plat'] == 'ONT'], x='svtype', y='pcrt', hue='assembler',
                hue_order=['flye', 'shasta'], palette=[ASSMBLERCOLOR['flye'], '#eb5362'], ax=axes[1])

    for i, ax in enumerate(axes):
        ax.set_ylim(0, 90)
        ax.set_yticks(np.linspace(0, 90, 4))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 90, 4)], fontsize=12)
        ax.legend(title='Assembler', loc='lower left')
        if i == 0:
            ax.set_ylabel('% of assembler unique SV', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)
        ax.set_xlabel('')

    fig.tight_layout()
    # plt.show()

    fig.savefig(f'{figdir}/figs5b.pdf')

def plot_s5c(figdir):
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}
    datasets = ['hifi_18kb', 'ont_30kb']

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        aligner_unique_types = []

        for caller in ASMCALLERS:
            for assembler in plat_assemblers[plat]:
                df_uniques = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{assembler}-{caller}.aligner-unique.tsv', sep='\t')
                unique_types_dict = {'lra':{}, 'minimap2': {}}
                unique_total = {'lra': 0, 'minimap2': 0}
                for idx, row in df_uniques.iterrows():
                    svtype, aligner = row['TYPE_MATCH'], row['ALIGNER']
                    unique_total[aligner] += 1
                    if svtype in unique_types_dict[aligner]:
                        unique_types_dict[aligner][svtype] += 1
                    else:
                        unique_types_dict[aligner][svtype] = 1

                for aligner, count in unique_types_dict.items():
                    for key, val in count.items():
                        aligner_unique_types.append((aligner, key, val / unique_total[aligner] * 100, caller, assembler))


        df_aligner_unique_type = pd.DataFrame(aligner_unique_types, columns=['aligner', 'svtype', 'pcrt', 'caller', 'assembler'])

        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        sns.barplot(data=df_aligner_unique_type, x='svtype', y='pcrt', hue='aligner', hue_order=['lra', 'minimap2'],
                    palette=['#1b9e77', '#d95f02'], capsize=.2, ax=ax)

        ax.set_ylabel('% of contig aligner specific SV types', fontsize=13)
        ax.set_ylim([0, 80])
        ax.set_yticks(np.linspace(0, 80, 5))
        ax.set_yticklabels([int(ele) for ele in np.linspace(0, 80, 5)], fontsize=12)


        ax.set_xlabel('')
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['INS', 'DEL', 'INV'], fontsize=13)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

        fig.tight_layout()
        fig.savefig(f'{figdir}/figs5c-{plat}.pdf')
        # plt.show()


def prepare_data():
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        outdir = f'{WORKDIR}/{plat}/{dataset}_figs5_tmpfile'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for caller in ASMCALLERS:
            tmp_file = f'{outdir}/{caller}.{dataset}.txt'
            tmp_file_writer = open(tmp_file, 'w')
            for assembler in plat_assemblers[plat]:
                vcf_file = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                if dataset == 'ont_9kb' and caller == 'pav':
                    vcf_file = f'{WORKDIR}/{plat}/lra_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                print(vcf_file, file=tmp_file_writer)
            tmp_file_writer.close()

            merged_out_vcf = f'{outdir}/{caller}.{dataset}.jasmine.merged.vcf'
            print(f'Obtain {caller} {dataset} assembler uniques and overlaps ...')
            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {outdir}/jasmine.callers.log'
            os.system(cmd)
            # print(f'Annotating {caller} {dataset} ...')
            get_assembler_compare_info(plat_assemblers[plat], merged_out_vcf, outdir, caller, simple_reps, rmsk, sds)
            os.remove(tmp_file)


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

    df_matched = pd.DataFrame(matched_list,columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{caller}.assembler-concordant.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'ASSEMBLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{caller}.assembler-unique.tsv', sep='\t', header=True, index=False)

    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END','SVLEN', 'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{caller}.assembler-concordant.exTD.tsv', header=True, sep='\t', index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list,columns=['ID_MATCH', 'TYPE_MATCH', 'ASSEMBLER', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{caller}.assembler-unique.exTD.tsv', sep='\t', header=True, index=False)

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

def main():
    print('\nPrepare data for Extended Data Fig 5 =====')
    print(f'Intermediate file directory: {WORKDIR}/dataset_figs5_tmpfile')

    prepare_data()

    print('\n==== Creating Extended Data Fig 5 =====')

    if not os.path.exists(f'{FIGDIR}/FigS5'):
        os.mkdir(f'{FIGDIR}/FigS5')

    plot_s5a(f'{FIGDIR}/FigS5')
    print(f'Figures saved to {FIGDIR}/FigS5/figs5a.pdf')

    plot_s5b(f'{FIGDIR}/FigS5')
    print(f'Figures saved to {FIGDIR}/FigS5/figs5b.pdf')

    plot_s5c(f'{FIGDIR}/FigS5')
    print(f'Figures saved to {FIGDIR}/FigS5/figs5c.pdf')

if __name__ == '__main__':
    main()