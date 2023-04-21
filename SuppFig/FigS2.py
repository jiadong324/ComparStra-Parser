#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/15

'''
import pandas as pd
import math
import os
import sys
import pysam
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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


def prepare_data():

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    plat_assembler = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    # dataset_dict = {'ONT': ['ont_9kb', 'ont_19kb', 'ont_30kb']}

    for platform, datasets in DATASET_DICT.items():

        merged_outdir = f'{WORKDIR}/{platform}/figs2_tmpfile'
        if not os.path.exists(merged_outdir):
            os.mkdir(merged_outdir)

        for caller in READCALLERS:
            for aligner in READALIGNERS:

                tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:
                    vcf_file = f'{WORKDIR}/{platform}/{aligner}_{dataset}/filtered/HG002.{caller}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{WORKDIR}/{platform}/figs2_tmpfile/{caller}.{aligner}.jasmine.merged.vcf'
                print(f'Producing {caller} {aligner} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {merged_outdir}/read_jasmine.log'
                os.system(cmd)

                os.remove(tmp_file)
                get_dataset_compare_info(datasets, merged_out_vcf, merged_outdir, caller, aligner, platform, simple_reps, rmsk, sds)


        for caller in ASMCALLERS:
            for assembler in plat_assembler[platform]:

                tmp_file = f'{merged_outdir}/{caller}.{assembler}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:
                    vcf_file = f'{WORKDIR}/{platform}/minimap2_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{merged_outdir}/{caller}.{assembler}.jasmine.merged.vcf'
                print(f'Producing {caller} {assembler} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {merged_outdir}/assm_jasmine.log'
                os.system(cmd)

                os.remove(tmp_file)

                get_dataset_compare_info(datasets, merged_out_vcf, merged_outdir, caller, assembler, platform, simple_reps, rmsk, sds)

def plot_s2a(figdir):

    plat_assembler = {'HiFi':['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    plats = ['HiFi', 'ONT']

    all_supp_pcrt = {'HiFi': [], 'ONT': []}

    for fig_idx, plat in enumerate(plats):

        for caller in READCALLERS:
            for aligner in READALIGNERS:
                merged_vcf = f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{aligner}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                for supp, count in supp_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, aligner, caller, 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                merged_vcf = f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{assembler}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)

                for supp, count in supp_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, assembler, caller, 'Assembly'))

    stra_order = ['Assembly', 'Read']
    stra_color = [STRACOLORS[ele] for ele in stra_order]
    stra_legends = [Patch(label=ele, color=STRACOLORS[ele]) for ele in stra_order]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 3))
    for col_idx, plat in enumerate(plats):
        this_ax = axes[col_idx]
        df_all_pcrt = pd.DataFrame(all_supp_pcrt[plat], columns=['pcrt', 'supp', 'aa', 'caller', 'stra'])
        sns.barplot(data=df_all_pcrt, x='supp', y='pcrt', hue='stra', hue_order=stra_order, palette=stra_color, capsize=.15, ax=this_ax)
        this_ax.set_title(plat, fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of WGS-SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.legend(handles=stra_legends)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        ax.set_xlabel('# of datasets', fontsize=13)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs2a-1.pdf')

    all_supp_pcrt = {'HiFi': [], 'ONT': []}

    for fig_idx, plat in enumerate(plats):

        for caller in READCALLERS:
            for aligner in READALIGNERS:
                matched_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{aligner}.{plat}-concordant.info.exTD.tsv', header=[0], sep='\t')
                unique_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{aligner}.{plat}.unique.exTD.tsv', header=[0], sep='\t')

                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1

                for supp, count in supps_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, aligner, caller, 'Read'))
                all_supp_pcrt[plat].append((len(unique_info) / merged_total * 100, 1, aligner, caller, 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                matched_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{assembler}.{plat}-concordant.info.exTD.tsv',header=[0], sep='\t')
                unique_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{assembler}.{plat}.unique.exTD.tsv', header=[0],sep='\t')
                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1

                for supp, count in supps_dict.items():
                    all_supp_pcrt[plat].append((count / merged_total * 100, supp, assembler, caller, 'Assembly'))
                all_supp_pcrt[plat].append((len(unique_info) / merged_total * 100, 1, assembler, caller, 'Assembly'))

    stra_order = ['Assembly', 'Read']
    stra_color = [STRACOLORS[ele] for ele in stra_order]
    stra_legends = [Patch(label=ele, color=STRACOLORS[ele]) for ele in stra_order]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 3))
    for col_idx, plat in enumerate(plats):
        this_ax = axes[col_idx]
        df_all_pcrt = pd.DataFrame(all_supp_pcrt[plat], columns=['pcrt', 'supp', 'aa', 'caller', 'stra'])
        sns.barplot(data=df_all_pcrt, x='supp', y='pcrt', hue='stra', hue_order=stra_order,
                    palette=stra_color, capsize=.15, ax=this_ax)

        this_ax.set_title(plat, fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of ExTD-SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.legend(handles=stra_legends)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        ax.set_xlabel('# of datasets', fontsize=13)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs2a-2.pdf')
    # plt.show()

def plot_s2b(figdir):

    plat_assembler = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    plats = ['HiFi', 'ONT']

    wgs_repro_pcrt = []
    extd_repro_pcrt = []

    for fig_idx, plat in enumerate(plats):

        for caller in READCALLERS:
            for aligner in READALIGNERS:
                merged_vcf = f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{aligner}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                wgs_repro_pcrt.append((supp_dict[3] / merged_total * 100, aligner, TOOLMAP[caller], plat, 'Read'))

                matched_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{aligner}.{plat}-concordant.info.exTD.tsv', header=[0], sep='\t')
                unique_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{aligner}.{plat}.unique.exTD.tsv', header=[0], sep='\t')

                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1
                extd_repro_pcrt.append((supps_dict[3] / merged_total * 100, aligner, TOOLMAP[caller], plat, 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                merged_vcf = f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{assembler}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)
                wgs_repro_pcrt.append((supp_dict[3] / merged_total * 100, assembler, TOOLMAP[caller], plat, 'Assembly'))

                matched_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{assembler}.{plat}-concordant.info.exTD.tsv',header=[0], sep='\t')
                unique_info = pd.read_csv(f'{WORKDIR}/{plat}/figs2_tmpfile/{caller}.{assembler}.{plat}.unique.exTD.tsv', header=[0],sep='\t')
                merged_total = len(matched_info) + len(unique_info)
                supps_dict = {}
                for idx, row in matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in supps_dict:
                        supps_dict[supp] += 1
                    else:
                        supps_dict[supp] = 1
                extd_repro_pcrt.append((supps_dict[3] / merged_total * 100, assembler, TOOLMAP[caller], plat, 'Assembly'))

    df_wgs_repro = pd.DataFrame(wgs_repro_pcrt, columns=['pcrt', 'aa', 'caller', 'plat', 'stra'])
    df_extd_repro = pd.DataFrame(extd_repro_pcrt, columns=['pcrt', 'aa', 'caller', 'plat', 'stra'])

    read_wgs_repro, read_extd_repro = df_wgs_repro[df_wgs_repro['stra'] == 'Read'], df_extd_repro[df_extd_repro['stra'] == 'Read']
    assm_wgs_repro, assm_extd_repro = df_wgs_repro[df_wgs_repro['stra'] == 'Assembly'], df_extd_repro[df_extd_repro['stra'] == 'Assembly']

    assm_legends = [Line2D([0], [0], label='SVIM-asm', color=TOOLCOLORS['svimasm'], lw=2),
                   Line2D([0], [0], label='PAV', color=TOOLCOLORS['pav'], lw=2),
                   Line2D([0], [0], label='WGS', color='black', lw=2),
                   Line2D([0], [0], label='ExTD', color='black', ls='--', lw=2)]

    fig_name = 1
    for idx, plat in enumerate(plats):
        fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(4, 4))
        assm_plat_wgs, assm_plat_extd = assm_wgs_repro[assm_wgs_repro['plat'] == plat], assm_extd_repro[assm_extd_repro['plat'] == plat]

        sns.lineplot(data=assm_plat_wgs, x='aa', y='pcrt', hue='caller', hue_order=['PAV', 'SVIM-asm'], palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']],
                     lw=2, marker='o', markersize=8, ax=axes[0])

        sns.lineplot(data=assm_plat_extd, x='aa', y='pcrt', hue='caller', hue_order=['PAV', 'SVIM-asm'], palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']],
                     lw=2, marker='o', ls='--', markersize=8, ax=axes[1])

        for i, ax in enumerate(axes):

            if i == 0:
                ax.set_ylabel('% of datasets concordant SVs', fontsize=13)
            else:
                ax.set_ylabel('')
            ax.set_xlabel('')

            ax.set_xticks([0, 1])
            if plat == 'HiFi':
                ax.set_xticklabels(['hifiasm', 'flye'], fontsize=13, rotation=60)
            else:
                ax.set_xticklabels(['flye', 'shasta'], fontsize=13, rotation=60)
            ax.set_ylim(20, 100)
            ax.set_yticks(np.linspace(20, 100, 5))
            ax.set_yticklabels([int(val) for val in np.linspace(20, 100, 5)], fontsize=12)

            if i == 0:
                ax.legend(handles=assm_legends)
            else:
                ax.legend('', frameon=False)

            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        fig.tight_layout()
        fig.savefig(f'{figdir}/figs2b-{fig_name}.pdf')
        fig_name += 1

    for plat in plats:
        fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))

        sns.lineplot(data=read_wgs_repro[read_wgs_repro['plat'] == plat], x='aa', y='pcrt', hue='caller',
                    hue_order=[TOOLMAP[ele] for ele in READCALLERS], palette=[TOOLCOLORS[ele] for ele in READCALLERS],
                    lw=2, marker='o', markersize=8, ax=axes[0])

        sns.lineplot(data=read_extd_repro[read_extd_repro['plat'] == plat], x='aa', y='pcrt', hue='caller',
                    hue_order=[TOOLMAP[ele] for ele in READCALLERS], palette=[TOOLCOLORS[ele] for ele in READCALLERS],
                    lw=2, marker='o', ls='--', markersize=8, ax=axes[1])

        for idx, ax in enumerate(axes):
            if idx == 0:
                ax.set_ylabel('% of datasets concordant SVs', fontsize=13)
            else:
                ax.set_ylabel('')
            ax.set_xlabel('')

            ax.legend('', frameon=False)

            ax.set_xticks(np.arange(len(READALIGNERS)))
            ax.set_xticklabels(READALIGNERS, rotation=60, fontsize=13)

            ax.set_ylim(20, 100)
            ax.set_yticks(np.linspace(20, 100, 5))
            ax.set_yticklabels([int(val) for val in np.linspace(20, 100, 5)], fontsize=12)

            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        fig.tight_layout()
        fig.savefig(f'{figdir}/figs2b-{fig_name}.pdf')
        fig_name += 1

    # plt.show()


def get_dataset_compare_info(datasets, merged_vcf, compare_outdir, caller, tool, platform, simple_reps, rmsk, sds):

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

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if supp > 1:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_datasets = datasets[supp_vec.index('1')]
                unique_list.append((merged_id, merged_type, unique_datasets, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))
                if region_label != 'Tandem Repeats':
                    extd_unique_list.append((merged_id, merged_type, unique_datasets, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

    ## Comparing SV at whole genome scale among datasets
    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN', 'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}-concordant.info.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns =['ID_MATCH', 'TYPE_MATCH', 'DATASET', '#CHROM', 'POS', 'END', 'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}.unique.tsv', sep='\t', header=True, index=False)

    ## Comparing SVs outside of tandem repeat regions among datasets
    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_matched.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}-concordant.info.exTD.tsv', header=True, sep='\t',index=False)

    df_extd_uniques = pd.DataFrame(extd_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'DATASET', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_extd_uniques.to_csv(f'{compare_outdir}/{caller}.{tool}.{platform}.unique.exTD.tsv', sep='\t', header=True, index=False)


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

    sorted_annotations = sorted(annotations, key=lambda x: x[1], reverse=True)

    return sorted_annotations[0]


def main():
    print('\n==== Preparing data for Extended Data Fig 2 =====')

    print(f'Intermediate file directory: {WORKDIR}/platform/figs2_tmpfile')

    prepare_data()

    if not os.path.exists(f'{FIGDIR}/FigS2'):
        os.mkdir(f'{FIGDIR}/FigS2')

    print('\n==== Creating Extended Data Fig 2a and 2b =====')

    plot_s2a(f'{FIGDIR}/FigS2')
    print(f'Figures saved to {FIGDIR}/FigS2/figs2a-1.pdf; {FIGDIR}/FigS2/figs2a-2.pdf')

    plot_s2b(f'{FIGDIR}/FigS2')
    print(f'Figures saved to {FIGDIR}/FigS2/figs2b-1.pdf; {FIGDIR}/FigS2/figs2b-2.pdf; {FIGDIR}/FigS2/figs2b-3.pdf; {FIGDIR}/FigS2/figs2b-4.pdf')

if __name__ == '__main__':
    main()
