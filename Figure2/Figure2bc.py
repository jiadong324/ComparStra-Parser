#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/16

'''


import math
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pysam
from matplotlib.patches import Patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Constant import *
from Helpers.Functions import get_overlaps, get_survivor_supp

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

    for platform, datasets in DATASET_DICT.items():

        merged_outdir = f'{WORKDIR}/{platform}/fig2b2c_tmpfile'
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

                merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
                print(f'Producing {caller} {aligner} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {merged_outdir}/read_dataset_repro.jasmine.log'
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
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {merged_outdir}/assm_dataset_repro.jasmine.log'
                os.system(cmd)
                os.remove(tmp_file)

                get_dataset_compare_info(datasets, merged_out_vcf, merged_outdir, caller, assembler, platform, simple_reps, rmsk, sds)


def plot_2b2c(figdir):

    plat_assembler = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}
    supp_pcrt = []
    plats = ['HiFi', 'ONT']

    for fig_idx, plat in enumerate(plats):
        for caller in READCALLERS:
            for aligner in READALIGNERS:
                wgs_merged_vcf = f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{aligner}.jasmine.merged.vcf'
                wgs_supp_dict, wgs_merged_total = get_survivor_supp(wgs_merged_vcf)
                supp_pcrt.append((wgs_supp_dict[3] / wgs_merged_total * 100, aligner, TOOLMAP[caller], f'{plat}-WGS', 'Read'))

                extd_matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{aligner}.{plat}-concordant.info.exTD.tsv',header=[0], sep='\t')
                extd_unique_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{aligner}.{plat}.unique.exTD.tsv', header=[0],sep='\t')

                extd_merged_total = len(extd_matched_info) + len(extd_unique_info)
                extd_supps_dict = {}
                for idx, row in extd_matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in extd_supps_dict:
                        extd_supps_dict[supp] += 1
                    else:
                        extd_supps_dict[supp] = 1

                supp_pcrt.append((extd_supps_dict[3] / extd_merged_total * 100, aligner, TOOLMAP[caller], f'{plat}-ExTD', 'Read'))

        for col_idx, assembler in enumerate(plat_assembler[plat]):
            for caller in ASMCALLERS:
                wgs_merged_vcf = f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{assembler}.jasmine.merged.vcf'
                wgs_supp_dict, wgs_merged_total = get_survivor_supp(wgs_merged_vcf)
                supp_pcrt.append((wgs_supp_dict[3] / wgs_merged_total * 100, assembler, TOOLMAP[caller], f'{plat}-WGS', 'Assembly'))

                extd_matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{assembler}.{plat}-concordant.info.exTD.tsv',header=[0], sep='\t')
                extd_unique_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{assembler}.{plat}.unique.exTD.tsv', header=[0],sep='\t')
                extd_merged_total = len(extd_matched_info) + len(extd_unique_info)
                extd_supps_dict = {}
                for idx, row in extd_matched_info.iterrows():
                    supp = int(row['SUPP'])
                    if supp in extd_supps_dict:
                        extd_supps_dict[supp] += 1
                    else:
                        extd_supps_dict[supp] = 1

                supp_pcrt.append((extd_supps_dict[3] / extd_merged_total * 100, assembler, TOOLMAP[caller], f'{plat}-ExTD', 'Assembly'))

    stra_order = ['Assembly', 'Read']
    stra_color = [STRACOLORS[ele] for ele in stra_order]
    stra_legends = [Patch(label=ele, color=STRACOLORS[ele]) for ele in stra_order]

    df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'aa', 'caller', 'plat', 'stra'])
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    sns.barplot(data=df_supp_pcrt, x='plat', y='pcrt', hue='stra',
                hue_order=stra_order, palette=stra_color, capsize=0.15, ax=ax)

    ax.set_ylabel('% of concordant SVs', fontsize=13)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    ax.legend(handles=stra_legends)
    ax.set_xlabel('')
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(['WGS-SVs', 'ExTD-SVs', 'WGS-SVs', 'ExTD-SVs'], fontsize=13)

    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig2b.pdf')

    shift_labels = ['0', '0,10', '>10']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']
    bpstd = []

    for k, plat in enumerate(['HiFi', 'ONT']):

        wgs_start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
        extd_start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}

        wgs_caller_total = [1 for k in range(len(callers))]
        extd_caller_total = [1 for k in range(len(callers))]

        for asm_caller in ASMCALLERS:
            caller_index = callers.index(asm_caller)
            for assembler in plat_assembler[plat]:

                wgs_matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{asm_caller}.{assembler}.{plat}-concordant.info.tsv', sep='\t', header=[0])
                for idx, row in wgs_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        wgs_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            wgs_start_std_dict['0'][caller_index] += 1

                extd_matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{asm_caller}.{assembler}.{plat}-concordant.info.exTD.tsv', sep='\t', header=[0])
                for idx, row in extd_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        extd_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            extd_start_std_dict['0'][caller_index] += 1

        for caller in READCALLERS:
            caller_index = callers.index(caller)
            for aligner in READALIGNERS:
                wgs_matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{aligner}.{plat}-concordant.info.tsv', sep='\t', header=[0])
                for idx, row in wgs_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        wgs_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            wgs_start_std_dict['0'][caller_index] += 1

                extd_matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{aligner}.{plat}-concordant.info.exTD.tsv', sep='\t', header=[0])
                for idx, row in extd_matched_info.iterrows():
                    if int(row['SUPP']) == 3:
                        extd_caller_total[caller_index] += 1
                        if abs(float(row['START_STD'])) == 0:
                            extd_start_std_dict['0'][caller_index] += 1

        extd_start_shift = [i / j * 100 for i, j in zip(extd_start_std_dict['0'], extd_caller_total)]
        for i, pcrt in enumerate(extd_start_shift):
            if i < len(READCALLERS):
                bpstd.append((pcrt, plat, callers[i], 'Read', 0, 'ExTD'))
            else:
                bpstd.append((pcrt, plat, callers[i], 'Assembly', 0, 'ExTD'))

        wgs_start_shift = [i / j * 100 for i, j in zip(wgs_start_std_dict['0'], wgs_caller_total)]

        for i, pcrt in enumerate(wgs_start_shift):
            if i < len(READCALLERS):
                bpstd.append((pcrt, plat, callers[i], 'Read', 0, 'WGS'))
            else:
                bpstd.append((pcrt, plat, callers[i], 'Assembly', 0, 'WGS'))

    df_bpstd = pd.DataFrame(bpstd, columns=['pcrt', 'plat', 'caller', 'stra', 'shift', 'region'])
    # df_extd_bpstd = pd.DataFrame(extd_bpstd, columns=['pcrt', 'plat', 'caller', 'shift', 'stra'])
    hifi_bpstd = df_bpstd[df_bpstd['plat'] == 'HiFi']
    ont_bpstd = df_bpstd[df_bpstd['plat'] == 'ONT']

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(4, 4))

    sns.boxplot(data=hifi_bpstd, x='region', y='pcrt', hue='stra', hue_order=['Read', 'Assembly'],
                palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], ax=axes[0])
    axes[0].set_title('HiFi', fontsize=13)

    sns.boxplot(data=ont_bpstd, x='region', y='pcrt', hue='stra', hue_order=['Read', 'Assembly'],
                palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], ax=axes[1])
    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of BSD-0 SVs', fontsize=13)
            ax.legend(title='')
        else:
            ax.set_ylabel('')
            ax.legend('', frameon=False)

        ax.set_xlabel('')
        # ax.set_xticks([0, 1])
        # ax.set_xticklabels(['WGS', 'ExTD'], fontsize=13)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig2c-1.pdf')


    xaxis_labels = ['pav', 'svimasm', 'pbsv', 'sniffles', 'svision', 'svim', 'cutesv']
    barwidth = 0.35

    hifi_xticks = np.arange(len(xaxis_labels))
    ont_xticks = [r + barwidth + 0.06 for r in hifi_xticks]
    rs = [hifi_xticks, ont_xticks]
    xticks = [r + (barwidth + 0.06) / 2 for r in hifi_xticks]


    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    for plat_idx, plat in enumerate(['HiFi', 'ONT']):
        shift_dict = {shift: [0 for k in range(len(xaxis_labels))] for shift in shift_labels}
        totals = [0 for k in range(len(xaxis_labels))]

        for caller in READCALLERS:
            caller_index = xaxis_labels.index(caller)
            for row_idx, aligner in enumerate(READALIGNERS):
                matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{aligner}.{plat}-concordant.info.tsv', sep='\t',
                    header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region, rep_pcrt = abs(float(row['START_STD'])), abs(
                        float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE'], float(
                        row['PCRT'])

                    totals[caller_index] += 1

                    if start_std == 0:
                        shift_dict['0'][caller_index] += 1
                    elif start_std <= 10:
                        shift_dict['0,10'][caller_index] += 1
                    else:
                        shift_dict['>10'][caller_index] += 1

        assemblers = plat_assembler[plat]
        for caller in ASMCALLERS:
            caller_index = xaxis_labels.index(caller)
            for assembler in assemblers:
                matched_info = pd.read_csv(f'{WORKDIR}/{plat}/fig2b2c_tmpfile/{caller}.{assembler}.{plat}-concordant.info.tsv', sep='\t',header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])),row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                    totals[caller_index] += 1

                    if start_std == 0:
                        shift_dict['0'][caller_index] += 1
                    elif start_std <= 10:
                        shift_dict['0,10'][caller_index] += 1
                    else:
                        shift_dict['>10'][caller_index] += 1

        for l, shift in enumerate(shift_labels):

            this_shift = [i / j * 100 for i, j in zip(shift_dict[shift], totals)]

            if l == 0:
                if plat == 'HiFi':
                    ax.bar(rs[plat_idx], this_shift, color=SHIFTCOLORDICT[shift_labels[l]], width=barwidth,
                           edgecolor='black', label=shift)
                else:
                    ax.bar(rs[plat_idx], this_shift, color=SHIFTCOLORDICT[shift_labels[l]], width=barwidth,
                           edgecolor='black', label=shift, hatch='\\\\')
            else:
                bottom = []
                for m in range(0, l):
                    bottom.append([i / j * 100 for i, j in zip(shift_dict[shift_labels[m]], totals)])
                bottom_sum = [sum(x) for x in zip(*bottom)]
                if plat == 'HiFi':
                    ax.bar(rs[plat_idx], this_shift, bottom=bottom_sum, color=SHIFTCOLORDICT[shift_labels[l]],
                           width=barwidth, edgecolor='black', label=shift)
                else:
                    ax.bar(rs[plat_idx], this_shift, bottom=bottom_sum, color=SHIFTCOLORDICT[shift_labels[l]],
                           width=barwidth, edgecolor='black', label=shift, hatch='\\\\')

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=2)

    ax.set_xlabel('')
    ax.set_xticks(xticks)
    ax.set_xticklabels([TOOLMAP[caller] for caller in xaxis_labels], fontsize=13)

    ax.legend('', frameon=False)

    ax.set_ylabel('% of concordant SVs', fontsize=13)
    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig2c-2.pdf')

    plt.show()

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
    print('\nPrepare data for Figure 2b and 2c =====')
    print(f'Intermediate file directory: {WORKDIR}/fig2b2c_tmpfile')

    prepare_data()

    if not os.path.exists(f'{FIGDIR}/Fig2'):
        os.mkdir(f'{FIGDIR}/Fig2')

    print('\n==== Creating Figure2b and 2c =====')
    plot_2b2c(f'{FIGDIR}/Fig2')

    print(f'Figures saved to {FIGDIR}/Fig2/fig2b.pdf; {FIGDIR}/Fig2/fig2c.pdf;')

if __name__ == '__main__':
    main()