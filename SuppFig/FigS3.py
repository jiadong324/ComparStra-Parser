#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/18

'''

import sys
import os
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
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    for dataset in ['hifi_18kb', 'ont_30kb']:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        compare_dir = f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile'
        if not os.path.exists(compare_dir):
            os.mkdir(compare_dir)

        pw_read_aligner_compare(compare_dir, READALIGNERS, plat, dataset, simple_reps, rmsk, sds)
        compare_between_ctg_aligner(compare_dir, plat, dataset, plat_assemblers, simple_reps, rmsk, sds)
        compare_between_ctg_assembler(compare_dir, plat, dataset, plat_assemblers, simple_reps, rmsk, sds)

def plot_s3a(figdir):
    datasets = ['hifi_18kb', 'ont_30kb']
    plat_assemblers = {'HiFi': ['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}
    matched_pcrt = []

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in ASMCALLERS:
            for aligner in ['minimap2', 'lra']:
                df_matched = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{aligner}-{caller}.assembler-concordant.tsv',sep='\t')
                df_unique = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{aligner}-{caller}.assembler-unique.tsv', sep='\t')
                assm = plat_assemblers[plat][0] + '-' + plat_assemblers[plat][1]
                matched_pcrt.append((caller, f'{assm} ({plat})', plat, 100 * len(df_matched) / (len(df_unique) + len(df_matched))))

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for caller in READCALLERS:
            for i in range(len(READALIGNERS)):
                for j in range(len(READALIGNERS)):
                    if j <= i:
                        continue
                    df_matched = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{caller}/{caller}.{READALIGNERS[i]}-{READALIGNERS[j]}.concordant.tsv', sep='\t')
                    df_unique = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{caller}/{caller}.{READALIGNERS[i]}-{READALIGNERS[j]}.unique.tsv', sep='\t')
                    matched_pcrt.append((caller, f'{READALIGNERS[i]}-{READALIGNERS[j]}', plat, 100 * len(df_matched) / (len(df_unique) + len(df_matched))))

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    df_matched_pcrt = pd.DataFrame(matched_pcrt, columns=['caller', 'aligner', 'plat', 'pcrt'])

    my_pal = {aligner: '#fc8d62' if 'flye' in aligner else '#8da0cb' for aligner in df_matched_pcrt['aligner'].unique()}
    sns.violinplot(data=df_matched_pcrt, x='pcrt', y='aligner', palette=my_pal, scale='width', ax=ax)

    ax.set_xlabel('Percent of concordant SVs', fontsize=13)
    ax.set_xlim([25, 100])
    ax.set_xticks(np.linspace(25, 100, 4))
    ax.set_xticklabels([int(ele) for ele in np.linspace(25, 100, 4)], fontsize=12)
    ax.set_ylabel('')

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs3a-1.pdf')

    matched_pcrt = []
    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in ASMCALLERS:
            for assembler in plat_assemblers[plat]:
                df_matched = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{assembler}-{caller}.aligner-concordant.tsv', sep='\t')
                df_unique = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{assembler}-{caller}.aligner-unique.tsv', sep='\t')
                matched_pcrt.append((caller, 'ctg-minimap2-lra', 100 * len(df_matched) / (len(df_matched) + len(df_unique))))

                if plat == 'ONT':
                    print(caller, assembler, 100 * len(df_matched) / (len(df_matched) + len(df_unique)))

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in READCALLERS:
            for i in range(len(READALIGNERS)):
                for j in range(len(READALIGNERS)):
                    if j <= i:
                        continue
                    df_matched = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{caller}/{caller}.{READALIGNERS[i]}-{READALIGNERS[j]}.concordant.tsv',sep='\t')
                    df_unique = pd.read_csv(f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{caller}/{caller}.{READALIGNERS[i]}-{READALIGNERS[j]}.unique.tsv',sep='\t')
                    matched_pcrt.append((caller, f'{READALIGNERS[i]}-{READALIGNERS[j]}', 100 * len(df_matched) / (len(df_unique) + len(df_matched))))

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    df_matched_pcrt = pd.DataFrame(matched_pcrt, columns=['caller', 'aligner', 'pcrt'])

    my_pal = {aligner: '#fc8d62' if 'ctg' in aligner else '#8da0cb' for aligner in
              df_matched_pcrt['aligner'].unique()}
    sns.violinplot(data=df_matched_pcrt, x='pcrt', y='aligner', palette=my_pal, scale='width', ax=ax)

    ax.set_xlabel('Percent of concordant SVs', fontsize=13)
    ax.set_xlim([25, 100])
    ax.set_xticks(np.linspace(25, 100, 4))
    ax.set_xticklabels([int(ele) for ele in np.linspace(25, 100, 4)], fontsize=12)
    ax.set_ylabel('')

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/figs3a-2.pdf')
    # plt.show()


def pw_read_aligner_compare(compare_dir, aligners, plat, dataset, simple_reps, rmsk, sds):

    for caller in READCALLERS:
        outdir = f'{compare_dir}/{caller}'
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for i in range(len(aligners)):
            for j in range(len(aligners)):
                if i == j:
                    continue
                aligner_unique, aligner_concor = [], []
                tmp_file = f'{outdir}/{caller}.{aligners[i]}-{aligners[j]}.txt'
                tmp_file_writer = open(tmp_file, 'w')
                vcf_file1 = f'{WORKDIR}/{plat}/{aligners[i]}_{dataset}/filtered/HG002.{caller}.vcf'
                vcf_file2 = f'{WORKDIR}/{plat}/{aligners[j]}_{dataset}/filtered/HG002.{caller}.vcf'

                print(f'{vcf_file1}', file=tmp_file_writer)
                print(f'{vcf_file2}', file=tmp_file_writer)

                tmp_file_writer.close()
                merged_out_vcf = f'{outdir}/{caller}.{aligners[i]}-{aligners[j]}.jasmine.merged.vcf'
                print(f'Producing {caller} {aligners[i]}-{aligners[j]} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 --dup_to_ins ' \
                      f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1 > {outdir}/jasmine.pw.aligners.log'

                os.system(cmd)
                os.remove(tmp_file)
                # print(f'Processing {dataset} {caller} {aligners[i]}-{aligners[j]} ...')
                with open(merged_out_vcf, 'r') as f:
                    for line in f:
                        if '#' in line:
                            continue

                        entries = line.strip().split('\t')
                        info_tokens = entries[7].split(';')
                        merged_id = entries[2]
                        info_dict = {}

                        for token in info_tokens:
                            if '=' in token:
                                info_dict[token.split('=')[0]] = token.split('=')[1]

                        end, merged_type = int(info_dict['END']), info_dict['SVTYPE']

                        region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 50, simple_reps, rmsk, sds)

                        supp, supp_vec = int(info_dict['SUPP']), info_dict['SUPP_VEC']
                        if supp == 2:
                            aligner_concor.append((merged_id, merged_type, f'{aligners[i]}-{aligners[j]}', entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))
                        else:
                            index = supp_vec.index('1')
                            this_aligners = [aligners[i], aligners[j]]
                            aligner_unique.append((merged_id, merged_type, this_aligners[index], entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

                df_matched = pd.DataFrame(aligner_concor, columns=['ID_MATCH', 'TYPE_MATCH', 'ALIGNER', '#CHROM', 'POS', 'END', 'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])

                df_matched.to_csv(f'{outdir}/{caller}.{aligners[i]}-{aligners[j]}.concordant.tsv', header=True, sep='\t', index=False)

                df_uniques = pd.DataFrame(aligner_unique, columns=['ID_MATCH', 'TYPE_MATCH', 'ALIGNER', '#CHROM', 'POS', 'END',
                                                   'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])

                df_uniques.to_csv(f'{outdir}/{caller}.{aligners[i]}-{aligners[j]}.unique.tsv', sep='\t', header=True, index=False)


def compare_between_ctg_aligner(compare_dir,  plat, dataset, plat_assemblers, simple_reps, rmsk, sds):

    aligners = ['lra', 'minimap2']

    for caller in ASMCALLERS:
        for assembler in plat_assemblers[plat]:
            aligner_unique, aligner_concor = [], []

            merge_txt = f'{compare_dir}/{assembler}-{caller}.aligners.merged.txt'
            merge_txt_writer = open(merge_txt, 'w')

            for aligner in aligners:
                vcf_file = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                print(vcf_file, file=merge_txt_writer)
            merge_txt_writer.close()

            merged_out_vcf = f'{compare_dir}/{assembler}-{caller}.{dataset}.aligners.merged.vcf'
            print(f'Producing {caller} {assembler} aligners merged calls ...')
            cmd = f'{JASMINE} file_list={merge_txt} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {compare_dir}/jasmine.ctg_aligner.log'
            os.system(cmd)
            os.remove(merge_txt)

            with open(merged_out_vcf, 'r') as f:
                for line in f:
                    if '#' in line:
                        continue

                    entries = line.strip().split('\t')
                    info_tokens = entries[7].split(';')
                    merged_id = entries[2]
                    info_dict = {}

                    for token in info_tokens:
                        if '=' in token:
                            info_dict[token.split('=')[0]] = token.split('=')[1]

                    end, merged_type = int(info_dict['END']), info_dict['SVTYPE']
                    supp, supp_vec = int(info_dict['SUPP']), info_dict['SUPP_VEC']
                    region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 50, simple_reps, rmsk, sds)

                    if supp == 2:
                        aligner_concor.append((merged_id, merged_type, caller, assembler, 'lra-minimap2', entries[0],
                                               int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))
                    else:
                        index = supp_vec.index('1')
                        aligner_unique.append((merged_id, merged_type, caller, assembler, aligners[index], entries[0], int(entries[1]),
                                               end, info_dict['SVLEN'], region_label, rptype, pcrt))

            df_matched = pd.DataFrame(aligner_concor, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', 'ASSEMBLER', 'ALIGNER', '#CHROM', 'POS', 'END',
                                                               'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])

            df_matched.to_csv(f'{compare_dir}/{assembler}-{caller}.aligner-concordant.tsv', header=True, sep='\t', index=False)

            df_uniques = pd.DataFrame(aligner_unique, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', 'ASSEMBLER', 'ALIGNER', '#CHROM', 'POS', 'END', 'SVLEN',
                                               'REGION_TYPE', 'RPTYPE', 'PCRT'])

            df_uniques.to_csv(f'{compare_dir}/{assembler}-{caller}.aligner-unique.tsv', sep='\t', header=True, index=False)

def compare_between_ctg_assembler(compare_dir, plat, dataset, plat_assemblers, simple_reps, rmsk, sds):

    for caller in ASMCALLERS:
        for aligner in ['lra', 'minimap2']:

            aligner_unique, aligner_concor = [], []

            merge_txt = f'{compare_dir}/{aligner}-{caller}.assemblers.merged.txt'
            merge_txt_writer = open(merge_txt, 'w')

            for assembler in plat_assemblers[plat]:
                vcf_file = f'{WORKDIR}/{plat}/{aligner}_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                print(vcf_file, file=merge_txt_writer)
            merge_txt_writer.close()

            merged_out_vcf = f'{compare_dir}/{aligner}-{caller}.{dataset}.assemblers.merged.vcf'
            print(f'Producing {caller} {aligner} assemblers merged calls ...')
            cmd = f'{JASMINE} file_list={merge_txt} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1 > {compare_dir}/jasmine.ctg_assembler.log'
            os.system(cmd)
            os.remove(merge_txt)

            with open(merged_out_vcf, 'r') as f:
                for line in f:
                    if '#' in line:
                        continue

                    entries = line.strip().split('\t')
                    info_tokens = entries[7].split(';')
                    merged_id = entries[2]
                    info_dict = {}

                    for token in info_tokens:
                        if '=' in token:
                            info_dict[token.split('=')[0]] = token.split('=')[1]

                    end, merged_type = int(info_dict['END']), info_dict['SVTYPE']
                    supp, supp_vec = int(info_dict['SUPP']), info_dict['SUPP_VEC']
                    region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 50, simple_reps, rmsk, sds)

                    if supp == 2:
                        aligner_concor.append((merged_id, merged_type, caller, 'flye-hifiasm', aligner, entries[0],
                                               int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))
                    else:
                        index = supp_vec.index('1')
                        aligner_unique.append((merged_id, merged_type, caller, plat_assemblers[plat][index], aligner, entries[0], int(entries[1]),
                                               end, info_dict['SVLEN'], region_label, rptype, pcrt))

            df_matched = pd.DataFrame(aligner_concor, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', 'ASSEMBLER', 'ALIGNER', '#CHROM', 'POS', 'END',
                                                               'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])

            df_matched.to_csv(f'{compare_dir}/{aligner}-{caller}.assembler-concordant.tsv', header=True, sep='\t', index=False)

            df_uniques = pd.DataFrame(aligner_unique, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', 'ASSEMBLER', 'ALIGNER', '#CHROM', 'POS', 'END', 'SVLEN',
                                               'REGION_TYPE', 'RPTYPE', 'PCRT'])

            df_uniques.to_csv(f'{compare_dir}/{aligner}-{caller}.assembler-unique.tsv', sep='\t', header=True, index=False)


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


def plot_s3b(figdir):

    unique_pcrt = []

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

            unique_pcrt.append((len(unique_info) / merged_total * 100, caller, plat, 'ExTD'))
            extd_supp4_pcrt[plat][this_dataset_index] += supp_dict[4] / merged_total * 100

            wgs_merged_vcf = f'{WORKDIR}/{plat}/{dataset}_fig2d2e_tmpfile/{caller}.{dataset}.jasmine.merged.vcf'
            wgs_supp_dict, wgs_merged_total = get_survivor_supp(wgs_merged_vcf)

            unique_pcrt.append((wgs_supp_dict[1] / wgs_merged_total * 100, caller, plat, 'WGS'))
            wgs_supp4_pcrt[plat][this_dataset_index] += wgs_supp_dict[4] / wgs_merged_total * 100

    df_unique = pd.DataFrame(unique_pcrt, columns=['pcrt', 'caller', 'plat', 'region'])

    plat_legends = [Patch(label='HiFi', color=PLATCOLORS['HiFi']), Patch(label='ONT', color=PLATCOLORS['ONT'])]

    fig, ax = plt.subplots(1, 1, figsize=(3, 4))

    sns.barplot(data=df_unique, x='region', y='pcrt', hue='plat', hue_order=['HiFi', 'ONT'], capsize=.15,
                palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=ax)


    ax.set_ylabel('% of aligner unique SVs', fontsize=13)

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
    fig.savefig(f'{figdir}/fig3b.pdf')
    # plt.show()


def plot_s3c(figdir):

    datasets = ['hifi_18kb', 'ont_30kb']

    svtypes = ['INS', 'DEL', 'INV', 'DUP']
    fig_name = 1
    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        fig, ax = plt.subplots(1, 1, figsize=(5, 4))

        svtypes_pcrt = []
        this_total = {aligner: 0 for aligner in READALIGNERS}
        svtypes_count = {aligner: {svtype: 0 for svtype in svtypes} for aligner in READALIGNERS}

        for i in range(len(READALIGNERS)):
            for j in range(len(READALIGNERS)):
                if j <= i:
                    continue
                for caller in READCALLERS:
                    unique_info_out = f'{WORKDIR}/{plat}/{dataset}_figs3_tmpfile/{caller}/{caller}.{READALIGNERS[i]}-{READALIGNERS[j]}.unique.tsv'
                    df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])

                    for idx, row in df_unique.iterrows():
                        svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['ALIGNER'], \
                                                            row['REGION_TYPE']

                        this_total[aligner] += 1
                        if svtype in svtypes:
                            svtypes_count[aligner][svtype] += 1

        for aligner, region_counts in svtypes_count.items():
            for svtype, counts in region_counts.items():
                svtypes_pcrt.append((aligner, dataset, svtype, counts / this_total[aligner] * 100))

        df_svtypes_pcrt = pd.DataFrame(svtypes_pcrt, columns=['aligner', 'dataset', 'svtype', 'pcrt'])

        sns.barplot(data=df_svtypes_pcrt, x='svtype', y='pcrt', hue='aligner', palette=['#1b9e77', '#d95f02', '#7570b3', '#e7298a'], ax=ax)

        ax.set_ylabel(f'% of read aligner specific SV types ({plat})', fontsize=13)
        ax.set_ylim([0, 80])
        ax.set_yticks(np.linspace(0, 80, 5))
        ax.set_yticklabels([int(ele) for ele in np.linspace(0, 80, 5)], fontsize=12)

        ax.set_xlabel('')
        ax.set_xticks([0, 1, 2, 3])
        ax.set_xticklabels(['INS', 'DEL', 'INV', 'DUP'], fontsize=13)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        fig.tight_layout()
        fig.savefig(f'{figdir}/figs3c-{fig_name}.pdf')
        fig_name += 1

    # plt.show()


def main():

    print('\n==== Preparing data for Extended Data Fig 3 =====')
    print(f'Intermediate file directory: {WORKDIR}/dataset_figs3_tmpfile')

    prepare_data()

    if not os.path.exists(f'{FIGDIR}/FigS3'):
        os.mkdir(f'{FIGDIR}/FigS3')

    print('\n==== Creating Extended Data Fig 3a, 3b and 3c =====')
    plot_s3a(f'{FIGDIR}/FigS3')
    print(f'Figures saved to {FIGDIR}/FigS3/figs3a-1.pdf; {FIGDIR}/FigS2/figs3a-2.pdf')

    plot_s3b(f'{FIGDIR}/FigS3')
    print(f'Figures saved to {FIGDIR}/FigS3/figs3b.pdf')

    plot_s3c(f'{FIGDIR}/FigS3')
    print(f'Figures saved to {FIGDIR}/FigS3/figs3c-1.pdf; {FIGDIR}/FigS3/figs3c-2.pdf')


if __name__ == '__main__':
    main()
