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
from matplotlib.patches import Patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Constant import *
from Helpers.Functions import get_overlaps

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


def plot_3d(figdir):

    assmbler_dict = {'HiFi': 'hifiasm', 'ONT': 'flye'}

    legend1 = [Patch(label='HiFi-18kb', facecolor=PLATCOLORS['HiFi'], edgecolor='w'),
                Patch(label='ONT-30kb', facecolor=PLATCOLORS['ONT'], edgecolor='w')]

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    for dataset in ['hifi_18kb', 'ont_30kb']:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        read_supp_counts = {'ins': {}, 'del': {}}
        read_totals = {'ins': 0, 'del': 0}

        assm_supp_counts = {'ins': {}, 'del': {}}
        assm_totals = {'ins': 0, 'del': 0}

        for svtype in ['ins', 'del']:
            for line in open(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/{dataset}_{svtype}_callers_minimap2_merged_suppinfo.txt'):
                supp_vec, counts = line.strip().split('\t')[0], int(line.strip().split('\t')[1])
                supp = sum([int(val) for val in supp_vec])
                read_totals[svtype] += counts
                if supp in read_supp_counts:
                    read_supp_counts[svtype][supp] += counts
                else:
                    read_supp_counts[svtype][supp] = counts

            for line in open(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/{dataset}_{svtype}_callers_{assmbler_dict[plat]}_merged_suppinfo.txt'):
                supp_vec, counts = line.strip().split('\t')[0], int(line.strip().split('\t')[1])
                supp = sum([int(val) for val in supp_vec])
                assm_totals[svtype] += counts
                if supp in assm_supp_counts:
                    assm_supp_counts[svtype][supp] += counts
                else:
                    assm_supp_counts[svtype][supp] = counts


        ax.scatter(read_supp_counts['ins'][5] / read_totals['ins'], read_supp_counts['del'][5] / read_totals['del'], color=PLATCOLORS[plat], edgecolor='k', marker='o', s=220)
        ax.scatter(assm_supp_counts['ins'][2] / assm_totals['ins'], assm_supp_counts['del'][2] / assm_totals['del'], color=PLATCOLORS[plat], edgecolor='k', marker='s', s=220)

        ax.scatter(read_supp_counts['ins'][1] / read_totals['ins'], read_supp_counts['del'][1] / read_totals['del'], color=PLATCOLORS[plat], edgecolor='k', marker='o', hatch='///', s=220)
        ax.scatter(assm_supp_counts['ins'][1] / assm_totals['ins'], assm_supp_counts['del'][1] / assm_totals['del'], color=PLATCOLORS[plat], edgecolor='k', marker='s', hatch='///', s=220)

        ax.set_ylim(-0.05, 0.85)
        ax.set_yticks(np.linspace(0, 0.8, 5))
        ax.set_yticklabels([int(val * 100) for val in np.linspace(0, 0.8, 5)], fontsize=12)
        ax.set_ylabel('% of deletions', fontsize=13)

        ax.set_xlim(-0.05, 0.85)
        ax.set_xticks(np.linspace(0, 0.8, 5))
        ax.set_xticklabels([int(val * 100) for val in np.linspace(0, 0.8, 5)], fontsize=12)
        ax.set_xlabel('% of insertions', fontsize=13)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.plot([-0.05, 0.85], [-0.05, 0.85], color='gray', ls='--', lw=1.5)
        ax.legend(handles=legend1)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig3d.pdf')
    # plt.show()

def plot_3e(figdir):
    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'], 'ONT': ['ont_9kb', 'ont_19kb' , 'ont_30kb']}

    hq_svnums = []

    for region_type in ['WGS', 'ExTD']:
        for plat, datasets in datasets_dict.items():
            for dataset in datasets:
                for svtype in ['ins', 'del']:
                    suppvec_dict = {}
                    for line in open(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/stra_hq_compare/{region_type}/{dataset}.{svtype}.suppvec_info.{region_type}.tsv'):
                        entries = line.strip().split('\t')
                        suppvec, count = entries[0], int(entries[1])
                        suppvec_dict[suppvec] = count
                    hq_svnums.append((region_type, SVTYPEMAP[svtype], 'Read', suppvec_dict['01'] + suppvec_dict['11'], dataset, f'{plat}-{SVTYPEMAP[svtype]}'))
                    hq_svnums.append((region_type, SVTYPEMAP[svtype], 'Assembly', suppvec_dict['10'] + suppvec_dict['11'], dataset, f'{plat}-{SVTYPEMAP[svtype]}'))

    df_svnum = pd.DataFrame(hq_svnums, columns=['region', 'svtype', 'stra', 'count', 'dataset', 'plat'])


    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(4, 6))

    sns.barplot(data=df_svnum[df_svnum['region'] == 'WGS'], x='plat', y='count', hue='stra', capsize=0.15,
                 palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], ax=axes[0])

    axes[0].set_ylabel('# of HC WGS-INS/DEL', fontsize=13)
    axes[0].set_yticks(np.linspace(3000, 12000, 4))
    axes[0].set_yticklabels([int(val) for val in np.linspace(3000, 12000, 4)], fontsize=12)

    sns.barplot(data=df_svnum[df_svnum['region'] == 'ExTD'], x='plat', y='count', hue='stra', capsize=0.15,
                 palette=[STRACOLORS['Read'], STRACOLORS['Assembly']], ax=axes[1])

    axes[1].set_ylabel('# of HC ExTD-INS/DEL', fontsize=13)
    axes[1].set_yticks(np.linspace(0, 3000, 4))
    axes[1].set_yticklabels([int(val) for val in np.linspace(0, 3000, 4)], fontsize=12)

    for i, ax in enumerate(axes):
        ax.legend(title='')
        ax.set_xlabel('')
        # ax.set_xticks([0, 1, 2, 3])
        # ax.set_xticklabels(['HiFi', 'ONT'], fontsize=13)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig3e.pdf')
    # plt.show()

def plot_3f(figdir):

    datasets_dict = {'HiFi': 'hifi_18kb', 'ONT': 'ont_30kb'}

    for region_type in ['WGS', 'ExTD']:
        assembly_unique = {'del': [], 'ins': []}
        overlaps = {'del': [], 'ins': []}
        read_unique = {'del': [], 'ins': []}

        read_pcrt = {'del': [], 'ins': []}
        assembly_pcrt = {'del': [], 'ins': []}

        for plat, dataset in datasets_dict.items():
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

                read_pcrt[svtype].append(suppvec_dict['11'] / (suppvec_dict['01'] + suppvec_dict['11']))
                assembly_pcrt[svtype].append(suppvec_dict['11'] / (suppvec_dict['10'] + suppvec_dict['11']))

        barwidth = 0.3

        legends1 = [Patch(label='Assembly specific', facecolor=STRACOLORS['Assembly'], edgecolor='k', alpha=0.7),
                   Patch(label='Captured by read', facecolor=STRACOLORS['Read'], edgecolor='k', alpha=0.7, hatch='///')]

        legends2 = [Patch(label='Read specific', facecolor=STRACOLORS['Read'], edgecolor='k', alpha=0.7),
                    Patch(label='Captured by assembly', facecolor=STRACOLORS['Assembly'], edgecolor='k', alpha=0.7, hatch='///')]

        r1 = np.arange(2)
        r2 = [r + barwidth + 0.04 for r in r1]
        rs = [r1, r2]
        xticks = [r + (barwidth + 0.04) / 2 for r in r1]

        ## Assembly percentage
        fig, ax = plt.subplots(1, 1, figsize=(5, 3))
        for idx, svtype in enumerate(['del', 'ins']):
            ax.barh(rs[idx], [1 for i in range(2)], height=barwidth, edgecolor='black', color=STRACOLORS['Assembly'], alpha=0.7)
            ax.barh(rs[idx], assembly_pcrt[svtype], height=barwidth, hatch='///', edgecolor='k', facecolor=STRACOLORS['Read'], alpha=0.7)

        ax.legend(handles=legends1, loc='lower left')
        ax.set_xlabel(f'{region_type}-SVs from Assembly', fontsize=13)
        ax.set_yticks(xticks)
        ax.set_yticklabels(['HiFi-18kb', 'ONT-30kb'], fontsize=13)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_xlim([0.5, 1.01])
        ax.set_xticks(np.linspace(0.5, 1, 6))
        ax.set_xticklabels([f'{int(val * 100)}%' for val in np.linspace(0.5, 1, 6)], fontsize=13)

        fig.tight_layout()
        fig.savefig(f'{figdir}/fig3f-{region_type}-1.pdf')

        ## Read percentage
        fig, ax = plt.subplots(1, 1, figsize=(5, 3))
        for idx, svtype in enumerate(['del', 'ins']):
            ax.barh(rs[idx], [1 for i in range(2)], height=barwidth, edgecolor='black', color=STRACOLORS['Read'], alpha=0.7)
            ax.barh(rs[idx], read_pcrt[svtype], height=barwidth, hatch='///', edgecolor='k', facecolor=STRACOLORS['Assembly'], alpha=0.7)

        ax.legend(handles=legends2, loc='lower left')
        ax.set_xlabel(f'{region_type}-SVs from Read', fontsize=13)
        ax.set_yticks(xticks)
        ax.set_yticklabels(['HiFi-18kb', 'ONT-30kb'], fontsize=13)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_xlim([0.5, 1.01])
        ax.set_xticks(np.linspace(0.5, 1, 6))
        ax.set_xticklabels([f'{int(val * 100)}%' for val in np.linspace(0.5, 1, 6)], fontsize=13)

        fig.tight_layout()
        fig.savefig(f'{figdir}/fig3f-{region_type}-2.pdf')

    # plt.show()

def obtain_reads_hq_insdel():
    caller_supp = 5

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for svtype in ['ins', 'del']:
            sr4_caller_vcfs = f'{WORKDIR}/{plat}/{dataset}_callers_sr4_vcf_path.txt'
            sr4_caller_writer = open(sr4_caller_vcfs, 'w')

            for caller in READCALLERS:
                vcf_path = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{caller}.{svtype}.vcf'

                print(vcf_path, file=sr4_caller_writer)

            sr4_caller_writer.close()
            merged_caller_dir = f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile'

            if not os.path.exists(merged_caller_dir):
                os.mkdir(merged_caller_dir)

            caller_sr4_merged_vcf = f'{merged_caller_dir}/{dataset}_{svtype}_callers_minimap2_merged.vcf'
            cmd = f'{JASMINE} file_list={sr4_caller_vcfs} out_file={caller_sr4_merged_vcf} max_dist=1000 --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1 > {merged_caller_dir}/read_hq.jasmine.log'

            os.system(cmd)
            os.remove(sr4_caller_vcfs)
            obtain_read_confident_calls(caller_sr4_merged_vcf, merged_caller_dir, dataset, svtype, caller_supp, simple_reps, rmsk, sds)

def obtain_assm_hq_insdel():

    caller_supp = 2

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    assembler_dict = {'HiFi': 'hifiasm', 'ONT': 'flye'}

    for dataset in DATASETS:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        assembler = assembler_dict[plat]
        for svtype in ['ins', 'del']:

            merged_caller_dir = f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile'

            if not os.path.exists(merged_caller_dir):
                os.mkdir(merged_caller_dir)

            caller_vcf_path = f'{merged_caller_dir}/{dataset}_callers_sa2_vcf_path.txt'
            caller_vcf_path_writer = open(caller_vcf_path, 'w')

            for asmcaller in ASMCALLERS:
                vcf_file = f'{WORKDIR}/{plat}/minimap2_{dataset}/filtered/HG002.{asmcaller}.{assembler}.{svtype}.vcf'
                if asmcaller == 'pav' and plat == 'ONT':
                    vcf_file = f'{WORKDIR}/{plat}/lra_{dataset}/filtered/HG002.{asmcaller}.{assembler}.{svtype}.vcf'
                print(vcf_file,file=caller_vcf_path_writer)

            caller_vcf_path_writer.close()

            print(f' ==== Merge {svtype} on {dataset} ==== ')
            caller_sr2_merged_vcf = f'{merged_caller_dir}/{dataset}_{svtype}_callers_{assembler}_merged.vcf'

            cmd = f'{JASMINE} file_list={caller_vcf_path} out_file={caller_sr2_merged_vcf} max_dist=1000 --normalize_chrs --dup_to_ins ' \
                  f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1 > {merged_caller_dir}/assm_hq.jasmine.log'

            os.system(cmd)
            os.remove(caller_vcf_path)
            obtain_assm_confident_calls(caller_sr2_merged_vcf, merged_caller_dir, dataset, svtype, assembler, caller_supp, simple_reps, rmsk, sds)

def compare_stra_hq_insdel():

    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')
    assembler_dict = {'HiFi': 'hifiasm', 'ONT': 'flye'}

    # dataset_dict = {'ONT': ['ont_9kb', 'ont_19kb', 'ont_30kb']}

    for plat, datasets in DATASET_DICT.items():
        assembler = assembler_dict[plat]
        compare_dir = f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/stra_hq_compare'
        if not os.path.exists(compare_dir):
            os.mkdir(compare_dir)

        for dataset in datasets:
            for svtype in ['ins', 'del']:
                wgs_stra_compar_dir = f'{compare_dir}/WGS'
                if not os.path.exists(wgs_stra_compar_dir):
                    os.mkdir(wgs_stra_compar_dir)

                stra_compare_vcf_path = f'{wgs_stra_compar_dir}/stra_compare_path.txt'
                stra_compare_vcf_path_writer = open(stra_compare_vcf_path, 'w')
                print(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/{dataset}_{svtype}_callers_{assembler}_merged.sc2.vcf', file=stra_compare_vcf_path_writer)
                print(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/{dataset}_{svtype}_callers_minimap2_merged.sc5.vcf', file=stra_compare_vcf_path_writer)

                stra_compare_vcf_path_writer.close()

                stra_compare_merged_vcf = f'{wgs_stra_compar_dir}/{dataset}_{svtype}_assm_read_merged.WGS.vcf'

                cmd = f'{JASMINE} file_list={stra_compare_vcf_path} out_file={stra_compare_merged_vcf} max_dist=1000 --dup_to_ins ' \
                      f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1 > {wgs_stra_compar_dir}/wgs_hq_compare.jasmine.log'
                os.system(cmd)
                os.remove(stra_compare_vcf_path)

                get_hq_compare_info(dataset, svtype, stra_compare_merged_vcf, wgs_stra_compar_dir, 'WGS', simple_reps, rmsk, sds)

                extd_stra_compar_dir = f'{compare_dir}/ExTD'
                if not os.path.exists(extd_stra_compar_dir):
                    os.mkdir(extd_stra_compar_dir)

                extd_stra_compare_vcf_path = f'{extd_stra_compar_dir}/extd_stra_compare_path.txt'
                extd_stra_compare_vcf_path_writer = open(extd_stra_compare_vcf_path, 'w')
                print(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/{dataset}_{svtype}_callers_{assembler}_merged.sc2.extd.vcf',
                      file=extd_stra_compare_vcf_path_writer)
                print(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/{dataset}_{svtype}_callers_minimap2_merged.sc5.extd.vcf',
                      file=extd_stra_compare_vcf_path_writer)

                extd_stra_compare_vcf_path_writer.close()
                extd_stra_compare_merged_vcf = f'{extd_stra_compar_dir}/{dataset}_{svtype}_assm_read_merged.ExTD.vcf'
                cmd = f'{JASMINE} file_list={extd_stra_compare_vcf_path} out_file={extd_stra_compare_merged_vcf} max_dist=1000 --dup_to_ins ' \
                      f'genome_file={HG19REF} samtools_path={SAMTOOLS} spec_len=50 spec_reads=1 > {extd_stra_compar_dir}/extd_hq_compare.jasmine.log'

                os.system(cmd)
                os.remove(extd_stra_compare_vcf_path)

                get_hq_compare_info(dataset, svtype, extd_stra_compare_merged_vcf, extd_stra_compar_dir, 'ExTD', simple_reps, rmsk, sds)


def get_hq_compare_info(dataset, svtype, merged_vcf, compare_outdir, label, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    svs_by_regions = {repeat: [0, 0, 0] for repeat in GENOMICREGIONS}
    supp_vec_dict = {}

    if not os.path.exists(compare_outdir):
        os.mkdir(compare_outdir)

    read_unique_bed = open(f'{compare_outdir}/{dataset}.{svtype}.read-uniques.{label}.bed', 'w')
    assm_unique_bed = open(f'{compare_outdir}/{dataset}.{svtype}.assm-uniques.{label}.bed', 'w')
    assm_unique_vcf = open(f'{compare_outdir}/{dataset}.{svtype}.assm-uniques.{label}.vcf', 'w')
    assm_unique_counter = 0
    read_unique_counter = 0

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=assm_unique_vcf)
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
            sv_length = info_dict['SVLEN']

            if supp_vec in supp_vec_dict:
                supp_vec_dict[supp_vec] += 1
            else:
                supp_vec_dict[supp_vec] = 1


            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if supp == 2:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_list.append((merged_id, merged_type, supp_vec, entries[0], int(entries[1]), int(info_dict['END']), info_dict['SVLEN'], region_label, rptype, pcrt))
                if supp_vec == '01':
                    read_unique_counter += 1
                    svs_by_regions[region_label][2] += 1

                    if svtype == 'ins':
                        print(f'{entries[0]}\t{entries[1]}\t{entries[1]}\t{merged_id}\tINS_{sv_length}', file=read_unique_bed)
                    else:
                        print(f'{entries[0]}\t{entries[1]}\t{end}\t{merged_id}\tDEL', file=read_unique_bed)

                if supp_vec == '10':
                    assm_unique_counter += 1
                    svs_by_regions[region_label][0] += 1
                    print(line.strip(), file=assm_unique_vcf)
                    if svtype == 'ins':
                        print(f'{entries[0]}\t{entries[1]}\t{entries[1]}\t{merged_id}\tINS_{sv_length}\t', file=assm_unique_bed)
                    else:
                        print(f'{entries[0]}\t{entries[1]}\t{end}\t{merged_id}\tDEL', file=assm_unique_bed)

    assm_unique_bed.close()
    read_unique_bed.close()
    assm_unique_vcf.close()

    print(f'Read specific {svtype}: {read_unique_counter}')
    print(f'Assembly specific {svtype}: {assm_unique_counter}')

    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN','START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(f'{compare_outdir}/{dataset}.{svtype}.concordants.info.{label}.tsv', header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', '#CHROM', 'POS', 'END', 'SVLEN','REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(f'{compare_outdir}/{dataset}.{svtype}.uniques.info.{label}.tsv', sep='\t', header=True, index=False)

    suppvec_writer = open(f'{compare_outdir}/{dataset}.{svtype}.suppvec_info.{label}.tsv', 'w')
    for supp_vec, count in supp_vec_dict.items():
        print(f'{supp_vec}\t{count}', file=suppvec_writer)
    suppvec_writer.close()


def obtain_read_confident_calls(merged_vcf, workdir, dataset, svtype, supp_callers, simple_reps, rmsk, sds):
    suppvec_dict = {}
    merged_total = 0
    extd_supp_vec_dict = {}

    suppvec_info = {}
    extd_suppvec_info = {}

    extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged.extd.vcf', 'w')
    scs_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged.sc{supp_callers}.vcf', 'w')
    scs_extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged.sc{supp_callers}.extd.vcf', 'w')

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=extd_merged_vcf)
                print(line.strip(), file=scs_merged_vcf)
                print(line.strip(), file=scs_extd_merged_vcf)
                continue

            entries = line.strip().split('\t')
            info_tokens = entries[7].split(";")
            info_dict = {}
            merged_total += 1

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            supp_vec = info_dict['SUPP_VEC']
            callers = []
            for i, val in enumerate(supp_vec):
                if val == '1':
                    callers.append(TOOLMAP[READCALLERS[i]])

            merged_id = entries[2]

            suppvec_info[merged_id] = ','.join(callers)

            if supp_vec in suppvec_dict:
                suppvec_dict[supp_vec] += 1
            else:
                suppvec_dict[supp_vec] = 1

            merged_type = info_dict['SVTYPE']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if region_label != 'Tandem Repeats':
                print(line.strip(), file=extd_merged_vcf)

                extd_suppvec_info[merged_id] = ','.join(callers)

                if supp_vec in extd_supp_vec_dict:
                    extd_supp_vec_dict[supp_vec] += 1
                else:
                    extd_supp_vec_dict[supp_vec] = 1

            if supp >= supp_callers:
                print(line.strip(), file=scs_merged_vcf)
                if region_label != 'Tandem Repeats':
                    print(line.strip(), file=scs_extd_merged_vcf)

    supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged_suppinfo.txt', 'w')
    for supp, count in suppvec_dict.items():
        print(f'{supp}\t{count}', file=supp_writer)

    supp_writer.close()

    extd_supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_minimap2_merged_suppinfo.extd.txt', 'w')
    for supp, count in extd_supp_vec_dict.items():
        print(f'{supp}\t{count}', file=extd_supp_writer)

    extd_supp_writer.close()

def obtain_assm_confident_calls(merged_vcf, workdir, dataset, svtype, assembler, supp_callers, simple_reps, rmsk, sds):
    suppvec_dict = {}
    merged_total = 0
    extd_supp_vec_dict = {}

    extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged.extd.vcf', 'w')
    scs_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged.sc{supp_callers}.vcf', 'w')
    scs_extd_merged_vcf = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged.sc{supp_callers}.extd.vcf', 'w')

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=extd_merged_vcf)
                print(line.strip(), file=scs_merged_vcf)
                print(line.strip(), file=scs_extd_merged_vcf)
                continue

            entries = line.strip().split('\t')
            info_tokens = entries[7].split(";")
            info_dict = {}
            merged_total += 1

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            supp_vec = info_dict['SUPP_VEC']

            if supp_vec in suppvec_dict:
                suppvec_dict[supp_vec] += 1
            else:
                suppvec_dict[supp_vec] = 1

            merged_type = info_dict['SVTYPE']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

            if region_label != 'Tandem Repeats':
                print(line.strip(), file=extd_merged_vcf)

                if supp_vec in extd_supp_vec_dict:
                    extd_supp_vec_dict[supp_vec] += 1
                else:
                    extd_supp_vec_dict[supp_vec] = 1

            if supp >= supp_callers:
                print(line.strip(), file=scs_merged_vcf)
                if region_label != 'Tandem Repeats':
                    print(line.strip(), file=scs_extd_merged_vcf)

    supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged_suppinfo.txt', 'w')
    for supp, count in suppvec_dict.items():
        print(f'{supp}\t{count}', file=supp_writer)

    supp_writer.close()

    extd_supp_writer = open(f'{workdir}/{dataset}_{svtype}_callers_{assembler}_merged_suppinfo.extd.txt', 'w')
    for supp, count in extd_supp_vec_dict.items():
        print(f'{supp}\t{count}', file=extd_supp_writer)

    extd_supp_writer.close()

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


def main():
    print('\n===== Prepare data for Figure 3d, 3e and 3f =====')
    print(f'Intermediate file directory: {WORKDIR}/platform/fig3d3e3f_tmpfile')

    print(f'\tObtain read based high-confident calls')
    obtain_reads_hq_insdel()
    print(f'\tObtain assembly based high-confident calls')
    obtain_assm_hq_insdel()
    print(f'\tComparing high-confident calls of two strategies')
    compare_stra_hq_insdel()

    print('\n==== Creating Figure 3d, 3e and 3f =====')

    if not os.path.exists(f'{FIGDIR}/Fig3'):
        os.mkdir(f'{FIGDIR}/Fig3')

    plot_3d(f'{FIGDIR}/Fig3')
    print(f'Figures saved to {FIGDIR}/Fig3/fig3d.pdf')

    plot_3e(f'{FIGDIR}/Fig3')
    print(f'Figures saved to {FIGDIR}/Fig3/fig3e.pdf')

    plot_3f(f'{FIGDIR}/Fig3')
    print(f'Figures saved to {FIGDIR}/Fig3/fig3f-WGS-1.pdf; {FIGDIR}/Fig3/fig3f-WGS-2.pdf; {FIGDIR}/Fig3/fig3f-ExTD-1.pdf; {FIGDIR}/Fig3/fig3f-ExTD-2.pdf')


if __name__ == '__main__':
    main()