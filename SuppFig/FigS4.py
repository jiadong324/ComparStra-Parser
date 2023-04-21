#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/18

'''
import os
import sys
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Functions import get_overlaps
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


def plot_s4a(figdir):


    dataset = 'hifi_18kb'

    for i, caller in enumerate(READCALLERS):
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        wgs_type_dict = {}

        df_wgs = pd.read_csv(f'{WORKDIR}/HiFi/figs4_tmpfile/{dataset}_{caller}-concordant.info.WGS.tsv', sep='\t')
        for idx, row in df_wgs.iterrows():
            svtype = row['TYPE']
            if svtype in wgs_type_dict:
                wgs_type_dict[svtype] += 1
            else:
                wgs_type_dict[svtype] = 1

        wgs_log_num, wgs_pcrt = [], []
        xticks = []
        incorrect = 0
        for key, val in wgs_type_dict.items():
            if ',' in key:
                incorrect += val
            else:
                wgs_log_num.append(math.log10(val))
                wgs_pcrt.append(val / len(df_wgs))
                xticks.append(key)

        xticks.append('Disc')
        wgs_log_num.append(math.log10(incorrect))
        wgs_pcrt.append(incorrect / len(df_wgs))

        barwidth = 0.5

        bars = ax.bar(np.arange(len(xticks)), wgs_log_num, width=barwidth, color='#1b9e77')
        for j, bar in enumerate(bars):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height + 0.01, f'{round(wgs_pcrt[j] * 100, 2)}%', ha='center', va='bottom')

        bars[-1].set_color('#e7298a')
        ax.set_xticks(np.arange(len(xticks)))
        ax.set_xticklabels(xticks, fontsize=13, rotation=60)
        ax.set_yticks([math.log10(1), math.log10(10), math.log10(100), math.log10(1000), math.log10(10000)])
        ax.set_yticklabels(['0', '10', '100', r'$10^3$', r'$10^4$'], fontsize=12)

        ax.set_ylabel(f'# of SVs from {caller}', fontsize=13)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        fig.tight_layout()
        fig.savefig(f'{figdir}/figs4a-{caller}.pdf')

    # plt.show()

def plot_s4b(figdir):

    aligners = ['minimap2', 'ngmlr', 'lra', 'winnowmap']
    dataset = 'hifi_18kb'

    for i, aligner in enumerate(aligners):
        fig, axes = plt.subplots(1, 2, figsize=(7, 4))
        wgs_type_dict = {}

        df_wgs = pd.read_csv(f'{WORKDIR}/HiFi/figs4_tmpfile/{dataset}_{aligner}-concordant.info.WGS.tsv', sep='\t')
        for idx, row in df_wgs.iterrows():
            svtype = row['TYPE']
            if svtype == 'DEL,INS':
                svtype = 'INS,DEL'
            if svtype == 'DUP,INS':
                svtype = 'INS,DUP'

            if svtype in wgs_type_dict:
                wgs_type_dict[svtype] += 1
            else:
                wgs_type_dict[svtype] = 1

        wgs_log_num, wgs_pcrt = [], []
        xticks = []
        incorrect = 0
        incorrect_dict = {'Rare': 0}
        for key, val in wgs_type_dict.items():
            if ',' in key:
                incorrect += val
                if val > 50:
                    incorrect_dict[key] = val
                else:
                    incorrect_dict['Rare'] += val
            else:
                wgs_log_num.append(math.log10(val))
                wgs_pcrt.append(val / len(df_wgs))
                xticks.append(key)

        xticks.append('Disc')
        wgs_log_num.append(math.log10(incorrect))
        wgs_pcrt.append(incorrect / len(df_wgs))

        barwidth = 0.5

        bars = axes[0].bar(np.arange(len(xticks)), wgs_log_num, width=barwidth, color = '#1b9e77')
        for j, bar in enumerate(bars):
            height = bar.get_height()
            axes[0].text(bar.get_x() + bar.get_width() / 2., height + 0.01, f'{round(wgs_pcrt[j] * 100, 2)}%', ha='center', va='bottom')

        bars[-1].set_color('#e7298a')
        axes[0].set_xticks(np.arange(len(xticks)))
        axes[0].set_xticklabels(xticks, fontsize=13)
        axes[0].set_yticks([math.log10(1), math.log10(10), math.log10(100), math.log10(1000), math.log10(10000)])
        axes[0].set_yticklabels(['0', '10', '100', r'$10^3$', r'$10^4$'], fontsize=12)


        axes[0].set_ylabel(f'NO. of SVs ({aligner})', fontsize=13)

        axes[0].spines['top'].set_visible(False)
        axes[0].spines['right'].set_visible(False)


        sorted_incorrect = sorted(incorrect_dict.items(), key=lambda x:x[1])
        bars = axes[1].barh(np.arange(len(sorted_incorrect)), [math.log10(ele[1]) for ele in sorted_incorrect], height=barwidth, color='#377eb8')
        bars[-1].set_color('#d95f02')
        axes[1].set_yticks(np.arange(len(sorted_incorrect)))
        axes[1].set_yticklabels([ele[0] for ele in sorted_incorrect])
        axes[1].set_xticks([math.log10(1), math.log10(10), math.log10(100), math.log10(1000)])
        axes[1].set_xticklabels(['0', '10', '100', r'$10^3$'], fontsize=12)
        axes[1].set_xlabel('NO. of SVs', fontsize=13)
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['right'].set_visible(False)

        fig.tight_layout()
        fig.savefig(f'{figdir}/figs4b-{aligner}.pdf')

    # plt.show()


def prepare_data():
    dataset = 'hifi_18kb'
    simple_reps = pysam.Tabixfile(SIMREP, 'r')
    rmsk = pysam.Tabixfile(RMSK, 'r')
    sds = pysam.Tabixfile(SD, 'r')

    if not os.path.exists(f'{WORKDIR}/HiFi/figs4_tmpfile'):
        os.mkdir(f'{WORKDIR}/HiFi/figs4_tmpfile')

    for caller in READCALLERS:
        merge_aligner_txt = f'{WORKDIR}/HiFi/figs4_tmpfile/{caller}_{dataset}.filtered_vcf.txt'
        merge_aligner_writer = open(merge_aligner_txt, 'w')
        id_svtype_dict = {}
        for aligner in READALIGNERS:
            vcf_path = f'{WORKDIR}/HiFi/{aligner}_{dataset}/filtered/HG002.{caller}.vcf'
            print(vcf_path, file=merge_aligner_writer)
            this_caller_id_dict = get_sv_id_type_dict(vcf_path)
            id_svtype_dict[aligner] = this_caller_id_dict

        merge_aligner_writer.close()
        merged_vcf = f'{WORKDIR}/HiFi/figs4_tmpfile/{caller}_{dataset}.filtered.alltypes.merged.vcf'
        cmd = f'{JASMINE} file_list={merge_aligner_txt} out_file={merged_vcf} max_dist=1000 pec_len=50 spec_reads=1 --ignore_type > {WORKDIR}/HiFi/figs4_tmpfile/caller.jasmine.log'
        os.system(cmd)
        os.remove(merge_aligner_txt)

        get_aligner_concordant(merged_vcf, id_svtype_dict, f'{WORKDIR}/HiFi/figs4_tmpfile', dataset, caller, READALIGNERS, simple_reps, rmsk, sds)

    for aligner in READALIGNERS:

        merged_txt = f'{WORKDIR}/HiFi/figs4_tmpfile/{aligner}_{dataset}.filtered_vcf.txt'
        merged_txt_writer = open(merged_txt, 'w')

        id_type_dict = {}

        for caller in READCALLERS:
            vcf_path = f'{WORKDIR}/HiFi/{aligner}_{dataset}/filtered/HG002.{caller}.vcf'
            print(vcf_path, file=merged_txt_writer)

            this_caller_id_dict = get_sv_id_type_dict(vcf_path)
            id_type_dict[caller] = this_caller_id_dict

        merged_txt_writer.close()

        merged_vcf = f'{WORKDIR}/HiFi/figs4_tmpfile/{aligner}_{dataset}.filtered.alltypes.merged.vcf'
        cmd = f'{JASMINE} file_list={merged_txt} out_file={merged_vcf} max_dist=1000 pec_len=50 spec_reads=1 --ignore_type > {WORKDIR}/HiFi/figs4_tmpfile/aligner.jasmine.log'
        os.system(cmd)
        os.remove(merged_txt)

        get_caller_concordant(merged_vcf, id_type_dict, f'{WORKDIR}/HiFi/figs4_tmpfile', dataset, aligner, simple_reps,rmsk, sds)

def get_caller_concordant(merged_vcf, idtype_dict, outdir, dataset, aligner, simple_reps, rmsk, sds):

    matched_list = []
    extd_matched_list = []

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            merged_id = entries[2]
            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]
            supp, idlist = int(info_dict['SUPP']), info_dict['IDLIST'].split(',')
            end = int(info_dict['END'])

            if supp == 5:
                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
                svtype_set = set()
                for i, sv_id in enumerate(idlist):
                    svtype_set.add(idtype_dict[READCALLERS[i]][sv_id])
                svtype_str = ','.join(list(svtype_set))
                matched_list.append((entries[0], int(entries[1]), end, info_dict['SVLEN'], merged_id, svtype_str))
                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((entries[0], int(entries[1]), end, info_dict['SVLEN'], merged_id, svtype_str))

    df_matched = pd.DataFrame(matched_list,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'ID_MATCH', 'TYPE'])
    df_matched.to_csv(f'{outdir}/{dataset}_{aligner}-concordant.info.WGS.tsv', header=True, sep='\t', index=False)

    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['#CHROM', 'POS', 'END', 'SVLEN', 'ID_MATCH', 'TYPE'])
    df_extd_matched.to_csv(f'{outdir}/{dataset}_{aligner}-concordant.info.ExTD.tsv', header=True, sep='\t', index=False)



def get_aligner_concordant(merged_vcf, idtype_dict, outdir, dataset, caller, aligners, simple_reps, rmsk, sds):

    matched_list = []
    extd_matched_list = []

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            merged_id = entries[2]
            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]
            supp, idlist = int(info_dict['SUPP']), info_dict['IDLIST'].split(',')
            end = int(info_dict['END'])

            if supp == 4:
                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
                svtype_set = set()
                for i, sv_id in enumerate(idlist):
                    svtype_set.add(idtype_dict[aligners[i]][sv_id])
                svtype_str = ','.join(list(svtype_set))
                matched_list.append((entries[0], int(entries[1]), end, info_dict['SVLEN'], merged_id, svtype_str))
                if region_label != 'Tandem Repeats':
                    extd_matched_list.append((entries[0], int(entries[1]), end, info_dict['SVLEN'], merged_id, svtype_str))

    df_matched = pd.DataFrame(matched_list,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'ID_MATCH', 'TYPE'])
    df_matched.to_csv(f'{outdir}/{dataset}_{caller}-concordant.info.WGS.tsv', header=True, sep='\t', index=False)

    df_extd_matched = pd.DataFrame(extd_matched_list, columns=['#CHROM', 'POS', 'END', 'SVLEN', 'ID_MATCH', 'TYPE'])
    df_extd_matched.to_csv(f'{outdir}/{dataset}_{caller}-concordant.info.ExTD.tsv', header=True, sep='\t', index=False)

def get_sv_id_type_dict(vcf_in):

    id_type_dict = {}

    with open(vcf_in, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split("\t")
            chrom, start, sv_id = entries[0], int(entries[1]), entries[2]
            info_tokens = entries[7].split(";")
            info_dict = {}

            for token in info_tokens:
                if "=" not in token:
                    continue
                info_dict[token.split("=")[0]] = token.split("=")[1].replace(">", "")

            id_type_dict[sv_id] = info_dict['SVTYPE']

    return id_type_dict

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

    sorted_annotations = sorted(annotations, key=lambda x: x[1], reverse=True)

    return sorted_annotations[0]


def main():
    print('\nPrepare data for Extended Data Fig 4 =====')
    print(f'Intermediate file directory: {WORKDIR}/HiFi/figs4_tmpfile')

    prepare_data()

    if not os.path.exists(f'{FIGDIR}/FigS4'):
        os.mkdir(f'{FIGDIR}/FigS4')

    print('\n==== Creating Extended Data Fig 4 =====')


    plot_s4a(f'{FIGDIR}/FigS4')
    print(f'Figures saved to {FIGDIR}/FigS4/figs4a.pdf')

    plot_s4b(f'{FIGDIR}/FigS4')
    print(f'Figures saved to {FIGDIR}/FigS4/figs4b.pdf')


if __name__ == '__main__':
    main()

