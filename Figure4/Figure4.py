#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/4/14

'''

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

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

def plot_4a(figdir):

    samples = ['HG003', 'HG004']
    datasets_dict = {'HiFi': 'hifi_18kb', 'ONT': 'ont_30kb'}

    invalid_counter = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}
    valid_total = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}
    uniques_total = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}

    for plat, dataset in datasets_dict.items():
        for svtype_idx, svtype in enumerate(['ins', 'del']):
            trio_scores = []
            trio_valid_qs = {}
            total_valid_svs = 0
            for sample_idx, sample in enumerate(samples):
                with open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/vapor_{svtype}_{sample}.tsv','r') as f:
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

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    width = 0.5
    for idx, plat in enumerate(['HiFi', 'ONT']):
        print('# DEL invalid: {0}/{1} ({2}) {3}'.format(invalid_counter[plat]['del'], valid_total[plat]['del'],
                                                        invalid_counter[plat]['del'] / valid_total[plat][
                                                            'del'] * 100, uniques_total[plat]['del']))
        print('# INS invalid: {0}/{1} ({2}) {3}'.format(invalid_counter[plat]['ins'], valid_total[plat]['ins'],
                                                        invalid_counter[plat]['ins'] / valid_total[plat][
                                                            'ins'] * 100, uniques_total[plat]['ins']))
        ax = axes[idx]
        ax.bar([0, 1], [uniques_total[plat]['del'], uniques_total[plat]['ins']], facecolor='#bfbdbe', edgecolor='w',
               width=width)
        ax.bar([0, 1], [valid_total[plat]['del'], valid_total[plat]['ins']], facecolor='#bfcad6', edgecolor='w',
               width=width)
        ax.bar([0, 1], [invalid_counter[plat]['del'], invalid_counter[plat]['ins']], facecolor='#fbb4ae',
               edgecolor='w', width=width)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['DEL', 'INS'], fontsize=13)

        if idx == 0:
            ax.set_ylabel('Number of read SVs', fontsize=13)

        ax.set_ylim(0, 800)
        ax.set_yticks(np.linspace(0, 800, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 800, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig4a.pdf')
    # plt.show()


def plot_4b(figdir):
    dataset_dict = {'HiFi': 'hifi_18kb', 'ONT': 'ont_30kb'}
    samples = ['HG003', 'HG004']

    invalid_counter = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}
    valid_total = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}
    uniques_total = {'HiFi': {'del': 0, 'ins': 0}, 'ONT': {'del': 0, 'ins': 0}}

    for plat, dataset in dataset_dict.items():
        for svtype in ['del', 'ins']:
            svs_info = {}
            for line in open(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/stra_hq_compare/ExTD/{dataset}_{svtype}_assm_read_merged.ExTD.vcf','r'):
                if '#' in line:
                    continue
                entries = line.strip().split('\t')
                info_tokens = entries[7].split(";")
                info_dict = {}
                for token in info_tokens:
                    if '=' in token:
                        info_dict[token.split('=')[0]] = token.split('=')[1]
                svs_info[f'{entries[0]}-{entries[1]}'] = (entries[2], int(info_dict['SVLEN']), entries[9].split(':')[0])

            df_unique_svs = pd.read_csv(f'{WORKDIR}/{plat}/fig3d3e3f_tmpfile/stra_hq_compare/ExTD/{dataset}.{svtype}.assm-uniques.ExTD.bed', sep='\t',
                names=['chrom', 'start', 'end', 'id', 'svtype'])
            uniques_total[plat][svtype] = len(df_unique_svs)

            sv_valid = {}
            for sample_idx, sample in enumerate(samples):
                with open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/ttmars_combined_res.{svtype}.{sample}.txt','r') as f:
                    for line in f:
                        entries = line.strip().split('\t')
                        valid, chrom, start, end = entries[3], entries[4], entries[5], entries[6]
                        sv_id = f'{chrom}-{start}-{end}'

                        if sv_id not in sv_valid:
                            sv_valid[sv_id] = []
                            sv_valid[sv_id].append(valid)
                        else:
                            sv_valid[sv_id].append(valid)

            for sv_id, values in sv_valid.items():
                if len(values) == 1:
                    continue
                valid_total[plat][svtype] += 1
                if values[0] == 'False' and values[1] == 'False':
                    invalid_counter[plat][svtype] += 1

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    width = 0.5
    for idx, plat in enumerate(['HiFi', 'ONT']):
        print('# DEL invalid: {0}/{1} ({2}) {3}'.format(invalid_counter[plat]['del'], valid_total[plat]['del'],
                                                        invalid_counter[plat]['del'] / valid_total[plat]['del'] * 100,
                                                        uniques_total[plat]['del']))
        print('# INS invalid: {0}/{1} ({2}) {3}'.format(invalid_counter[plat]['ins'], valid_total[plat]['ins'],
                                                        invalid_counter[plat]['ins'] / valid_total[plat]['ins'] * 100,
                                                        uniques_total[plat]['ins']))
        ax = axes[idx]
        ax.bar([0, 1], [uniques_total[plat]['del'], uniques_total[plat]['ins']], facecolor='#bfbdbe', edgecolor='w',
               width=width)
        ax.bar([0, 1], [valid_total[plat]['del'], valid_total[plat]['ins']], facecolor='#bfcad6', edgecolor='w',
               width=width)
        ax.bar([0, 1], [invalid_counter[plat]['del'], invalid_counter[plat]['ins']], facecolor='#fbb4ae', edgecolor='w',
               width=width)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['DEL', 'INS'], fontsize=13)

        if idx == 0:
            ax.set_ylabel('# of SVs', fontsize=13)

        ax.set_ylim(0, 800)
        ax.set_yticks(np.linspace(0, 800, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 800, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig4b.pdf')
    # plt.show()

def plot_4c(figdir):
    dataset_dict = {'HiFi': 'hifi_18kb', 'ONT': 'ont_30kb'}
    columns = ['chrom', 'start', 'end', 'id', 'svlen', 'gt', 'svtype', 'region', 'subtype', 'pcrt']

    invalid_regions = []

    for plat, dataset in dataset_dict.items():
        assm_invalid = pd.read_csv(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/ttmars_trio_invalid_insdel.tsv', sep='\t',
            names=columns)

        region_count = {'Assembly': {reg: {'INS': 0, 'DEL': 0} for reg in GENOMICREGIONS},
                        'Read': {reg: {'INS': 0, 'DEL': 0} for reg in GENOMICREGIONS}}
        invalid_total = {'Assembly': {'INS': 0, 'DEL': 0}, 'Read': {'INS': 0, 'DEL': 0}}

        for i, row in assm_invalid.iterrows():
            svtype, svlen, region = row['svtype'], abs(int(row['svlen'])), row['region']
            region_count['Assembly'][region][svtype] += 1
            invalid_total['Assembly'][svtype] += 1

        for line in open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/vapor_trio_invalid_insdel.tsv', 'r'):
            entries = line.strip().split('\t')
            svtype, svlen, region = entries[3].split('.')[1], abs(int(entries[4])), entries[6]
            region_count['Read'][region][svtype] += 1
            invalid_total['Read'][svtype] += 1

        for stra, regions in region_count.items():
            for reg, count in regions.items():
                if reg == 'Tandem Repeats':
                    continue
                for svtype, val in count.items():
                    invalid_regions.append((stra, svtype, reg, val / invalid_total[stra][svtype] * 100, plat))

    df_invalid_regions = pd.DataFrame(invalid_regions, columns=['stra', 'svtype', 'reg', 'pcrt', 'plat'])
    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    for idx, plat in enumerate(['HiFi', 'ONT']):
        ax = axes[idx]
        # ax.set_title(plat, fontsize=13)
        plat_invalid = df_invalid_regions[df_invalid_regions['plat'] == plat]
        sns.lineplot(data=plat_invalid, x='reg', y='pcrt', hue='stra', style='svtype', ci=None,
                     hue_order=['Assembly', 'Read'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']], lw=2,
                     markers=['o', 'X'], markersize=9, ax=ax)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.set_ylabel('% of invalid SVs', fontsize=13)

        ax.set_xlabel('')
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['RM', 'SD', 'SR'], fontsize=13)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if idx == 1:
            ax.legend('')

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig4c.pdf')
    # plt.show()

def plot_4d(figdir):
    plat = 'HiFi'
    columns = ['chrom', 'start', 'end', 'id', 'svtype', 'svlen', 'region', 'reptype', 'pcrt', 'mapq', 'start_diff_std',
               'size_std', 'size_mean', 'sigs']

    dissize = open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/ttmars_trio_invalid_insdel.discsize.tsv','w')
    dispos = open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/ttmars_trio_invalid_insdel.discpos.tsv', 'w')
    nosigs = open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/ttmars_trio_invalid_insdel.nosigs.tsv', 'w')
    inconclu = open(f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/ttmars_trio_invalid_insdel.inconclu.tsv','w')

    print(f'Chrom\tStart\tEnd\tID\tLength\tSig size Std\tSig pos Std\tSig size mean', file=dissize)
    size_outlier_svinfo = {}
    sig_annot = f'{WORKDIR}/{plat}/ExTD_hq_trio_validation/ttmars_trio_invalid_insdel.annot.tsv'

    df_sig_annot = pd.read_csv(sig_annot, sep='\t', names=columns)
    has_sig_counter = 0
    discsize = 0
    discpos = 0
    others = 0
    nosigs_counter = 0
    for idx, row in df_sig_annot.iterrows():
        chrom, start, end, svtype = row['chrom'], row['start'], row['end'], row['svtype']
        svid, svlen, est_size, size_std, start_std = row['id'], abs(int(row['svlen'])), int(row['size_mean']), int(
            row['size_std']), int(row['start_diff_std'])
        if int(row['mapq']) >= 20 and row['sigs'] >= 1:
            has_sig_counter += 1
            if size_std > 50 or start_std > 50:
                print(f'{chrom}\t{start}\t{end}\t{svid}\t{svlen}\t{size_std}\t{start_std}\t{est_size}', file=dispos)
                discpos += 1
                continue
            if int(row['size_mean']) < 0.7 * svlen or int(row['size_mean']) > (2 - 0.7) * svlen:
                size_outlier_svinfo[row['id']] = [svlen, int(row['size_mean'])]
                print(f'{chrom}\t{start}\t{end}\t{svid}\t{svlen}\t{size_std}\t{start_std}\t{est_size}', file=dissize)
                discsize += 1
                continue
            others += 1
            print(f'{chrom}\t{start}\t{end}\t{svid}\t{svlen}\t{size_std}\t{start_std}\t{est_size}', file=inconclu)
        else:
            print(f'{chrom}\t{start}\t{end}\t{svid}\t{svlen}\t{size_std}\t{start_std}\t{est_size}', file=nosigs)
            nosigs_counter += 1

    dissize.close()
    dispos.close()
    nosigs.close()
    inconclu.close()
    print(has_sig_counter)

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    size = 0.3
    ax.pie([60, 50], radius=1, colors=['#cccccc', '#9ecae1'],
           labels=[60, 50], labeldistance=0.85, wedgeprops=dict(width=size, edgecolor='w'))
    # [np-sigs, csv, inv, fewer-sig, others, novel]
    ax.pie([60, 24, 9, 7, 8, 2], radius=1 - size,
           colors=['#cccccc', '#fdae6b', '#a1d99b', '#bcbddc', '#7fcdbb', '#c994c7'],
           labels=[60, 24, 9, 7, 8, 2], labeldistance=0.7, wedgeprops=dict(width=size, edgecolor='w'))

    fig.tight_layout()
    fig.savefig(f'{figdir}/fig4d.pdf')
    # plt.show()


def main():
    print('\n==== Creating Figure 4 =====')

    if not os.path.exists(f'{FIGDIR}/Fig4'):
        os.mkdir(f'{FIGDIR}/Fig4')

    plot_4a(f'{FIGDIR}/Fig4')
    print(f'Figures saved to {FIGDIR}/Fig4/fig4a.pdf')

    plot_4b(f'{FIGDIR}/Fig4')
    print(f'Figures saved to {FIGDIR}/Fig4/fig4b.pdf')

    plot_4c(f'{FIGDIR}/Fig4')
    print(f'Figures saved to {FIGDIR}/Fig4/fig4c.pdf')

    plot_4d(f'{FIGDIR}/Fig4')
    print(f'Figures saved to {FIGDIR}/Fig4/fig4d.pdf')

if __name__ == '__main__':
    main()
