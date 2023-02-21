#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/2/21

'''
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


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

def figure3a_3b(workdir, aligners, datasets):

    wgs_stra_pcrts = []
    extd_stra_pcrts = []

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for col_idx, aligner in enumerate(aligners):
        stra_matched = [[0 for i in range(len(datasets))] for j in range(3)]
        for dataset_idx, dataset in enumerate(datasets):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    for assembler in plat_assemblers[plat]:
                        suppvec_info = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.jasmine.suppinfo.tsv'
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

                        df_extd_matched = pd.read_csv(f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.exTD.tsv', header=[0], sep='\t')
                        df_extd_uniques = pd.read_csv(f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.uniques.info.exTD.tsv', header=[0], sep='\t')

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
    ax.set_xticks(np.arange(len(datasets)))
    ax.set_xticklabels([PLATMAP[ele] for ele in datasets], fontsize=13, rotation=90)

    ax.legend(handles=stra_legends)

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()

    ## Figure 3b
    fig, axes = plt.subplots(2, 4, sharex='col', sharey='row', figsize=(7, 4))

    fig_idx = 0
    for assembler in plat_assemblers['HiFi']:
        axes[0][fig_idx].set_title(assembler, fontsize=13)
        wgs_hifi_pcrt = df_wgs_stra_pcrt[df_wgs_stra_pcrt['plat'] == 'HiFi']
        sns.lineplot(data=wgs_hifi_pcrt[wgs_hifi_pcrt['assembler'] == assembler], y='pcrt', x='dataset', hue='aligner',
                     hue_order=aligners, palette=[ALIGNERCOLOR[ele] for ele in aligners], lw=2, marker='o', ci=None,
                     ax=axes[0][fig_idx])

        extd_hifi_pcrt = df_extd_stra_pcrt[df_extd_stra_pcrt['plat'] == 'HiFi']
        sns.lineplot(data=extd_hifi_pcrt[extd_hifi_pcrt['assembler'] == assembler], y='pcrt', x='dataset',
                     hue='aligner',
                     hue_order=aligners, palette=[ALIGNERCOLOR[ele] for ele in aligners], lw=2, ls='--', marker='o',
                     ci=None, ax=axes[1][fig_idx])

        fig_idx += 1

    for assembler in plat_assemblers['ONT']:
        axes[0][fig_idx].set_title(assembler, fontsize=13)
        wgs_ont_pcrt = df_wgs_stra_pcrt[df_wgs_stra_pcrt['plat'] == 'ONT']
        sns.lineplot(data=wgs_ont_pcrt[wgs_ont_pcrt['assembler'] == assembler], y='pcrt', x='dataset', hue='aligner',
                     hue_order=aligners, palette=[ALIGNERCOLOR[ele] for ele in aligners], lw=2, marker='o', ci=None,
                     ax=axes[0][fig_idx])

        extd_ont_pcrt = df_extd_stra_pcrt[df_extd_stra_pcrt['plat'] == 'ONT']
        sns.lineplot(data=extd_ont_pcrt[extd_ont_pcrt['assembler'] == assembler], y='pcrt', x='dataset', hue='aligner',
                     hue_order=aligners, palette=[ALIGNERCOLOR[ele] for ele in aligners], lw=2, ls='--', marker='o',
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
    plt.show()

def figure3c(workdir, dataset):

    aligners = ['minimap2', 'winnowmap', 'ngmlr', 'lra']

    bpstd_pcrt_info = {'HiFi': [], 'ONT': []}
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    plat = 'HiFi'
    if 'ont' in dataset:
        plat = 'ONT'
    for col_idx, aligner in enumerate(aligners):
        for assembler in plat_assemblers[plat]:
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    ins_bpstd_dict = {'Acc': 0, 'Inacc': 0}
                    del_bpstd_dict = {'Acc': 0, 'Inacc': 0}

                    ins_total = 0
                    del_total = 0

                    matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.exTD.tsv'
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

                    bpstd_pcrt_info[plat].append(('INS', ins_bpstd_dict['Acc'] / ins_total * 100, assembler, aligner, caller, asm_caller))
                    bpstd_pcrt_info[plat].append(('DEL', del_bpstd_dict['Acc'] / del_total * 100, assembler, aligner, caller, asm_caller))

    for plat, bpstd_pcrt in bpstd_pcrt_info.items():
        df_bpstd = pd.DataFrame(bpstd_pcrt, columns=['svtype', 'pcrt', 'assembler', 'aligner', 'caller', 'asmcaller'])
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
            axes[1][i].set_xticks(np.arange(len(CALLERS)))
            axes[1][i].set_xticklabels([TOOLMAP[ele] for ele in CALLERS], rotation=90, fontsize=13)

            axes[1][i].spines['left'].set_visible(False)
            axes[1][i].spines['right'].set_visible(False)
            axes[1][i].legend('', frameon=False)

            axes[1][i].set_ylim(80, 100)
            axes[1][i].set_ylabel('% of BSD-10 DEL', fontsize=13)
            axes[1][i].set_yticks(np.linspace(80, 100, 3))
            axes[1][i].set_yticklabels([int(ele) for ele in np.linspace(80, 100, 3)], fontsize=12)

        fig.tight_layout()

    plt.show()

def figure3d(workdir, datasets):

    assmbler_dict = {'HiFi': 'hifiasm', 'ONT': 'flye'}

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        read_supp_counts = {'ins': {}, 'del': {}}
        read_totals = {'ins': 0, 'del': 0}

        assm_supp_counts = {'ins': {}, 'del': {}}
        assm_totals = {'ins': 0, 'del': 0}

        for svtype in ['ins', 'del']:
            for line in open(f'{workdir}/{plat}/read_callers_merged/{dataset}_{svtype}_callers_minimap2_merged_suppinfo.txt'):
                supp_vec, counts = line.strip().split('\t')[0], int(line.strip().split('\t')[1])
                supp = sum([int(val) for val in supp_vec])
                read_totals[svtype] += counts
                if supp in read_supp_counts:
                    read_supp_counts[svtype][supp] += counts
                else:
                    read_supp_counts[svtype][supp] = counts

            for line in open(f'{workdir}/{plat}/assm_callers_merged/{dataset}_{svtype}_callers_{assmbler_dict[plat]}_merged_suppinfo.txt'):
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

    fig.tight_layout()
    plt.show()

def figure3e(workdir):
    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'], 'ONT': ['ont_9kb', 'ont_19kb' , 'ont_30kb']}

    hq_svnums = []

    for region_type in ['WGS', 'ExTD']:
        for plat, datasets in datasets_dict.items():
            for dataset in datasets:
                for svtype in ['ins', 'del']:
                    suppvec_dict = {}
                    for line in open(f'{workdir}/{plat}/stra_compare_hq/{region_type}/{dataset}.{svtype}.suppvec_info.{region_type}.tsv'):
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
    plt.show()

def figure3f(workdir):

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
                for line in open(f'{workdir}/{plat}/stra_compare_hq/{region_type}/{dataset}.{svtype}.suppvec_info.{region_type}.tsv'):
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

    plt.show()